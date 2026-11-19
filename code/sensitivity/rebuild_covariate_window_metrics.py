#!/usr/bin/env python3
"""Rebuild both covariate-window metric sets from local frozen records.

The selection and metric functions preserve the accepted historical window
implementations. The final analysis cohort supplies the locked 1,783-person
ID set; the source run may contain one additional record excluded later.
Only aggregate run summaries are printed. Metric files and optional selection
audits remain in the local data directory.
"""
from __future__ import annotations

from collections import Counter, defaultdict
import gzip
import hashlib
import json
import os
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "code" / "cohort_derivation"))
from entry_anchored import (  # noqa: E402
    augment_author_id_verification,
    boolean,
    date_interval,
    integer,
    iso_date,
    preentry_window_start,
    read_table,
    string_list,
    text,
    validate_preentry_facts,
)

METRIC_FIELDS = [
    "preentry_pubs",
    "preentry_citations",
    "preentry_hindex",
    "preentry_coauthors",
    "intl_collab_share",
]
MODES = ["index_inclusive", "fixed_2000_2004"]


def sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def required_path(environment_name: str) -> Path:
    value = os.environ.get(environment_name)
    if not value:
        raise ValueError(f"Set {environment_name} to the corresponding local input file")
    return Path(value).resolve()


def select_variant(entry_row: pd.Series, facts: pd.DataFrame, mode: str):
    """Select signed facts using the accepted index-inclusive/fixed rules."""
    entry = text(entry_row.entry_id)
    entry_date = iso_date(entry_row.entry_date, f"{entry}.entry_date")
    window_start = preentry_window_start(entry_date)
    entry_work = text(entry_row.entry_work_id)
    subset = facts.loc[facts.entry_id.astype(str) == entry]
    accepted, audit = [], []
    for _, row in subset.iterrows():
        work = text(row.work_id)
        record = row.to_dict()
        reason = ""
        if text(row.verification_status) != "verified":
            reason = "UNVERIFIED_OR_CONFLICTING"
        elif not text(row.publication_date):
            reason = "EXACT_DATE_REQUIRED"
        else:
            date_start, date_end = date_interval(
                row.publication_date, row.date_precision, "publication_date"
            )
            if mode == "index_inclusive":
                if work == entry_work or boolean(row.is_entry_work, "is_entry_work"):
                    if date_start > entry_date:
                        reason = "ON_OR_AFTER_ENTRY_DATE"
                elif date_end < window_start:
                    reason = "BEFORE_ROLLING_5Y_WINDOW"
                elif date_start >= entry_date:
                    reason = "ON_OR_AFTER_ENTRY_DATE"
                elif not (date_start >= window_start and date_end < entry_date):
                    reason = "BOUNDARY_UNCERTAIN"
            elif mode == "fixed_2000_2004":
                if boolean(row.is_entry_work, "is_entry_work") or work == entry_work:
                    reason = "ENTRY_WORK_EXCLUDED"
                elif (
                    date_start < pd.Timestamp("2000-01-01")
                    or date_end >= pd.Timestamp("2005-01-01")
                ):
                    reason = "OUTSIDE_FIXED_2000_2004"
            else:
                raise ValueError(mode)
        if reason:
            record["selection_status"] = "EXCLUDED"
            record["selection_reason"] = reason
            audit.append(record)
            continue
        record["selection_status"] = "INCLUDED"
        record["selection_reason"] = (
            "INDEX_INCLUSIVE_ROLLING_5Y"
            if mode == "index_inclusive"
            else "FIXED_2000_2004"
        )
        record["source_scope"] = (
            "FROZEN_ARCHIVE"
            if text(row.present_in_frozen_profile) == "yes"
            else "VERIFIED_EXTERNAL_SUPPLEMENT"
        )
        accepted.append(record)
        audit.append(record)
    return accepted, audit


def metrics(entry_row: pd.Series, accepted: list[dict]) -> dict:
    """Compute signed-facts metrics, including completeness-gated fields."""
    entry = text(entry_row.entry_id)
    citations = [integer(w["cited_by_count"], "cited_by_count") for w in accepted]
    sources = sorted({text(w["citation_source"]).lower() for w in accepted})
    if len(sources) > 1:
        raise ValueError(f"{entry}: mixed citation sources {sources}")
    ordered = sorted(citations, reverse=True)
    hindex = sum(1 for rank, cited in enumerate(ordered, 1) if cited >= rank)
    aliases = {text(entry_row.fa_id)} | set(string_list(entry_row.author_id_aliases))
    aliases.discard("")
    warnings = []
    if any(text(w["coauthor_coverage"]) != "complete" for w in accepted):
        coauthor_n = None
        warnings.append("COAUTHOR_COVERAGE_INCOMPLETE")
    else:
        coauthors = set()
        for work in accepted:
            coauthors.update(string_list(work["coauthor_identifiers"]))
            coauthors.discard(text(work["byline_name"]))
        coauthors.difference_update(aliases)
        coauthors.discard("")
        coauthor_n = len(coauthors)
    if any(text(w["country_coverage"]) != "complete" for w in accepted):
        intl_share = None
        warnings.append("COUNTRY_COVERAGE_INCOMPLETE")
    else:
        flags = [int(len(set(string_list(w["author_countries"]))) > 1) for w in accepted]
        intl_share = round(sum(flags) / len(flags), 4) if flags else None
    return {
        "preentry_pubs": len(accepted),
        "preentry_citations": sum(citations),
        "preentry_hindex": hindex,
        "preentry_coauthors": coauthor_n,
        "intl_collab_share": intl_share,
        "citation_sources": ";".join(sources),
        "frozen_archive_works": sum(
            w["source_scope"] == "FROZEN_ARCHIVE" for w in accepted
        ),
        "verified_external_works": sum(
            w["source_scope"] == "VERIFIED_EXTERNAL_SUPPLEMENT" for w in accepted
        ),
        "metric_warnings": ";".join(warnings),
    }


def frozen_work_record(work: dict) -> dict:
    """Create the compact frozen-record representation used by the source run."""
    return {
        "date": work.get("date") or "",
        "year": work.get("year"),
        "cited": work.get("cited"),
        "intl": len(
            {
                country
                for author in work["authorships"]
                for country in (author.get("ctry") or [])
            }
        )
        > 1,
        "coauth": {
            author.get("id") or author.get("name")
            for author in work["authorships"]
        },
    }


def build_merge_map(
    cohort: pd.DataFrame,
    verification: pd.DataFrame,
    entry_merges: pd.DataFrame,
) -> dict:
    """Reconstruct verified author-ID aliases using the source-run procedure."""
    verification = augment_author_id_verification(verification, cohort, entry_merges)
    parent = {}

    def find(value):
        parent.setdefault(value, value)
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = parent[value]
        return value

    for _, group in verification.groupby("pair_key"):
        ids = group.fa_id.tolist()
        for other in ids[1:]:
            parent[find(ids[0])] = find(other)
    groups = defaultdict(set)
    for value in parent:
        groups[find(value)].add(value)
    canon = set(cohort.fa_id)
    merge_map = {}
    for values in groups.values():
        target = next((value for value in values if value in canon), sorted(values)[0])
        for value in values:
            merge_map[value] = target
    return merge_map


def load_person_works(pull_paths: list[Path], canon: set, merge_map: dict):
    """Load outcome first and pre-entry second, retaining first-seen works."""
    person_works = defaultdict(dict)
    for pull in pull_paths:
        with gzip.open(pull, "rt", encoding="utf-8") as handle:
            for line in handle:
                work = json.loads(line)
                hit = [
                    (merge_map.get(author.get("id"), author.get("id")), author["pos"])
                    for author in work["authorships"]
                    if author.get("id")
                    and (author.get("id") in merge_map or author.get("id") in canon)
                ]
                hit = [(value, position) for value, position in hit if value in canon]
                if not hit:
                    continue
                slim = frozen_work_record(work)
                for value, position in hit:
                    slot = person_works[value].setdefault(work["work_id"], [slim, []])
                    slot[1].append(position)
    return person_works


def baseline_metrics(
    entry_row: pd.Series,
    works: dict,
    mode: str,
    merge_map: dict,
    facts: pd.DataFrame,
):
    """Compute frozen/raw-path metrics with the accepted window rules."""
    entry_date = iso_date(entry_row.entry_date, "entry_date")
    window_start = preentry_window_start(entry_date)
    covar_pre = []
    for work_id, (work, positions) in works.items():
        date = work.get("date") or f"{work['year']}-07-01"
        if mode == "index_inclusive" and work_id == text(entry_row.entry_work_id):
            verified_entry_date = facts.loc[
                (facts.entry_id.astype(str) == text(entry_row.entry_id))
                & (facts.work_id.astype(str) == work_id)
                & (facts.verification_status.astype(str) == "verified"),
                "publication_date",
            ]
            if len(verified_entry_date) == 1 and text(verified_entry_date.iloc[0]):
                date = text(verified_entry_date.iloc[0])
        timestamp = iso_date(date, "work.date")
        if mode == "index_inclusive":
            is_entry_work = work_id == text(entry_row.entry_work_id)
            if (is_entry_work and timestamp <= entry_date) or (
                not is_entry_work and window_start <= timestamp < entry_date
            ):
                covar_pre.append((work_id, work))
        elif (
            mode == "fixed_2000_2004"
            and pd.Timestamp("2000-01-01")
            <= timestamp
            < pd.Timestamp("2005-01-01")
        ):
            covar_pre.append((work_id, work))
    if mode == "index_inclusive":
        entry_work_id = text(entry_row.entry_work_id)
        if entry_work_id and not any(work_id == entry_work_id for work_id, _ in covar_pre):
            entry_fact = facts.loc[
                (facts.entry_id.astype(str) == text(entry_row.entry_id))
                & (facts.work_id.astype(str) == entry_work_id)
                & (facts.verification_status.astype(str) == "verified")
            ]
            if len(entry_fact) == 1:
                fact = entry_fact.iloc[0]
                fact_date = iso_date(fact.publication_date, "entry_work.publication_date")
                if fact_date <= entry_date:
                    coauthors = set(string_list(fact.coauthor_identifiers))
                    coauthors.discard(text(fact.byline_name))
                    coauthors.discard("")
                    covar_pre.append(
                        (
                            entry_work_id,
                            {
                                "cited": integer(fact.cited_by_count, "cited_by_count"),
                                "intl": len(set(string_list(fact.author_countries))) > 1,
                                "coauth": coauthors,
                            },
                        )
                    )
            elif len(entry_fact) > 1:
                raise ValueError(f"Duplicate verified entry-work facts: {entry_row.entry_id}")
    citations = [work.get("cited") or 0 for _, work in covar_pre]
    ordered = sorted(citations, reverse=True)
    hindex = sum(1 for rank, cited in enumerate(ordered, 1) if cited >= rank)
    coauth = {x for _, work in covar_pre for x in work.get("coauth", set())}
    coauth = {x for x in coauth if merge_map.get(x, x) != entry_row.fa_id}
    intl_flags = [int(bool(work.get("intl"))) for _, work in covar_pre]
    accepted = [
        {"work_id": work_id, "source_scope": "FROZEN_ARCHIVE"}
        for work_id, _ in covar_pre
    ]
    return {
        "preentry_pubs": len(covar_pre),
        "preentry_citations": sum(citations),
        "preentry_hindex": hindex,
        "preentry_coauthors": len(coauth),
        "intl_collab_share": round(sum(intl_flags) / len(intl_flags), 4)
        if intl_flags
        else None,
        "citation_sources": "openalex",
        "frozen_archive_works": len(covar_pre),
        "verified_external_works": 0,
        "metric_warnings": "",
    }, accepted


def read_entry_map(path: Path) -> dict:
    if path.suffix.lower() == ".xlsx":
        entry_map = pd.read_excel(path, sheet_name="Assignments", dtype=str)
    else:
        entry_map = read_table(path)
    if not {"fa_id", "entry_screening_id"}.issubset(entry_map.columns):
        raise ValueError("Entry-ID map requires fa_id and entry_screening_id")
    if not entry_map.fa_id.is_unique:
        raise ValueError("Entry-ID map fa_id values must be unique")
    return entry_map.set_index("fa_id").entry_screening_id.to_dict()


def validate_final_cohort(analysis: pd.DataFrame) -> pd.Series:
    if len(analysis) != 1783 or not analysis.fa_id.is_unique:
        raise ValueError("The final analysis cohort must contain 1,783 unique researchers")
    primary = analysis.primary_cohort.astype(bool)
    if int(primary.sum()) != 1536:
        raise ValueError("The primary cohort must contain 1,536 researchers")
    for event_type, expected in [(1, 866), (2, 534), (0, 136)]:
        observed = int(analysis.loc[primary, "event_type"].eq(event_type).sum())
        if observed != expected:
            raise ValueError(f"event_type={event_type}: expected {expected}, found {observed}")
    calibration = analysis.loc[
        primary
        & analysis.core_selected.eq(1)
        & analysis.p_female.notna()
        & analysis.gender_gold.isin(["F", "M"])
    ]
    if len(calibration) != 495:
        raise ValueError("The primary calibration set must contain 495 researchers")
    return primary


def main():
    analysis_path = Path(
        os.environ.get("ANALYSIS_COHORT_FILE", ROOT / "data" / "analysis_cohort.csv")
    ).resolve()
    source_run = Path(
        os.environ.get("WINDOW_SOURCE_RUN_DIR", ROOT / "data" / "window_source_run")
    ).resolve()
    formal_cohort_path = Path(
        os.environ.get(
            "FORMAL_COHORT_REVISED_FILE", source_run / "data" / "formal_cohort_revised.csv"
        )
    ).resolve()
    verification_path = Path(
        os.environ.get(
            "AUTHOR_ID_VERIFICATION_SHEET",
            source_run / "data" / "author_id_verification_sheet.csv",
        )
    ).resolve()
    outcome_pull = Path(
        os.environ.get(
            "OUTCOME_PULL_FILE", source_run / "data" / "outcome_pull_2026-09-07.jsonl.gz"
        )
    ).resolve()
    preentry_pull = Path(
        os.environ.get(
            "PREENTRY_PULL_FILE", source_run / "data" / "preentry_pull_2026-09-07.jsonl.gz"
        )
    ).resolve()
    facts_path = required_path("PREENTRY_FACTS")
    merges_path = required_path("ENTRY_AUTHOR_ID_MERGES")
    entry_map_path = required_path("ENTRY_ID_MAP")
    output_dir = Path(
        os.environ.get("WINDOW_METRICS_OUTPUT_DIR", ROOT / "data")
    ).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    analysis = pd.read_csv(analysis_path)
    primary = validate_final_cohort(analysis)
    facts = validate_preentry_facts(read_table(facts_path))
    cohort = read_table(formal_cohort_path)
    required_cohort_columns = {
        "fa_id",
        "entry_id",
        "entry_work_id",
        "entry_date",
        "final_decision",
        "covariate_algorithm",
        "author_id_aliases",
    }
    missing = sorted(required_cohort_columns - set(cohort.columns))
    if missing:
        raise ValueError(f"formal_cohort_revised: missing columns {missing}")
    if not cohort.fa_id.is_unique:
        raise ValueError("formal_cohort_revised: fa_id must be unique")
    if not set(analysis.fa_id).issubset(set(cohort.fa_id)):
        raise ValueError("Every final-cohort fa_id must occur in formal_cohort_revised")
    cohort_by_id = cohort.set_index("fa_id", drop=False)

    verification = pd.read_csv(
        verification_path, dtype=str, keep_default_na=False, encoding="utf-8-sig"
    )
    merge_map = build_merge_map(cohort, verification, read_table(merges_path))
    person_works = load_person_works(
        [outcome_pull, preentry_pull], set(cohort.fa_id), merge_map
    )
    entry_map = read_entry_map(entry_map_path)

    summary = {
        "analyst": "ZZ H",
        "status": "PASS_WINDOW_METRICS_REBUILT",
        "full_cohort_N": len(analysis),
        "primary_cohort_N": int(primary.sum()),
        "calibration_N": 495,
        "source_run_cohort_N": len(cohort),
        "variants": {},
    }
    for mode in MODES:
        rows, audit_rows = [], []
        routes = Counter()
        for fa_id in analysis.fa_id:
            entry_row = cohort_by_id.loc[fa_id].copy()
            if not text(entry_row.entry_id):
                if fa_id not in entry_map:
                    raise ValueError("A required entry-ID mapping is absent")
                entry_row["entry_id"] = entry_map[fa_id]
            signed_route = (
                entry_row.final_decision == "retain"
                and entry_row.covariate_algorithm == "signed_preentry_facts"
            )
            if signed_route:
                routes["signed_preentry_facts"] += 1
                entry_facts = facts.loc[facts.entry_id == entry_row.entry_id]
                accepted, audit = select_variant(entry_row, entry_facts, mode)
                values = metrics(entry_row, accepted)
            else:
                routes["frozen_raw_records"] += 1
                works = person_works.get(fa_id, {})
                values, accepted = baseline_metrics(
                    entry_row, works, mode, merge_map, facts
                )
                audit = []
            rows.append({"fa_id": fa_id, **{field: values[field] for field in METRIC_FIELDS}})
            audit_rows.append(
                {
                    "fa_id": fa_id,
                    "mode": mode,
                    "route": "signed_preentry_facts" if signed_route else "frozen_raw_records",
                    **values,
                    "accepted_work_ids": ";".join(w["work_id"] for w in accepted),
                    "selection_audit_rows": len(audit),
                }
            )

        metric_frame = pd.DataFrame(rows)
        if len(metric_frame) != 1783 or not metric_frame.fa_id.is_unique:
            raise ValueError("Reconstructed metric output must contain 1,783 unique IDs")
        if set(metric_frame.fa_id) != set(analysis.fa_id):
            raise ValueError("Reconstructed metric IDs differ from the final analysis cohort")
        output = output_dir / f"window_metrics_{mode}.csv"
        metric_frame.to_csv(output, index=False, encoding="utf-8-sig")
        audit = pd.DataFrame(audit_rows)
        audit_path = output_dir / f"window_selection_audit_{mode}.csv"
        audit.to_csv(audit_path, index=False, encoding="utf-8-sig")

        metric_by_id = metric_frame.set_index("fa_id")
        primary_ids = analysis.loc[primary, "fa_id"]
        primary_metrics = metric_by_id.loc[primary_ids]
        main_pubs = analysis.set_index("fa_id").preentry_pubs
        delta = metric_frame.set_index("fa_id").preentry_pubs - main_pubs
        if mode == "index_inclusive" and not delta.eq(1).all():
            raise ValueError("Index-inclusive preentry_pubs must increase by exactly one")
        summary["variants"][mode] = {
            "output": str(output),
            "sha256": sha256(output),
            "selection_audit": str(audit_path),
            "routes": dict(routes),
            "primary_no_preentry_work_N": int(
                primary_metrics.preentry_pubs.eq(0).sum()
            ),
            "primary_missing_intl_share_N": int(
                primary_metrics.intl_collab_share.isna().sum()
            ),
            "preentry_pubs_delta_vs_main": {
                str(key): int(value) for key, value in delta.value_counts().sort_index().items()
            },
            "new_model_run": False,
        }

    summary_path = output_dir / "window_metrics_reconstruction.json"
    summary_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
