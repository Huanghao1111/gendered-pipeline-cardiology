from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import pandas as pd

SCHEMA_VERSION = "1.0.0"
WINDOW_YEARS = 5
ENTRY_PATTERN = re.compile(r"^ENTRY-\d{4}$")
SHA_PATTERN = re.compile(r"^[0-9a-f]{64}$")
POSITIONS = {"first", "middle", "last"}
PRESENCE = {"yes", "no", "uncertain"}
VERIFICATION = {"verified", "unverified", "conflicting"}
COVERAGE = {"complete", "partial", "unknown"}


def is_graduation_author_count_eligible(n_auth: int) -> bool:
    """Graduation requires a multi-author work with no more than 30 authors."""
    return 2 <= n_auth <= 30


SIGNED_REQUIRED = [
    "entry_id", "fa_id", "fa_name", "entry_work_id", "entry_date", "entry_year",
    "baseline_first_pub_year", "baseline_academic_age",
    "earliest_verified_publication_year", "earliest_evidence_url", "change_reason",
    "final_decision", "final_exclusion_reason", "final_adjudicator",
    "final_adjudication_date", "final_signature", "source_workbook_sha256",
]
PREENTRY_REQUIRED = [
    "entry_id", "work_id", "title", "publication_date", "date_precision", "publication_type",
    "author_position", "present_in_frozen_profile", "byline_name", "is_entry_work",
    "verification_status", "evidence_url", "evidence_access_date",
    "cited_by_count", "citation_source", "citation_snapshot_date",
    "coauthor_identifiers", "coauthor_coverage", "author_countries", "country_coverage",
]
PREENTRY_OPTIONAL = [
    "doi", "pmid", "electronic_date", "issue_date", "journal", "notes",
]


def sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path, dtype=object, keep_default_na=False, encoding="utf-8-sig")
    if suffix == ".xlsx":
        return pd.read_excel(path, dtype=object, keep_default_na=False)
    if suffix == ".jsonl":
        return pd.read_json(path, lines=True, dtype=False, keep_default_date=False)
    if suffix == ".json":
        return pd.DataFrame(json.loads(path.read_text(encoding="utf-8")))
    raise ValueError(f"Unsupported table format: {path}")


def text(value) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    return str(value).strip()


def integer(value, field: str) -> int:
    try:
        if isinstance(value, bool) or text(value) == "":
            raise ValueError
        number = float(value)
        if not number.is_integer():
            raise ValueError
        return int(number)
    except Exception as exc:
        raise ValueError(f"{field} must be an integer; found {value!r}") from exc


def iso_date(value, field: str, allow_blank: bool = False) -> pd.Timestamp | None:
    raw = text(value)
    if not raw:
        if allow_blank:
            return None
        raise ValueError(f"{field} is required in YYYY-MM-DD format")
    if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", raw):
        raise ValueError(f"{field} must use full ISO YYYY-MM-DD; found {raw!r}")
    try:
        return pd.Timestamp(raw)
    except Exception as exc:
        raise ValueError(f"{field} is not a valid date: {raw!r}") from exc


def date_interval(value, precision, field: str) -> tuple[pd.Timestamp, pd.Timestamp]:
    """Return the full possible interval for day/month/year precision."""
    raw = text(value)
    precision = text(precision)
    if precision == "day":
        date = iso_date(raw, field)
        return date, date
    if precision == "month":
        if not re.fullmatch(r"\d{4}-\d{2}", raw):
            raise ValueError(f"{field} must use YYYY-MM when date_precision=month")
        start = pd.Timestamp(raw + "-01")
        return start, start + pd.offsets.MonthEnd(0)
    if precision == "year":
        if not re.fullmatch(r"\d{4}", raw):
            raise ValueError(f"{field} must use YYYY when date_precision=year")
        return pd.Timestamp(raw + "-01-01"), pd.Timestamp(raw + "-12-31")
    raise ValueError(f"{field}: date_precision must be day/month/year")


def preentry_window_start(entry_date: pd.Timestamp) -> pd.Timestamp:
    return entry_date - pd.DateOffset(years=WINDOW_YEARS)


def boolean(value, field: str) -> bool:
    if isinstance(value, bool):
        return value
    raw = text(value).lower()
    if raw in {"true", "yes", "1"}:
        return True
    if raw in {"false", "no", "0"}:
        return False
    raise ValueError(f"{field} must be true/false; found {value!r}")


def url(value, field: str) -> str:
    raw = text(value)
    if not re.match(r"^https?://", raw):
        raise ValueError(f"{field} must be an http(s) URL")
    return raw


def string_list(value) -> list[str]:
    if isinstance(value, (list, tuple, set)):
        return [text(x) for x in value if text(x)]
    raw = text(value)
    if not raw:
        return []
    return [x.strip() for x in re.split(r"[;|]", raw) if x.strip()]


def require_columns(df: pd.DataFrame, required: list[str], label: str) -> None:
    missing = sorted(set(required) - set(df.columns))
    if missing:
        raise ValueError(f"{label}: missing required columns {missing}")


def validate_signed_entry_table(df: pd.DataFrame, test_mode: bool = False) -> pd.DataFrame:
    require_columns(df, SIGNED_REQUIRED, "signed_entry_table")
    result = df.copy()
    if result["entry_id"].duplicated().any():
        raise ValueError("signed_entry_table: entry_id must be unique")
    included_decisions = result["final_decision"].isin(["retain", "uncertain_baseline_unchanged"])
    if result.loc[included_decisions, "fa_id"].duplicated().any():
        raise ValueError("signed_entry_table: fa_id must be unique among included records; excluded duplicate-entry rows may share the retained canonical ID")
    for idx, row in result.iterrows():
        entry = text(row["entry_id"])
        if not ENTRY_PATTERN.fullmatch(entry):
            raise ValueError(f"Invalid entry_id: {entry}")
        iso_date(row["entry_date"], f"{entry}.entry_date")
        iso_date(row["final_adjudication_date"], f"{entry}.final_adjudication_date")
        entry_year = integer(row["entry_year"], f"{entry}.entry_year")
        if not 2005 <= entry_year <= 2010:
            raise ValueError(f"{entry}: entry_year outside 2005–2010")
        baseline_year = integer(row["baseline_first_pub_year"], f"{entry}.baseline_first_pub_year")
        baseline_age = integer(row["baseline_academic_age"], f"{entry}.baseline_academic_age")
        if baseline_age != entry_year - baseline_year:
            raise ValueError(f"{entry}: baseline_academic_age does not match frozen values")
        earliest = integer(row["earliest_verified_publication_year"], f"{entry}.earliest_verified_publication_year")
        new_age = entry_year - earliest
        if new_age < 0:
            raise ValueError(f"{entry}: earliest verified year is after entry_year")
        decision = text(row["final_decision"])
        if decision not in {"retain", "exclude", "uncertain_baseline_unchanged"}:
            raise ValueError(f"{entry}: final_decision must be retain, exclude, or uncertain_baseline_unchanged, not {decision!r}")
        reason = text(row["final_exclusion_reason"])
        if decision == "uncertain_baseline_unchanged":
            if earliest != baseline_year or new_age != baseline_age:
                raise ValueError(f"{entry}: uncertain_baseline_unchanged must preserve the frozen year and academic age")
            result.at[idx, "computed_academic_age"] = baseline_age
            continue
        if decision == "exclude" and not reason:
            raise ValueError(f"{entry}: exclusion reason is required")
        if decision == "retain" and new_age > 5:
            raise ValueError(f"{entry}: retain is inconsistent with academic age >5")
        if decision == "exclude" and new_age > 5 and reason != "academic_age_gt_5":
            raise ValueError(f"{entry}: exclusion reason must be academic_age_gt_5")
        if decision == "exclude" and new_age <= 5 and reason == "academic_age_gt_5":
            raise ValueError(f"{entry}: academic-age exclusion is inconsistent with age <=5")
        adjudicator = text(row["final_adjudicator"])
        if not test_mode and adjudicator not in {"Yong Liu", "刘勇"}:
            raise ValueError(f"{entry}: final adjudicator must be Yong Liu / 刘勇")
        if not text(row["final_signature"]):
            raise ValueError(f"{entry}: final_signature is required")
        url(row["earliest_evidence_url"], f"{entry}.earliest_evidence_url")
        if not text(row["change_reason"]):
            raise ValueError(f"{entry}: change_reason is required")
        if not SHA_PATTERN.fullmatch(text(row["source_workbook_sha256"]).lower()):
            raise ValueError(f"{entry}: invalid source_workbook_sha256")
        if not text(row["entry_work_id"]):
            raise ValueError(f"{entry}: entry_work_id is required")
        result.at[idx, "computed_academic_age"] = new_age
    return result


def academic_age_audit(df: pd.DataFrame) -> pd.DataFrame:
    records = []
    for _, row in df.iterrows():
        old_year = integer(row["baseline_first_pub_year"], "baseline_first_pub_year")
        new_year = integer(row["earliest_verified_publication_year"], "earliest_verified_publication_year")
        entry_year = integer(row["entry_year"], "entry_year")
        records.append({
            "entry_id": text(row["entry_id"]),
            "fa_id": text(row["fa_id"]),
            "entry_year": entry_year,
            "baseline_first_pub_year": old_year,
            "earliest_verified_publication_year": new_year,
            "baseline_academic_age": entry_year - old_year,
            "recalculated_academic_age": entry_year - new_year,
            "first_publication_year_changed": old_year != new_year,
            "change_reason": text(row["change_reason"]),
            "source_entry_id": text(row["entry_id"]),
            "evidence_url": text(row["earliest_evidence_url"]),
            "final_decision": text(row["final_decision"]),
            "final_exclusion_reason": text(row["final_exclusion_reason"]),
        })
    return pd.DataFrame(records)


def validate_preentry_facts(df: pd.DataFrame) -> pd.DataFrame:
    require_columns(df, PREENTRY_REQUIRED, "preentry_facts")
    result = df.copy()
    duplicates = result.duplicated(subset=["entry_id", "work_id"], keep=False)
    if duplicates.any():
        ids = result.loc[duplicates, ["entry_id", "work_id"]].to_dict("records")
        raise ValueError(f"preentry_facts: duplicate (entry_id, work_id) primary keys: {ids[:5]}")
    for idx, row in result.iterrows():
        entry = text(row["entry_id"])
        work = text(row["work_id"])
        label = f"{entry}/{work}"
        if not ENTRY_PATTERN.fullmatch(entry) or not work:
            raise ValueError(f"Invalid primary key: {label}")
        for field in ["title", "publication_type", "byline_name"]:
            if not text(row[field]):
                raise ValueError(f"{label}: {field} is required")
        if text(row["author_position"]) not in POSITIONS:
            raise ValueError(f"{label}: author_position must be first/middle/last")
        if text(row["present_in_frozen_profile"]) not in PRESENCE:
            raise ValueError(f"{label}: invalid present_in_frozen_profile")
        boolean(row["is_entry_work"], f"{label}.is_entry_work")
        status = text(row["verification_status"])
        if status not in VERIFICATION:
            raise ValueError(f"{label}: invalid verification_status")
        if status == "verified" or text(row["publication_date"]):
            date_interval(row["publication_date"], row["date_precision"], f"{label}.publication_date")
        for field in ["electronic_date", "issue_date"]:
            if field in result.columns:
                iso_date(row[field], f"{label}.{field}", allow_blank=True)
        if status == "verified":
            url(row["evidence_url"], f"{label}.evidence_url")
            iso_date(row["evidence_access_date"], f"{label}.evidence_access_date")
            integer(row["cited_by_count"], f"{label}.cited_by_count")
            if integer(row["cited_by_count"], f"{label}.cited_by_count") < 0:
                raise ValueError(f"{label}: cited_by_count cannot be negative")
            if not text(row["citation_source"]):
                raise ValueError(f"{label}: citation_source is required for verified rows")
            iso_date(row["citation_snapshot_date"], f"{label}.citation_snapshot_date")
        if text(row["coauthor_coverage"]) not in COVERAGE:
            raise ValueError(f"{label}: invalid coauthor_coverage")
        if text(row["country_coverage"]) not in COVERAGE:
            raise ValueError(f"{label}: invalid country_coverage")
    return result


def select_preentry_facts(entry_row: pd.Series, facts: pd.DataFrame) -> tuple[list[dict], list[dict]]:
    entry = text(entry_row["entry_id"])
    entry_date = iso_date(entry_row["entry_date"], f"{entry}.entry_date")
    window_start = preentry_window_start(entry_date)
    entry_work = text(entry_row["entry_work_id"])
    subset = facts.loc[facts["entry_id"].astype(str) == entry]
    accepted, audit = [], []
    for _, row in subset.iterrows():
        work = text(row["work_id"])
        record = row.to_dict()
        record["window_start_inclusive"] = window_start.strftime("%Y-%m-%d")
        record["window_end_exclusive"] = entry_date.strftime("%Y-%m-%d")
        reason = ""
        status = text(row["verification_status"])
        if status != "verified":
            reason = "UNVERIFIED_OR_CONFLICTING"
        elif boolean(row["is_entry_work"], "is_entry_work") or work == entry_work:
            reason = "ENTRY_WORK_EXCLUDED"
        elif not text(row["publication_date"]):
            reason = "EXACT_DATE_REQUIRED"
        else:
            date_start, date_end = date_interval(row["publication_date"], row["date_precision"], "publication_date")
            if date_end < window_start:
                reason = "BEFORE_ROLLING_5Y_WINDOW"
            elif date_start >= entry_date:
                reason = "ON_OR_AFTER_ENTRY_DATE"
            elif not (date_start >= window_start and date_end < entry_date):
                reason = "BOUNDARY_UNCERTAIN"
        if reason:
            record["selection_status"] = "EXCLUDED"
            record["selection_reason"] = reason
            audit.append(record)
            continue
        record["selection_status"] = "INCLUDED"
        record["selection_reason"] = "INSIDE_ROLLING_5Y_WINDOW"
        record["source_scope"] = (
            "FROZEN_ARCHIVE" if text(row["present_in_frozen_profile"]) == "yes"
            else "VERIFIED_EXTERNAL_SUPPLEMENT"
        )
        accepted.append(record)
        audit.append(record)
    return accepted, audit


def compute_preentry_covariates(entry_row: pd.Series, facts: pd.DataFrame) -> tuple[dict, list[dict]]:
    accepted, audit = select_preentry_facts(entry_row, facts)
    entry = text(entry_row["entry_id"])
    citations = [integer(w["cited_by_count"], "cited_by_count") for w in accepted]
    sources = sorted({text(w["citation_source"]).lower() for w in accepted})
    if len(sources) > 1:
        raise ValueError(f"{entry}: mixed citation sources {sources}; harmonise before metric reconstruction")
    ordered = sorted(citations, reverse=True)
    hindex = sum(1 for rank, cited in enumerate(ordered, 1) if cited >= rank)
    aliases = {text(entry_row["fa_id"])} | set(string_list(entry_row.get("author_id_aliases", "")))
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
        flags = [1 if len(set(string_list(w["author_countries"]))) > 1 else 0 for w in accepted]
        intl_share = round(sum(flags) / len(flags), 4) if flags else None
    metrics = {
        "preentry_pubs": len(accepted),
        "preentry_citations": sum(citations),
        "preentry_hindex": hindex,
        "preentry_coauthors": coauthor_n,
        "intl_collab_share": intl_share,
        "citation_sources": ";".join(sources),
        "frozen_archive_works": sum(w["source_scope"] == "FROZEN_ARCHIVE" for w in accepted),
        "verified_external_works": sum(w["source_scope"] == "VERIFIED_EXTERNAL_SUPPLEMENT" for w in accepted),
        "metric_warnings": ";".join(warnings),
    }
    return metrics, audit


def partition_original_for_entry_anchored(entry_row: pd.Series, works: dict):
    """Partition frozen works without changing the original outcome side.

    post and pre_for_outcome reproduce the frozen d<=entry_date logic. covar_pre
    applies the locked rolling 5-year [entry_date-5y, entry_date) window and
    explicitly removes the entry work, without altering the outcome side.
    """
    entry_date = iso_date(entry_row["entry_date"], "entry_date")
    window_start = preentry_window_start(entry_date)
    entry_work = text(entry_row.get("entry_work_id", ""))
    post, pre_for_outcome, covar_pre = [], [], []
    for work_id, (work, positions) in works.items():
        date = work.get("date") or f"{work['year']}-07-01"
        timestamp = iso_date(date, "work.date")
        if timestamp > entry_date:
            post.append((work, positions, date))
        else:
            pre_for_outcome.append((work, positions, date))
        if window_start <= timestamp < entry_date and work_id != entry_work:
            covar_pre.append((work, positions, date))
    return post, pre_for_outcome, covar_pre


def covariates_from_original_pre(covar_pre: list) -> dict:
    citations = [w.get("cited") or 0 for w, _, _ in covar_pre]
    ordered = sorted(citations, reverse=True)
    hindex = sum(1 for rank, cited in enumerate(ordered, 1) if cited >= rank)
    return {
        "preentry_pubs": len(covar_pre),
        "preentry_citations": sum(citations),
        "preentry_hindex": hindex,
        "preentry_coauthors_raw": len({x for w, _, _ in covar_pre for x in w.get("coauth", set())}),
        "intl_collab_share_raw": (
            round(sum(1 if w.get("intl") else 0 for w, _, _ in covar_pre) / len(covar_pre), 4)
            if covar_pre else None
        ),
    }


def normalize_author_id(value) -> str:
    raw = text(value).rstrip("/")
    if re.fullmatch(r"A\d+", raw):
        return "https://openalex.org/" + raw
    if re.fullmatch(r"https://openalex\.org/A\d+", raw):
        return raw
    raise ValueError(f"Invalid OpenAlex author ID: {value!r}")


def prepare_revised_cohort(frozen: pd.DataFrame, signed: pd.DataFrame) -> pd.DataFrame:
    """Apply signed changes to the complete frozen cohort without an inner join.

    The signed table is a change/adjudication table, not a complete cohort. Rows
    absent from it remain at baseline. Signed rows are linked by entry work, and
    retained rows may then adopt a newly verified canonical author ID.
    """
    require_columns(
        frozen,
        ["fa_id", "work_id", "journal", "entry_year", "entry_date", "fa_name",
         "institution", "corresponding", "first_pub_year", "primary_cohort"],
        "frozen_formal_cohort",
    )
    require_columns(signed, SIGNED_REQUIRED + (["author_id_aliases"] if "author_id_aliases" in signed.columns else []), "signed_entry_table")
    if frozen["work_id"].duplicated().any():
        raise ValueError("frozen_formal_cohort: work_id must be unique")
    if frozen["fa_id"].duplicated().any():
        raise ValueError("frozen_formal_cohort: fa_id must be unique")
    if signed["entry_work_id"].duplicated().any():
        raise ValueError("signed_entry_table: entry_work_id must be unique")

    signed_by_work = {text(row["entry_work_id"]): row for _, row in signed.iterrows()}
    used_signed_works = set()
    records = []
    for _, frozen_row in frozen.iterrows():
        row = frozen_row.copy()
        work_id = text(row["work_id"])
        frozen_fa = normalize_author_id(row["fa_id"])
        first_year = integer(row["first_pub_year"], "first_pub_year")
        entry_year = integer(row["entry_year"], "entry_year")
        row["entry_work_id"] = work_id
        row["_frozen_fa_id"] = frozen_fa
        row["baseline_first_pub_year"] = first_year
        row["baseline_academic_age"] = entry_year - first_year
        signed_row = signed_by_work.get(work_id)
        if signed_row is None:
            row["entry_id"] = ""
            row["fa_id"] = frozen_fa
            row["earliest_verified_publication_year"] = first_year
            row["computed_academic_age"] = entry_year - first_year
            row["final_decision"] = "baseline_unchanged"
            row["final_exclusion_reason"] = ""
            row["change_reason"] = "No signed change row; preserve the frozen baseline record."
            row["author_id_aliases"] = frozen_fa
            row["covariate_algorithm"] = "frozen_entry_anchored_archive"
            records.append(row)
            continue

        used_signed_works.add(work_id)
        entry = text(signed_row["entry_id"])
        if pd.Timestamp(signed_row["entry_date"]) != pd.Timestamp(row["entry_date"]):
            raise ValueError(f"{entry}: signed entry_date differs from the frozen entry work")
        if integer(signed_row["entry_year"], f"{entry}.entry_year") != entry_year:
            raise ValueError(f"{entry}: signed entry_year differs from the frozen entry work")
        decision = text(signed_row["final_decision"])
        covariate_algorithm = text(signed_row.get("covariate_algorithm", "signed_preentry_facts"))
        if covariate_algorithm not in {"signed_preentry_facts", "frozen_entry_anchored_archive"}:
            raise ValueError(f"{entry}: unsupported covariate_algorithm {covariate_algorithm!r}")
        row["covariate_algorithm"] = covariate_algorithm
        row["entry_id"] = entry
        row["final_decision"] = decision
        row["final_exclusion_reason"] = text(signed_row["final_exclusion_reason"])
        row["change_reason"] = text(signed_row["change_reason"])
        row["final_adjudicator"] = text(signed_row["final_adjudicator"])
        row["final_adjudication_date"] = text(signed_row["final_adjudication_date"])
        row["final_signature"] = text(signed_row["final_signature"])
        row["source_workbook_sha256"] = text(signed_row["source_workbook_sha256"])
        if decision == "exclude":
            records.append(row)
            continue
        if decision == "uncertain_baseline_unchanged":
            # The signed canonical ID is documentation only for these profiles.
            # BASELINE_UNCHANGED explicitly means "use only the frozen row".
            row["fa_id"] = frozen_fa
            row["earliest_verified_publication_year"] = first_year
            row["computed_academic_age"] = entry_year - first_year
            row["author_id_aliases"] = frozen_fa
            row["documented_canonical_author_id"] = normalize_author_id(signed_row["fa_id"])
            row["covariate_algorithm"] = "frozen_entry_anchored_archive"
        elif decision == "retain":
            canonical = normalize_author_id(signed_row["fa_id"])
            aliases = {frozen_fa, canonical}
            if "author_id_aliases" in signed_row.index:
                aliases.update(normalize_author_id(value) for value in string_list(signed_row["author_id_aliases"]))
            row["fa_id"] = canonical
            row["earliest_verified_publication_year"] = integer(
                signed_row["earliest_verified_publication_year"], f"{entry}.earliest_verified_publication_year"
            )
            row["computed_academic_age"] = entry_year - int(row["earliest_verified_publication_year"])
            row["author_id_aliases"] = ";".join(sorted(aliases))
        else:
            raise ValueError(f"{entry}: unsupported final_decision {decision!r}")
        records.append(row)

    missing_signed = set(signed_by_work) - used_signed_works
    if missing_signed:
        raise ValueError(f"Signed entry works absent from the frozen cohort: {sorted(missing_signed)[:10]}")
    result = pd.DataFrame(records)
    included = result.loc[result["final_decision"].isin(
        ["baseline_unchanged", "retain", "uncertain_baseline_unchanged"]
    )].copy()
    if included.empty:
        raise ValueError("No cohort rows remain after applying signed decisions")
    if included["fa_id"].duplicated().any():
        duplicates = sorted(included.loc[included["fa_id"].duplicated(keep=False), "fa_id"].unique())
        raise ValueError(f"Unresolved duplicate canonical author IDs remain: {duplicates[:10]}")
    if included["work_id"].duplicated().any():
        raise ValueError("Duplicate entry works remain after applying signed decisions")
    return result.reset_index(drop=True)


def augment_author_id_verification(
    verification: pd.DataFrame,
    revised_cohort: pd.DataFrame,
    entry_merges: pd.DataFrame,
) -> pd.DataFrame:
    """Add signed aliases and entry canonical-ID mappings to the frozen merge sheet."""
    require_columns(verification, ["pair_key", "fa_id"], "author_id_verification_sheet")
    require_columns(
        entry_merges,
        ["entry_screening_id", "old_fa_id", "new_canonical_author_id", "implementation_status"],
        "entry_author_id_merges",
    )
    additions = []
    seen = set(zip(verification["pair_key"].map(text), verification["fa_id"].map(text)))

    def add(pair_key: str, author_id: str):
        key = (text(pair_key), normalize_author_id(author_id))
        if key not in seen:
            additions.append({"pair_key": key[0], "fa_id": key[1]})
            seen.add(key)

    included = revised_cohort.loc[revised_cohort["final_decision"].isin(
        ["baseline_unchanged", "retain", "uncertain_baseline_unchanged"]
    )]
    for _, row in included.iterrows():
        entry = text(row["entry_id"]) or "BASELINE_" + text(row["work_id"]).rsplit("/", 1)[-1]
        pair_key = "R1_SIGNED_" + entry
        add(pair_key, row["fa_id"])
        for alias in string_list(row.get("author_id_aliases", "")):
            add(pair_key, alias)
    for _, row in entry_merges.iterrows():
        if text(row["implementation_status"]) != "READY_FOR_CONTROLLED_INPUT":
            raise ValueError(f"{text(row['entry_screening_id'])}: entry author-ID merge is not ready")
        pair_key = "R1_ENTRY_ALIAS_" + text(row["entry_screening_id"])
        add(pair_key, row["old_fa_id"])
        add(pair_key, row["new_canonical_author_id"])
    if additions:
        extra = pd.DataFrame(additions)
        for column in verification.columns:
            if column not in extra.columns:
                extra[column] = ""
        verification = pd.concat([verification, extra[verification.columns]], ignore_index=True)
    return verification


def require_run_authorization(
    signed_entry_path: str | Path,
    preentry_facts_path: str | Path,
    manifest_path: str | Path,
    test_mode: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    manifest = json.loads(Path(manifest_path).read_text(encoding="utf-8"))
    supported_versions = {SCHEMA_VERSION, "1.1.0-baseline-unchanged"}
    if manifest.get("schema_version") not in supported_versions:
        raise ValueError("Authorization manifest schema_version mismatch")
    if manifest.get("run_authorized") is not True:
        raise ValueError("Formal run is not authorized")
    if manifest.get("schema_version") == "1.1.0-baseline-unchanged" and manifest.get("inputs_complete") is not True:
        raise ValueError("Final implementation inputs are not complete")
    if manifest.get("mode") not in {"baseline_algorithm", "revised"}:
        raise ValueError("mode must be baseline_algorithm or revised")
    expected_window = {
        "years": WINDOW_YEARS,
        "anchor": "entry_date",
        "left_boundary": "inclusive",
        "right_boundary": "exclusive",
        "entry_work": "excluded_by_date_and_entry_work_id",
    }
    if manifest.get("preentry_window") != expected_window:
        raise ValueError(f"preentry_window must exactly match {expected_window}")
    expected_signed = text(manifest.get("signed_entry_table_sha256")).lower()
    expected_facts = text(manifest.get("preentry_facts_sha256")).lower()
    if not SHA_PATTERN.fullmatch(expected_signed) or not SHA_PATTERN.fullmatch(expected_facts):
        raise ValueError("Both input SHA-256 values are required")
    if sha256(signed_entry_path) != expected_signed or sha256(preentry_facts_path) != expected_facts:
        raise ValueError("Input SHA-256 mismatch; refusing to substitute files")
    signed = validate_signed_entry_table(read_table(signed_entry_path), test_mode=test_mode)
    facts = validate_preentry_facts(read_table(preentry_facts_path))
    adjudicators = set(signed["final_adjudicator"].map(text))
    if adjudicators != {text(manifest.get("authorized_by"))}:
        raise ValueError("Authorization signatory does not match the signed entry table")
    signature = text(manifest.get("authorization_signature"))
    if not signature or signature.startswith("PENDING"):
        raise ValueError("A completed authorization_signature is required")
    unknown_entries = set(facts["entry_id"].astype(str)) - set(signed["entry_id"].astype(str))
    if unknown_entries:
        raise ValueError(f"preentry_facts contains unknown ENTRY IDs: {sorted(unknown_entries)[:10]}")
    return signed, facts, manifest
