#!/usr/bin/env python3
"""Derive A0 cessation, ever-graduation, and re-entry variables."""
from __future__ import annotations

import gzip
import json
import pickle
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from outcome_eligibility import (
    load_adjudications, resolve_eligibility, write_review_log,
)

HERE = Path(__file__).resolve()
OUTROOT = HERE.parents[2]
SOURCE = OUTROOT / "data"
OUT = OUTROOT / "data"
OUT.mkdir(parents=True, exist_ok=True)
CENSOR = pd.Timestamp("2024-12-31")
CONSORTIUM = re.compile(
    r"(?i)\b(consortium|investigators?|collaborat\w*|committee|trialists?|"
    r"task force|working group|study group|coalition|alliance|registry|"
    r"group\b|on behalf of)\b"
)

co = pickle.load(open(SOURCE / "formal_cohort_revised.pkl", "rb"))
verification = pd.read_csv(
    SOURCE / "author_id_verification_sheet_revised.csv", encoding="utf-8-sig", dtype=str
).fillna("")
base = pd.read_csv(SOURCE / "analysis_cohort_base.csv", encoding="utf-8-sig")

# Reproduce the frozen shadow-to-canonical author merge map.
parent: dict[str, str] = {}


def find(x: str) -> str:
    parent.setdefault(x, x)
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


for _, group in verification.groupby("pair_key"):
    ids = group["fa_id"].tolist()
    for other in ids[1:]:
        parent[find(ids[0])] = find(other)

groups = defaultdict(set)
for author_id in parent:
    groups[find(author_id)].add(author_id)
canonical = set(co.fa_id)
merge_map: dict[str, str] = {}
for author_ids in groups.values():
    candidates = [x for x in author_ids if x in canonical]
    target = candidates[0] if candidates else sorted(author_ids)[0]
    for author_id in author_ids:
        merge_map[author_id] = target

person_works = defaultdict(dict)
for pull_name in ("outcome_pull_2026-09-07.jsonl.gz", "preentry_pull_2026-09-07.jsonl.gz"):
    with gzip.open(SOURCE / pull_name, "rt", encoding="utf-8") as handle:
        for line in handle:
            work = json.loads(line)
            hits = [
                (merge_map.get(a.get("id"), a.get("id")), a["pos"])
                for a in work["authorships"]
                if a.get("id") and (a.get("id") in merge_map or a.get("id") in canonical)
            ]
            hits = [(c, p) for c, p in hits if c in canonical]
            if not hits:
                continue
            slim = {
                "work_id": work["work_id"],
                "title": work.get("title") or "",
                "date": work.get("date") or "",
                "year": work.get("year"),
                "type": work.get("type") or "",
                "pmid": work.get("pmid") or "",
                "n_auth": len(work["authorships"]),
                "consortium": any(
                    a.get("name") and CONSORTIUM.search(a["name"])
                    for a in work["authorships"]
                ),
            }
            for canon_id, position in hits:
                slot = person_works[canon_id].setdefault(work["work_id"], [slim, []])
                slot[1].append(position)


def years_elapsed(date: pd.Timestamp, entry: pd.Timestamp) -> float:
    return (date - entry).days / 365.25


def first_gap(years: list[int], entry_year: int, minimum_difference: int):
    if not years:
        return None
    sequence = [entry_year] + [y for y in years if y > entry_year] + [2025]
    gaps = [(a, b) for a, b in zip(sequence, sequence[1:])]
    eligible = [(a, b) for a, b in gaps if b - a >= minimum_difference]
    return eligible[0] if eligible else None


def resolve_first_event(entry: pd.Timestamp, grad_date, cess_date, original=False):
    t_cens = years_elapsed(CENSOR, entry)
    t_grad = years_elapsed(grad_date, entry) if grad_date is not None else float("inf")
    t_cess = years_elapsed(cess_date, entry) if cess_date is not None else float("inf")
    if original and cess_date is not None:
        # Preserve the primary-analysis rounding rule before resolving competing events.
        t_cess = round(t_cess, 3)
    if t_grad <= t_cess and t_grad <= t_cens:
        event_type, time_years = 1, t_grad
    elif t_cess <= t_cens:
        event_type, time_years = 2, t_cess
    else:
        event_type, time_years = 0, t_cens
    return event_type, round(max(time_years, 1 / 365.25), 3)


decisions = load_adjudications(SOURCE / "outcome_eligibility_adjudicated.csv")
review_rows, unresolved = [], set()
rows = []
grad_time_differences = []
for _, record in co.iterrows():
    author_id = record["fa_id"]
    entry = pd.Timestamp(record["entry_date"])
    entry_year = int(record["entry_year"])
    works = person_works.get(author_id, {})
    post, pre = [], []
    for work_id, (work, positions) in works.items():
        date_text = work["date"] or f"{work['year']}-07-01"
        item = (work, positions, date_text, work_id)
        (post if date_text > record["entry_date"] else pre).append(item)

    candidates = []
    for work, positions, date_text, work_id in post:
        if "last" not in positions or work["type"] != "article" or not work["pmid"]:
            continue
        if not 2 <= work["n_auth"] <= 30 or work["consortium"]:
            continue
        eligible, required, decision = resolve_eligibility(work_id, work["title"], decisions)
        if required or decision:
            review_rows.append({"fa_id": author_id, "work_id": work_id, "date": date_text,
                                "title": work["title"], "review_required": required,
                                "final_decision": decision})
        if required and not decision:
            unresolved.add(work_id)
        if not eligible:
            continue
        candidates.append((date_text, work))
    candidates.sort(key=lambda x: x[0])
    grad_date = pd.Timestamp(candidates[0][0]) if candidates else None

    all_years = sorted({w["year"] for w, _, _, _ in post + pre if w["year"]})
    gap_a0 = first_gap(all_years, entry_year, 5)
    cess_a0 = pd.Timestamp(f"{gap_a0[0]}-12-31") if gap_a0 else None
    ev_a0, time_a0 = resolve_first_event(entry, grad_date, cess_a0, original=True)
    grad_observable = grad_date is not None and grad_date <= CENSOR
    ever_type = 1 if grad_observable else 0
    ever_time = round(
        max(years_elapsed(grad_date if grad_observable else CENSOR, entry), 1 / 365.25), 3
    )

    dated_post = sorted(
        {pd.Timestamp(date_text) for _, _, date_text, _ in post if pd.Timestamp(date_text) <= CENSOR}
    )
    last_post = dated_post[-1] if dated_post else entry

    def reentry_fields(event_type, cessation_date, graduation_date):
        if event_type != 2 or cessation_date is None:
            return 0, pd.NaT, 0
        later = [d for d in dated_post if cessation_date < d <= CENSOR]
        reentry_date = later[0] if later else pd.NaT
        after_grad = int(
            bool(later)
            and graduation_date is not None
            and cessation_date < graduation_date <= CENSOR
        )
        return int(bool(later)), reentry_date, after_grad

    re_a0, re_date_a0, re_grad_a0 = reentry_fields(ev_a0, cess_a0, grad_date)

    original = base.loc[base.fa_id == author_id].iloc[0]
    if grad_date is not None:
        grad_time_differences.append(
            abs(years_elapsed(grad_date, entry) - float(original["grad_time_years"]))
        )

    rows.append(
        {
            "fa_id": author_id,
            "event_type_A0": ev_a0,
            "time_years_A0": time_a0,
            "ever_event_A3": ever_type,
            "ever_time_A3": ever_time,
            "grad_date_derived": grad_date.strftime("%Y-%m-%d") if grad_date is not None else "",
            "grad_after_censor": int(grad_date is not None and grad_date > CENSOR),
            "cess_gap_start_year_A0": gap_a0[0] if gap_a0 else pd.NA,
            "cess_gap_next_year_A0": gap_a0[1] if gap_a0 else pd.NA,
            "last_publication_date_through_2024": last_post.strftime("%Y-%m-%d"),
            "last_publication_year_through_2024": int(last_post.year),
            "reentry_A0": re_a0,
            "reentry_date_A0": re_date_a0.strftime("%Y-%m-%d") if pd.notna(re_date_a0) else "",
            "reentry_graduation_A0": re_grad_a0,
        }
    )

write_review_log(review_rows, OUT / "outcome_eligibility_review.csv")
if unresolved:
    raise ValueError(f"{len(unresolved)} title-flagged outcomes require locked include/exclude "
                     "decisions in data/outcome_eligibility_adjudicated.csv")

variants = pd.DataFrame(rows)
merged = base.merge(variants, on="fa_id", how="left", validate="1:1")
assert len(merged) == len(base)
assert base.fa_id.is_unique
assert merged["event_type_A0"].notna().all()

# A0 must reproduce the base person-level event table exactly.
event_mismatch = merged.loc[
    (merged.event_type_A0 != merged.event_type)
    | (abs(merged.time_years_A0 - merged.time_years) > 1e-12),
    ["fa_id", "event_type", "event_type_A0", "time_years", "time_years_A0"],
]
assert event_mismatch.empty, f"A0/base event mismatch: {len(event_mismatch)}"

variants.to_csv(OUT / "event_timelines.csv", index=False, encoding="utf-8-sig")
merged.to_csv(OUT / "analysis_cohort.csv", index=False, encoding="utf-8-sig")
event_mismatch.to_csv(OUT / "A0_event_mismatches.csv", index=False, encoding="utf-8-sig")

primary = merged[merged.primary_cohort == True].copy()  # noqa: E712
summary = {
    "n_full": int(len(merged)),
    "n_primary": int(len(primary)),
    "A0_exact_event_match": bool(event_mismatch.empty),
    "A0_event_mismatch_n": int(len(event_mismatch)),
    "maximum_graduation_time_difference_years": max(grad_time_differences, default=0),
    "graduation_after_2024_n_full": int(merged.grad_after_censor.sum()),
    "event_counts_primary": {
        variant: {str(int(k)): int(v) for k, v in primary[f"event_type_{variant}"].value_counts().sort_index().items()}
        for variant in ("A0",)
    },
    "ever_graduation_primary_n": int(primary.ever_event_A3.sum()),
    "cessation_first_reentry_primary": {
        variant: {
            "cessation_first_n": int((primary[f"event_type_{variant}"] == 2).sum()),
            "reentry_n": int(primary[f"reentry_{variant}"].sum()),
            "reentry_graduation_n": int(primary[f"reentry_graduation_{variant}"].sum()),
        }
        for variant in ("A0",)
    },
    "post2019_last_publication_primary": {
        "n": int((primary.last_publication_year_through_2024 > 2019).sum()),
        "A0_cessation_first_n": int(
            ((primary.last_publication_year_through_2024 > 2019) & (primary.event_type_A0 == 2)).sum()
        ),
    },
    "last_publication_year_among_A0_censored_primary": {
        str(int(k)): int(v)
        for k, v in primary.loc[primary.event_type_A0 == 0, "last_publication_year_through_2024"]
        .value_counts().sort_index().items()
    },
}
with open(OUT / "derivation_audit.json", "w", encoding="utf-8") as handle:
    json.dump(summary, handle, ensure_ascii=False, indent=2)
print(json.dumps(summary, ensure_ascii=False, indent=2))
