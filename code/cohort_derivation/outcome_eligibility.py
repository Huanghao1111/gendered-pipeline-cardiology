"""Manual eligibility adjudication for programmatically screened outcomes.

OpenAlex ``article`` records are not automatically excluded by title. Title
patterns only identify records requiring manual review. The locked local
adjudication table supplies final include/exclude decisions by OpenAlex work ID.
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

REVIEW_PATTERN = re.compile(
    r"(?i)\b(editorial|correction|erratum|corrigendum|retraction|reply|"
    r"response\s+to|comment\s+on|letter|review|case\s+report|meta-analysis)\b"
)
REVIEW_COLUMNS = [
    "fa_id", "work_id", "date", "title", "review_required", "final_decision"
]
DECISION_COLUMNS = [
    "work_id", "final_decision", "exclusion_reason", "verification_source",
    "reviewer", "adjudication_status",
]


def needs_manual_review(title: str) -> bool:
    """Flag special or uncertain document types; never auto-exclude by title."""
    title = (title or "").strip()
    return not title or bool(REVIEW_PATTERN.search(title))


def load_adjudications(path: Path) -> dict[str, str]:
    """Read locked include/exclude decisions without distributing person data."""
    if not path.exists():
        return {}
    frame = pd.read_csv(path, encoding="utf-8-sig", dtype=str).fillna("")
    missing = [c for c in ("work_id", "final_decision") if c not in frame.columns]
    if missing:
        raise ValueError(f"Outcome adjudication table lacks: {', '.join(missing)}")
    for column in DECISION_COLUMNS:
        if column not in frame.columns:
            frame[column] = ""
    frame["work_id"] = frame["work_id"].str.strip()
    frame["final_decision"] = frame["final_decision"].str.strip().str.lower()
    if (frame["work_id"] == "").any():
        raise ValueError("Outcome adjudication table contains a blank work_id")
    invalid = sorted(set(frame["final_decision"]) - {"include", "exclude"})
    if invalid:
        raise ValueError(f"final_decision must be include/exclude: {invalid}")
    duplicated = frame[frame.duplicated("work_id", keep=False)]
    if not duplicated.empty:
        raise ValueError("Duplicate work_id in outcome adjudication table")
    excluded = frame.loc[
        (frame["final_decision"] == "exclude")
        & (frame["exclusion_reason"].str.strip() == "")
    ]
    if not excluded.empty:
        raise ValueError("Every excluded outcome requires an exclusion_reason")
    return dict(zip(frame["work_id"], frame["final_decision"]))


def resolve_eligibility(work_id: str, title: str, decisions: dict[str, str]):
    """Return eligible, review-required, and the recorded final decision."""
    required = needs_manual_review(title)
    decision = decisions.get(work_id, "")
    if decision == "exclude":
        return False, required, decision
    if required and decision != "include":
        return False, required, ""
    return True, required, decision


def write_review_log(rows: list[dict], path: Path) -> None:
    """Write a local audit log; the file remains excluded from public release."""
    frame = pd.DataFrame(rows, columns=REVIEW_COLUMNS)
    frame.drop_duplicates(["fa_id", "work_id"], inplace=True)
    frame.sort_values(["fa_id", "date", "work_id"], inplace=True)
    frame.to_csv(path, index=False, encoding="utf-8-sig")
