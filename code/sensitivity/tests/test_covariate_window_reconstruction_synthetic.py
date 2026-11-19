#!/usr/bin/env python3
"""Fully synthetic tests for covariate-window reconstruction and country routing."""
from pathlib import Path
import sys
import unittest

import pandas as pd

SENSITIVITY_DIR = Path(__file__).resolve().parents[1]
ROOT = SENSITIVITY_DIR.parents[1]
sys.path.insert(0, str(SENSITIVITY_DIR))
sys.path.insert(0, str(ROOT / "code" / "cohort_derivation"))

import rebuild_covariate_window_metrics as rebuild  # noqa: E402
from entry_anchored import compute_preentry_covariates  # noqa: E402


ENTRY_ID = "ENTRY-0001"
PERSON = "synthetic-person"
BYLINE = "Synthetic Author"


def fact(
    work_id,
    publication_date,
    citations,
    coauthors,
    countries,
    *,
    is_entry=False,
    country_coverage="complete",
):
    return {
        "entry_id": ENTRY_ID,
        "work_id": work_id,
        "title": f"Synthetic {work_id}",
        "publication_date": publication_date,
        "date_precision": "day",
        "publication_type": "article",
        "author_position": "first",
        "present_in_frozen_profile": "yes",
        "byline_name": BYLINE,
        "is_entry_work": "yes" if is_entry else "no",
        "verification_status": "verified",
        "evidence_url": "https://example.invalid/evidence",
        "evidence_access_date": "2026-01-01",
        "cited_by_count": citations,
        "citation_source": "OpenAlex",
        "citation_snapshot_date": "2026-01-01",
        "coauthor_identifiers": ";".join(coauthors),
        "coauthor_coverage": "complete",
        "author_countries": ";".join(countries),
        "country_coverage": country_coverage,
    }


def entry_row():
    return pd.Series(
        {
            "entry_id": ENTRY_ID,
            "fa_id": PERSON,
            "entry_work_id": "ENTRY-WORK",
            "entry_date": "2005-06-15",
            "author_id_aliases": "synthetic-alias",
            "final_decision": "retain",
            "covariate_algorithm": "signed_preentry_facts",
        }
    )


def signed_facts():
    return pd.DataFrame(
        [
            fact(
                "W1",
                "2001-01-01",
                5,
                [PERSON, "synthetic-alias", "C1"],
                ["US", "CN"],
            ),
            fact("W2", "2004-01-01", 2, [PERSON, "C1", "C2"], ["US"]),
            fact("WSAME", "2005-06-15", 100, [PERSON, "CX"], ["US", "FR"]),
            fact("W4", "2000-01-01", 1, [PERSON, "C1"], ["US"]),
            fact(
                "ENTRY-WORK",
                "2005-06-15",
                4,
                [PERSON, "C3"],
                ["US", "CN"],
                is_entry=True,
            ),
        ]
    )


def raw_work(work_id, date, citations, authorships):
    return {
        "work_id": work_id,
        "date": date,
        "year": int(date[:4]),
        "cited": citations,
        "authorships": authorships,
    }


class CovariateWindowSyntheticTests(unittest.TestCase):
    def test_signed_complete_both_windows_and_same_day_exclusion(self):
        facts = signed_facts()
        row = entry_row()
        index_accepted, _ = rebuild.select_variant(row, facts, "index_inclusive")
        fixed_accepted, _ = rebuild.select_variant(row, facts, "fixed_2000_2004")
        self.assertEqual(
            {work["work_id"] for work in index_accepted}, {"W1", "W2", "ENTRY-WORK"}
        )
        self.assertEqual(
            {work["work_id"] for work in fixed_accepted}, {"W1", "W2", "W4"}
        )
        self.assertNotIn("WSAME", {work["work_id"] for work in index_accepted})

        index_metrics = rebuild.metrics(row, index_accepted)
        self.assertEqual(index_metrics["preentry_pubs"], 3)
        self.assertEqual(index_metrics["preentry_citations"], 11)
        self.assertEqual(index_metrics["preentry_hindex"], 2)
        self.assertEqual(index_metrics["preentry_coauthors"], 3)
        self.assertEqual(index_metrics["intl_collab_share"], 0.6667)

        fixed_metrics = rebuild.metrics(row, fixed_accepted)
        self.assertEqual(fixed_metrics["preentry_pubs"], 3)
        self.assertEqual(fixed_metrics["preentry_citations"], 8)
        self.assertEqual(fixed_metrics["preentry_hindex"], 2)
        self.assertEqual(fixed_metrics["preentry_coauthors"], 2)
        self.assertEqual(fixed_metrics["intl_collab_share"], 0.3333)

    def test_signed_incomplete_country_is_missing_in_both_implementations(self):
        facts = signed_facts()
        facts.loc[facts.work_id.eq("W1"), "country_coverage"] = "partial"
        row = entry_row()
        accepted, _ = rebuild.select_variant(row, facts, "index_inclusive")
        variant_metrics = rebuild.metrics(row, accepted)
        primary_metrics, _ = compute_preentry_covariates(row, facts)
        self.assertIsNone(variant_metrics["intl_collab_share"])
        self.assertIsNone(primary_metrics["intl_collab_share"])
        self.assertIn("COUNTRY_COVERAGE_INCOMPLETE", variant_metrics["metric_warnings"])

    def test_frozen_route_uses_available_countries_without_missingness_gate(self):
        first = raw_work(
            "RAW1",
            "2004-01-01",
            3,
            [
                {"id": PERSON, "name": BYLINE, "pos": "first", "ctry": ["US"]},
                {"id": "X", "name": "X", "pos": "middle", "ctry": []},
                {"id": "Y", "name": "Y", "pos": "last", "ctry": ["CN"]},
            ],
        )
        second = raw_work(
            "RAW2",
            "2003-01-01",
            1,
            [
                {"id": PERSON, "name": BYLINE, "pos": "first", "ctry": ["US"]},
                {"id": "X", "name": "X", "pos": "last", "ctry": []},
            ],
        )
        works = {
            "RAW1": [rebuild.frozen_work_record(first), ["first"]],
            "RAW2": [rebuild.frozen_work_record(second), ["first"]],
        }
        self.assertTrue(works["RAW1"][0]["intl"])
        self.assertFalse(works["RAW2"][0]["intl"])
        facts_without_entry = signed_facts().loc[
            ~signed_facts().work_id.eq("ENTRY-WORK")
        ]
        for mode in rebuild.MODES:
            values, _ = rebuild.baseline_metrics(
                entry_row(), works, mode, {}, facts_without_entry
            )
            self.assertEqual(values["preentry_pubs"], 2)
            self.assertEqual(values["preentry_citations"], 4)
            self.assertEqual(values["preentry_hindex"], 1)
            self.assertEqual(values["preentry_coauthors"], 2)
            self.assertEqual(values["intl_collab_share"], 0.5)

    def test_frozen_index_recovers_entry_fact_but_not_same_day_nonentry(self):
        same_day = raw_work(
            "WSAME",
            "2005-06-15",
            100,
            [{"id": PERSON, "name": BYLINE, "pos": "first", "ctry": ["US"]}],
        )
        works = {"WSAME": [rebuild.frozen_work_record(same_day), ["first"]]}
        entry_fact = signed_facts().loc[
            signed_facts().work_id.eq("ENTRY-WORK")
        ].copy()
        entry_fact.loc[:, "coauthor_identifiers"] = f"{PERSON};C9"
        index_metrics, accepted = rebuild.baseline_metrics(
            entry_row(), works, "index_inclusive", {}, entry_fact
        )
        self.assertEqual([work["work_id"] for work in accepted], ["ENTRY-WORK"])
        self.assertEqual(index_metrics["preentry_pubs"], 1)
        self.assertEqual(index_metrics["preentry_citations"], 4)
        self.assertEqual(index_metrics["preentry_hindex"], 1)
        self.assertEqual(index_metrics["preentry_coauthors"], 1)
        self.assertEqual(index_metrics["intl_collab_share"], 1.0)

        fixed_metrics, fixed_accepted = rebuild.baseline_metrics(
            entry_row(), works, "fixed_2000_2004", {}, entry_fact
        )
        self.assertEqual(fixed_accepted, [])
        self.assertEqual(fixed_metrics["preentry_pubs"], 0)
        self.assertEqual(fixed_metrics["preentry_coauthors"], 0)
        self.assertIsNone(fixed_metrics["intl_collab_share"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
