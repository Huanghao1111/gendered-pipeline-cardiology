#!/usr/bin/env python3
"""Prepare analysis inputs for the two prespecified covariate windows."""
import json
import os
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
METRICS = [
    "preentry_pubs", "preentry_citations", "preentry_hindex",
    "preentry_coauthors", "intl_collab_share",
]
VARIANTS = {
    "index_inclusive": "INDEX_INCLUSIVE_METRICS_FILE",
    "fixed_2000_2004": "FIXED_2000_2004_METRICS_FILE",
}


def same_series(a, b):
    return a.eq(b) | (a.isna() & b.isna())


def main():
    main_path = Path(os.environ.get("ANALYSIS_COHORT_FILE", ROOT / "data/analysis_cohort.csv"))
    main = pd.read_csv(main_path)
    assert len(main) == 1783 and main.fa_id.is_unique
    primary = main.primary_cohort.astype(bool)
    assert int(primary.sum()) == 1536
    assert int(main.loc[primary, "event_type"].eq(1).sum()) == 866
    assert int(main.loc[primary, "event_type"].eq(2).sum()) == 534
    assert int(main.loc[primary, "event_type"].eq(0).sum()) == 136
    assert set(METRICS).issubset(main.columns)
    calibration = main.loc[
        primary & main.core_selected.eq(1) & main.p_female.notna()
        & main.gender_gold.isin(["F", "M"])
    ]
    assert len(calibration) == 495

    summary = {"analyst": "ZZ H", "main_analysis": str(main_path), "variants": {}}
    (ROOT / "data").mkdir(exist_ok=True)
    for variant, env_name in VARIANTS.items():
        source_path = Path(
            os.environ.get(env_name, ROOT / "data" / f"window_metrics_{variant}.csv")
        )
        source = pd.read_csv(source_path)
        assert source.fa_id.is_unique and set(source.fa_id) == set(main.fa_id)
        assert set(METRICS).issubset(source.columns)
        source_metrics = source.set_index("fa_id")[METRICS]

        data = main.copy()
        before = data[METRICS].copy()
        for field in METRICS:
            data[field] = data.fa_id.map(source_metrics[field])
        nonmetrics = [column for column in data.columns if column not in METRICS]
        pd.testing.assert_frame_equal(data[nonmetrics], main[nonmetrics], check_dtype=False)
        if variant == "index_inclusive":
            assert (data.preentry_pubs - before.preentry_pubs).eq(1).all()

        out = ROOT / "data" / f"analysis_cohort_{variant}.csv"
        data.to_csv(out, index=False)
        changed = {
            field: int((~same_series(data[field], before[field])).sum())
            for field in METRICS
        }
        summary["variants"][variant] = {
            "metrics_source": str(source_path),
            "output": str(out),
            "N": len(data),
            "primary_N": int(primary.sum()),
            "calibration_N": len(calibration),
            "metric_changes_vs_main": changed,
            "non_window_fields_identical_to_main": True,
        }

    output = ROOT / "data/covariate_window_inputs.json"
    output.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
