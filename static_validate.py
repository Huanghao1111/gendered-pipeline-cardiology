from pathlib import Path
import json
import os
import shutil

HERE = Path(__file__).resolve().parent
ROOT = HERE if (HERE / "README.md").is_file() and (HERE / "code").is_dir() else HERE / "code"
CODE = ROOT / "code"
read = lambda relative: (ROOT / relative).read_text(encoding="utf-8")

for name in [
    "code/primary/02_primary_analysis.R",
    "code/sensitivity/04_index_inclusive_window.R",
    "code/sensitivity/05_fixed_covariate_window.R",
    "code/sensitivity/06_full_cohort.R",
    "code/sensitivity/07_full_cohort_full_calibration.R",
    "code/sensitivity/01_definite_label_ipw_completecase.R",
]:
    text = read(name)
    assert "fail_rate <- sum(br_all$nfail) / (nrow(br_all) * M_BOOT)" in text, name
    assert "bootstrap_failed_reps" in text, name
    assert "stopifnot(!anyDuplicated(parts$b))" in text, name

windows = {
    "code/sensitivity/04_index_inclusive_window.R": (
        "data/analysis_cohort_index_inclusive.csv", "bootstrap_index_inclusive_parts"
    ),
    "code/sensitivity/05_fixed_covariate_window.R": (
        "data/analysis_cohort_fixed_2000_2004.csv", "bootstrap_fixed_2000_2004_parts"
    ),
}
for name, (input_file, cache) in windows.items():
    text = read(name)
    assert f'fread("{input_file}")' in text
    assert "bootstrap07_parts" not in text
    assert text.count(cache) >= 6
    assert "nrow(rs) == 495L" in text
    assert 'Sys.getenv("BOOT_M", 1)' in text

for name in ["code/sensitivity/06_full_cohort.R", "code/sensitivity/07_full_cohort_full_calibration.R"]:
    assert 'Sys.getenv("BOOT_M", 1)' in read(name)
    assert "M_b = 2" not in read(name)
assert "nrow(rs) == 495L" in read("code/primary/02_primary_analysis.R")
assert 'Sys.getenv("COUNTRY_TOPN", "3")' in read("code/sensitivity/01_definite_label_ipw_completecase.R")
assert "SJR-2024 Cardiology and Cardiovascular Medicine Q1" in read("config/config.yaml")
assert "core_reference_n: 800" in read("config/config.yaml")
assert "prepare_covariate_windows.py" in read("README.md")
assert "BOOT_M=1" in read("README.md") and "BOOT_M=2" in read("README.md")
assert "866 graduations, 534 cessation-first events, and 136 administrative censored records" in read("README.md")
assert "01b_primary_calibration_cv.R" in read("README.md")
primary_cv = read("code/primary/01b_primary_calibration_cv.R")
assert "nrow(cal) == 495L" in primary_cv
assert "SEED <- 20261112L" in primary_cv
assert "K <- 10L" in primary_cv and "REPEATS <- 10L" in primary_cv
assert "fit <- fit_calibration(cal[train])" in primary_cv
assert "sum(d0$event_type[primary] == 1) == 866L" in primary_cv
assert "sum(d0$event_type[primary] == 2) == 534L" in primary_cv
assert "sum(d0$event_type[primary] == 0) == 136L" in primary_cv

for path in CODE.rglob("*.R"):
    text = path.read_text(encoding="utf-8")
    assert "unique(rbindlist" not in text, str(path)
    assert "fail_rate <- sum(br$nfail)" not in text, str(path)
    assert "event_type_A0" not in text, str(path)
    assert "time_years_A0" not in text, str(path)

primary_guards = [
    "code/primary/02_primary_analysis.R",
    "code/descriptive/common.R",
    "code/secondary/common.R",
    "code/sensitivity/common.R",
    "code/sensitivity/04_index_inclusive_window.R",
    "code/sensitivity/05_fixed_covariate_window.R",
]
for name in primary_guards:
    text = read(name)
    assert "sum(d0$ev1) == 866L" in text, name
    assert "sum(d0$ev2) == 534L" in text, name
    assert "sum(d0$ev1 == 0L & d0$ev2 == 0L) == 136L" in text, name
for name in [
    "code/sensitivity/01_definite_label_ipw_completecase.R",
    "code/sensitivity/06_full_cohort.R",
    "code/sensitivity/07_full_cohort_full_calibration.R",
]:
    text = read(name)
    assert "d0[primary_cohort == TRUE, sum(ev1)] == 866L" in text, name
    assert "d0[primary_cohort == TRUE, sum(ev2)] == 534L" in text, name
    assert "d0[primary_cohort == TRUE, sum(ev1 == 0L & ev2 == 0L)] == 136L" in text, name
assert "c(136L, 866L, 534L)" in read("code/sensitivity/08_cessation_time_interaction.R")

rebuild = read("code/sensitivity/rebuild_covariate_window_metrics.py")
prepare = read("code/sensitivity/prepare_covariate_windows.py")
readme = read("README.md")
assert "def select_variant" in rebuild and "def baseline_metrics" in rebuild
assert 'MODES = ["index_inclusive", "fixed_2000_2004"]' in rebuild
assert 'source_run / "data" / "formal_cohort_revised.csv"' in rebuild
assert "set(analysis.fa_id).issubset(set(cohort.fa_id))" in rebuild
assert 'f"window_metrics_{mode}.csv"' in rebuild
assert 'ROOT / "data" / f"window_metrics_{variant}.csv"' in prepare
assert "rebuild_covariate_window_metrics.py" in readme
assert "separately reconstructed five-metric files" not in readme
assert "per-record mixed route" in readme
assert "frozen-set reproduction" in readme
assert (ROOT / "code/sensitivity/tests/test_covariate_window_reconstruction_synthetic.py").is_file()

cohort_checked = False
legacy_fields_identical = None
cohort_file = os.environ.get("ANALYSIS_COHORT_FILE")
if cohort_file:
    import numpy as np
    import pandas as pd

    data = pd.read_csv(cohort_file)
    primary = data["primary_cohort"].astype(bool)
    assert len(data) == 1783 and data["fa_id"].is_unique
    assert int(primary.sum()) == 1536
    assert data.loc[primary, "event_type"].value_counts().to_dict() == {1: 866, 2: 534, 0: 136}
    assert data.loc[primary, "time_years"].notna().all()
    assert data.loc[primary, "time_years"].gt(0).all()
    if "event_type_A0" in data and "time_years_A0" in data:
        assert data["event_type"].eq(data["event_type_A0"]).all()
        assert np.allclose(data["time_years"], data["time_years_A0"], rtol=0, atol=1e-12, equal_nan=True)
        legacy_fields_identical = True
    cohort_checked = True

print(json.dumps({
    "status": "PASS_STATIC_CHECKS",
    "R_scripts_executed": False,
    "Rscript_available": shutil.which("Rscript") is not None,
    "targeted_bootstrap_scripts_checked": 6,
    "window_inputs_and_caches_separated": True,
    "duplicate_bootstrap_ids_rejected": True,
    "legacy_R_event_fields_absent": True,
    "final_cohort_checked": cohort_checked,
    "final_primary_events": {"graduation": 866, "cessation_first": 534, "administrative_censoring": 136},
    "legacy_fields_identical_when_present": legacy_fields_identical,
}, indent=2))
