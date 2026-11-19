#!/usr/bin/env python3
"""Reproduce the signed dual-layer reference-standard selection.

Core: 800 people sampled within ten institution-region by entry-period
strata, using journal-implicit SHA256 ordering and equal-probability
systematic sampling. Special: the frozen high-recall CN/KR/JP candidate set.

The script writes new selection files only; it does not overwrite adjudicated
labels. Final analyses must still read the locked reference-standard sheet.
"""
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parent.parent
CFG = yaml.safe_load((ROOT / "config" / "config.yaml").read_text(encoding="utf-8"))
SEED = "EHJ-R1-reference-v2|20261111"
N_CORE = int(CFG["gender"].get("core_reference_n", 800))
EUROPE = {"AT", "BE", "BG", "CH", "CY", "CZ", "DE", "DK", "EE", "ES",
          "FI", "FR", "GB", "GR", "HR", "HU", "IE", "IS", "IT", "LT",
          "LU", "LV", "MT", "NL", "NO", "PL", "PT", "RO", "RS", "SE", "SI", "SK"}
EAST_ASIA = {"CN", "KR", "JP", "TW", "HK"}


def sha_text(value):
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def region(country):
    country = str(country or "").strip().upper()
    if not country or country in {"未知", "NAN", "NONE"}:
        return "Unknown"
    if country in {"US", "CA"}:
        return "North America"
    if country in EUROPE:
        return "Europe"
    if country in EAST_ASIA:
        return "East Asia"
    return "Other"


def read_cohort(path):
    if path.suffix.lower() == ".parquet":
        return pd.read_parquet(path)
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def main():
    cohort_path = Path(os.environ.get("REFERENCE_COHORT_FILE", ROOT / "data" / "analysis_cohort.csv"))
    cohort = read_cohort(cohort_path)
    required = {"fa_id", "inst_country", "entry_year"}
    missing = required - set(cohort.columns)
    if missing:
        raise ValueError(f"Missing reference-sampling fields: {sorted(missing)}")
    if not cohort.fa_id.is_unique:
        raise ValueError("fa_id must be unique in the reference-sampling frame")
    if "index_journal" not in cohort and "journal" in cohort:
        cohort["index_journal"] = cohort["journal"]
    for col in ["fa_name", "fa_orcid", "institution", "index_journal"]:
        if col not in cohort:
            cohort[col] = ""

    cohort["inst_region"] = cohort.inst_country.map(region)
    cohort["entry_period"] = cohort.entry_year.map(
        lambda value: "2005_2009" if int(value) <= 2009 else "2010")
    cohort["core_stratum"] = cohort.inst_region + "|" + cohort.entry_period

    if "east_asian_name_candidate" in cohort:
        candidate = cohort.east_asian_name_candidate
        if candidate.dtype == object:
            candidate = candidate.astype(str).str.lower().isin({"1", "true", "yes", "y"})
        cohort["east_asian_name_candidate"] = candidate.astype(bool)
    elif "special_selected" in cohort:
        cohort["east_asian_name_candidate"] = pd.to_numeric(cohort.special_selected).eq(1)
    else:
        raise ValueError(
            "Provide east_asian_name_candidate or special_selected from the frozen "
            "high-recall CN/KR/JP candidate screen; final name-origin labels are not a substitute.")

    counts = cohort.groupby("core_stratum", sort=True).size()
    raw = counts * N_CORE / len(cohort)
    allocation = raw.apply(lambda value: int(value // 1))
    remainder = N_CORE - int(allocation.sum())
    if remainder < 0:
        raise ValueError("Core allocation exceeds the target sample size")
    fractional_order = (raw - allocation).sort_values(ascending=False, kind="stable")
    for cell in fractional_order.index[:remainder]:
        allocation.loc[cell] += 1
    if int(allocation.sum()) != N_CORE:
        raise AssertionError("Core allocation does not sum to N_CORE")

    selected_parts = []
    selection_log = []
    for cell in counts.index:
        cell_df = cohort[cohort.core_stratum == cell].copy()
        region_name, entry_period = cell.split("|", 1)
        journal_parts = []
        for journal in ["EHJ", "Circulation", "JACC"]:
            part = cell_df[cell_df.index_journal == journal].copy()
            part["journal_hash"] = part.fa_id.map(
                lambda value: sha_text(f"{SEED}|{cell}|{journal}|{value}"))
            journal_parts.append(part.sort_values(["journal_hash", "fa_id"]))
        ordered = pd.concat(journal_parts, ignore_index=True)
        other = cell_df[~cell_df.index_journal.isin(["EHJ", "Circulation", "JACC"])]
        if len(other):
            raise ValueError(f"Unexpected index journal in {cell}: {sorted(other.index_journal.unique())}")
        N = len(ordered)
        n = int(allocation.loc[cell])
        if N != int(counts.loc[cell]) or not 0 < n <= N:
            raise AssertionError(f"Invalid allocation for {cell}: n={n}, N={N}")
        step = Fraction(N, n)
        start_hex = sha_text(f"{SEED}|{cell}|systematic-start")[:16]
        start = Fraction(int(start_hex, 16), 16**16 - 1) * step
        positions = [int(start + j * step) for j in range(n)]
        if len(set(positions)) != n or min(positions) < 0 or max(positions) >= N:
            raise AssertionError(f"Systematic positions are invalid for {cell}")
        chosen = ordered.iloc[positions].copy()
        chosen["systematic_position_1based"] = [position + 1 for position in positions]
        chosen["systematic_step"] = f"{N}/{n}"
        chosen["systematic_start"] = f"{start_hex} / (16^16 - 1) * ({N}/{n})"
        chosen["inclusion_probability"] = n / N
        chosen["design_weight"] = N / n
        selected_parts.append(chosen)
        selection_log.append({
            "cell": cell, "region": region_name, "entry_period": entry_period,
            "N_cohort": N, "N_sample": n, "raw_quota": float(raw.loc[cell]),
            "start_hash": start_hex, "positions_1based": [position + 1 for position in positions],
            "journal_counts": chosen.index_journal.value_counts().to_dict(),
        })

    core = pd.concat(selected_parts, ignore_index=True)
    if len(core) != N_CORE or not core.fa_id.is_unique:
        raise AssertionError("Core selection is incomplete or duplicated")
    special = cohort[cohort.east_asian_name_candidate].copy()
    core_ids = set(core.fa_id)
    special_ids = set(special.fa_id)
    if "core_selected" in cohort:
        locked_core = set(cohort.loc[pd.to_numeric(cohort.core_selected).eq(1), "fa_id"])
        if locked_core != core_ids:
            raise AssertionError("Generated core selection differs from the locked core_selected indicators")
    if "special_selected" in cohort:
        locked_special = set(cohort.loc[pd.to_numeric(cohort.special_selected).eq(1), "fa_id"])
        if locked_special != special_ids:
            raise AssertionError("Generated special selection differs from the locked special_selected indicators")

    annotation = cohort[cohort.fa_id.isin(core_ids | special_ids)].copy()
    annotation["core_selected"] = annotation.fa_id.isin(core_ids).astype(int)
    annotation["special_selected"] = annotation.fa_id.isin(special_ids).astype(int)
    annotation["sample_role"] = annotation.fa_id.map(
        lambda value: "both" if value in core_ids & special_ids else
        ("core_only" if value in core_ids else "special_only"))
    weight_map = core.set_index("fa_id").design_weight
    probability_map = core.set_index("fa_id").inclusion_probability
    annotation["design_weight"] = annotation.fa_id.map(weight_map)
    annotation["inclusion_probability"] = annotation.fa_id.map(probability_map)
    for col in ["annotator1", "annotator2", "adjudicated", "gender_final"]:
        if col not in annotation:
            annotation[col] = ""

    keep = ["fa_id", "fa_name", "fa_orcid", "institution", "inst_country",
            "index_journal", "entry_year", "entry_period", "inst_region", "core_stratum",
            "core_selected", "special_selected", "sample_role", "inclusion_probability",
            "design_weight", "annotator1", "annotator2", "adjudicated", "gender_final"]
    annotation[keep].sort_values(["sample_role", "core_stratum", "fa_id"]).to_csv(
        ROOT / "data" / "reference_standard_annotation_sheet.csv", index=False)
    core.to_csv(ROOT / "data" / "core_800_selection.csv", index=False)
    special.to_csv(ROOT / "data" / "special_candidates.csv", index=False)
    log = {
        "analyst": "ZZ H", "algorithm": "journal-implicit SHA256 ordering plus systematic sampling",
        "seed": SEED, "source_cohort_file": str(cohort_path), "source_cohort_N": len(cohort),
        "core_N": len(core), "special_N": len(special), "overlap_N": len(core_ids & special_ids),
        "selection_log": selection_log, "new_model_run": False,
    }
    (ROOT / "data" / "reference_standard_selection_log.json").write_text(
        json.dumps(log, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({key: value for key, value in log.items() if key != "selection_log"},
                     ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
