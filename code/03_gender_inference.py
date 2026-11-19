#!/usr/bin/env python3
"""Estimate name-based probabilities and archive the API responses.

Calibration and multiple imputation are implemented in
code/primary/02_primary_analysis.R. Unresolved names remain in the imputation
frame; no hard classifications are assigned here.
"""
import json, time
from pathlib import Path

import pandas as pd
import requests
import yaml

ROOT = Path(__file__).parent.parent
CFG = yaml.safe_load((ROOT / "config" / "config.yaml").read_text(encoding="utf-8"))
ARCH = ROOT / "data" / "gender_api_archive"
ARCH.mkdir(parents=True, exist_ok=True)

GENDERIZE = "https://api.genderize.io"
# NamSor / Gender-API endpoints require keys; set via environment variables.


def genderize_prob(first_name: str) -> tuple[float | None, dict]:
    r = requests.get(GENDERIZE, params={"name": first_name}, timeout=60)
    j = r.json()
    p = None
    if j.get("gender") and j.get("probability") is not None:
        p = j["probability"] if j["gender"] == "female" else 1 - j["probability"]
    return p, j


def main():
    cohort = pd.read_parquet(ROOT / "data" / "cohort.parquet")
    out = []
    for _, r in cohort.iterrows():
        first = str(r["fa_name"]).split()[0] if r["fa_name"] else ""
        p_f, raw = genderize_prob(first)
        (ARCH / f"{r['fa_id'].split('/')[-1]}.json").write_text(
            json.dumps({"api": "genderize", "query_date": time.strftime("%Y-%m-%d"),
                        "name": first, "response": raw}))
        out.append({"fa_id": r["fa_id"], "p_female_raw": p_f,
                    "source": "genderize"})
        time.sleep(0.1)  # rate discipline
    res = pd.DataFrame(out)
    res.to_parquet(ROOT / "data" / "gender_probs.parquet", index=False)
    print(res["p_female_raw"].describe())


if __name__ == "__main__":
    main()
