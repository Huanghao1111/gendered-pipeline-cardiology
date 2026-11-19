#!/usr/bin/env python3
"""Clean article metadata and record exclusions using config.yaml.

Only bibliographic metadata available at cohort entry are used.
"""
import re
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).parent.parent
CFG = yaml.safe_load((ROOT / "config" / "config.yaml").read_text(encoding="utf-8"))
CJK = re.compile(r"[぀-ヿ一-鿿가-힯]")


def load_raw() -> pd.DataFrame:
    files = sorted((ROOT / "data").glob("articles_raw_*.parquet"))
    assert files, "Place the frozen articles_raw_*.parquet extract in data/"
    return pd.read_parquet(files[-1])


def add_title(df: pd.DataFrame) -> pd.DataFrame:
    """Reconstruct plain title; also flag CJK-script titles/abstracts."""
    df["title_clean"] = df["title"].fillna("").str.strip()
    df["has_cjk"] = df["title_clean"].str.contains(CJK)
    return df


def title_rule_mask(df: pd.DataFrame) -> pd.Series:
    pats = CFG["article_cleaning"]["title_patterns_exclude"]
    m = pd.Series(False, index=df.index)
    for p in pats:
        m |= df["title_clean"].str.lower().str.contains(p, regex=True)
    return m


def anomalous_issue_mask(df: pd.DataFrame) -> pd.Series:
    """Flag supplement-labelled issues and issues at least factor times
    the journal-year median issue size."""
    factor = CFG["article_cleaning"]["issue_anomaly_factor"]
    rescued = set(CFG["article_cleaning"]["manual_rescue"])
    bib = df["biblio"].apply(lambda b: b if isinstance(b, dict) else {})
    df["_issue"] = bib.apply(lambda b: str(b.get("issue") or ""))
    df["_jy"] = df["_journal"] + ":" + df["publication_year"].astype(str)
    df["_jyi"] = df["_jy"] + ":" + df["_issue"]
    key = df["_jyi"]

    nonempty = df["_issue"] != ""
    issue_counts = (
        df.loc[nonempty].groupby(["_jy", "_jyi"], sort=False).size()
        .rename("articles").reset_index()
    )
    median_by_journal_year = issue_counts.groupby("_jy")["articles"].median()
    counts = df.loc[nonempty].groupby("_jyi")["_jyi"].transform("size")
    med = df.loc[nonempty, "_jy"].map(median_by_journal_year)

    flag = pd.Series(False, index=df.index)
    flag.loc[nonempty] = (
        df.loc[nonempty, "_issue"].str.contains("suppl", case=False, regex=False)
        | (counts >= factor * med)
    )
    flag &= ~key.isin(rescued)
    return flag


def main():
    df = load_raw()
    df = add_title(df)
    n0 = len(df)
    steps = []

    m = df["type"].isin(CFG["article_cleaning"]["doc_types_keep"])
    steps.append(("document type filter", int((~m).sum()))); df = df[m]

    m = ~anomalous_issue_mask(df)
    steps.append(("anomalous-issue (supplement) detection", int((~m).sum()))); df = df[m]

    m = ~title_rule_mask(df)
    steps.append(("title-pattern screening", int((~m).sum()))); df = df[m]

    if CFG["article_cleaning"]["cjk_abstract_screen"]:
        m = ~df["has_cjk"]
        steps.append(("CJK abstract-language screening", int((~m).sum()))); df = df[m]

    out = ROOT / "data" / "articles_clean.parquet"
    df.drop(columns=["_issue", "_jy", "_jyi"]).to_parquet(out, index=False)

    attr = pd.DataFrame(steps, columns=["step", "excluded"])
    attr["remaining"] = n0 - attr["excluded"].cumsum()
    attr.loc[-1] = ["initial screen", 0, n0]
    attr = attr.sort_index().reset_index(drop=True)
    attr.to_csv(ROOT / "data" / "attrition_table.csv", index=False)
    print(attr.to_string(index=False))
    print(f"clean records -> {out}")


if __name__ == "__main__":
    main()
