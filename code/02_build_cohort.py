#!/usr/bin/env python3
"""Step 3 — Build the author-level cohort.

First author per article (OpenAlex authorships[author_position=first]);
earliest qualifying paper per OpenAlex author ID; academic age =
entry year - year of first indexed publication (OpenAlex author record);
retain academic age <= 5 (config). Institution = first affiliation at entry.
Attrition counted at every step (author-level flow continues the
article-level table).
"""
import os
from pathlib import Path

import pandas as pd
import requests
import yaml

ROOT = Path(__file__).parent.parent
CFG = yaml.safe_load((ROOT / "config" / "config.yaml").read_text(encoding="utf-8"))
MAILTO = (f"mailto={os.environ['OPENALEX_MAILTO']}"
          if os.environ.get("OPENALEX_MAILTO") else "")


def first_author(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in df.iterrows():
        auths = r.get("authorships") or []
        fa = next((a for a in auths if a.get("author_position") == "first"), None)
        author = (fa or {}).get("author") or {}
        insts = [i.get("display_name") for i in (fa or {}).get("institutions", [])]
        rows.append({
            "work_id": r["id"], "journal": r["_journal"],
            "entry_year": r["publication_year"],
            "entry_date": r.get("publication_date"),
            "fa_id": author.get("id"),
            "fa_name": author.get("display_name"),
            "fa_orcid": author.get("orcid"),
            "institution": insts[0] if insts else None,
            "corresponding": bool((fa or {}).get("is_corresponding")),
            "single_author_entry": len(auths) == 1,
            "indeterminate_first_author": fa is None,
        })
    return pd.DataFrame(rows)


def first_publication_year(author_id: str) -> int | None:
    """Return the earliest publication year on the OpenAlex author record.

    The API fallback is for testing; reported analyses use the frozen snapshot.
    """
    url = (f"https://api.openalex.org/works?filter=author.id:{author_id}"
           f"&sort=publication_year:asc&per-page=1&{MAILTO}")
    j = requests.get(url, timeout=60).json()
    res = j.get("results") or []
    return res[0]["publication_year"] if res else None


def main():
    df = pd.read_parquet(ROOT / "data" / "articles_clean.parquet")
    fa = first_author(df)
    n_articles = len(df)

    locked_path = ROOT / "data" / "entry_eligibility_adjudicated.csv"
    if locked_path.exists():
        locked = pd.read_csv(locked_path)
        assert locked.work_id.is_unique
        fa = fa.merge(
            locked[["work_id", "final_decision", "final_exclusion_reason"]],
            on="work_id", how="left", validate="one_to_one",
        )
    else:
        fa["final_decision"] = ""
        fa["final_exclusion_reason"] = ""

    fa["entry_exclusion"] = ""
    fa.loc[fa.single_author_entry, "entry_exclusion"] = "single_author_entry"
    no_id = fa.fa_id.isna() & fa.entry_exclusion.eq("")
    fa.loc[no_id, "entry_exclusion"] = "indeterminate_first_author"
    locked_group = fa.final_decision.eq("exclude") & fa.final_exclusion_reason.eq("group_byline")
    locked_indeterminate = fa.final_decision.eq("exclude") & fa.final_exclusion_reason.eq("indeterminate_first_author")
    fa.loc[locked_group, "entry_exclusion"] = "group_byline"
    fa.loc[locked_indeterminate, "entry_exclusion"] = "indeterminate_first_author"

    excluded = fa.loc[fa.entry_exclusion.ne(""), ["work_id", "entry_exclusion"]]
    excluded.to_csv(ROOT / "data" / "entry_exclusions.csv", index=False)
    attrition = fa.entry_exclusion.value_counts()
    before = len(fa)
    fa = fa.loc[fa.entry_exclusion.eq("") & fa.fa_id.notna()].copy()
    fa = fa.sort_values(["fa_id", "entry_date"]).drop_duplicates("fa_id", keep="first")
    n_dup = before - int(attrition.sum()) - len(fa)

    fa["first_pub_year"] = fa["fa_id"].map(first_publication_year)
    fa["academic_age"] = fa["entry_year"] - fa["first_pub_year"]
    n_over_age = int((fa["academic_age"] > CFG["cohort"]["academic_age_max"]).sum())
    cohort = fa[fa["academic_age"] <= CFG["cohort"]["academic_age_max"]].copy()
    cohort["primary_cohort"] = cohort["entry_year"].isin(CFG["cohort"]["entry_years_primary"])

    cohort.to_parquet(ROOT / "data" / "cohort.parquet", index=False)
    print(f"articles {n_articles}; entry exclusions: {attrition.to_dict()}; "
          f"duplicate authors -{n_dup}; academic age > {CFG['cohort']['academic_age_max']} "
          f"-{n_over_age}; cohort {len(cohort)} "
          f"(primary subcohort {int(cohort['primary_cohort'].sum())})")


if __name__ == "__main__":
    main()
