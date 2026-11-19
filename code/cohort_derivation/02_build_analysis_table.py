#!/usr/bin/env python3
"""Assemble the base analysis cohort from linked records and gender sources.

Sources:
- cohort_linked_2026-09-07.pkl   (linked outcomes, pre-entry covariates, inst_tier)
- formal_cohort_2026-09-07.pkl   (entry paper, institution, first_pub_year)
- gender_probs_formal_2026-09-07.csv (calibrated-model input p_female)
- reference_standard_sheet_2026-09-07.csv (locked sampling and adjudication metadata)
- entry_institution_map_2026-09-07.pkl (fa_id -> institution metadata)
"""
import os
import pickle
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from entry_anchored import prepare_revised_cohort, require_run_authorization

ROOT = Path(os.environ.get('GENDERED_PIPELINE_RUN_DIR', Path.cwd())).resolve()
D = ROOT / 'data'
ANALYSIS_MODE = os.environ.get("ANALYSIS_MODE", "baseline_algorithm")
signed, facts, authorization = require_run_authorization(
    os.environ["SIGNED_ENTRY_TABLE"], os.environ["PREENTRY_FACTS"],
    os.environ["RUN_AUTHORIZATION_MANIFEST"],
    test_mode=os.environ.get("ENTRY_ANCHORED_TEST_MODE") == "YES",
)
if authorization["mode"] != ANALYSIS_MODE:
    raise ValueError("Authorization mode does not match ANALYSIS_MODE")

link = pd.read_pickle(D / 'cohort_linked_2026-09-07.pkl')
co = pd.read_pickle(D / 'formal_cohort_2026-09-07.pkl')
if ANALYSIS_MODE == "revised":
    co = prepare_revised_cohort(co, signed)
    co = co.loc[co['final_decision'].isin(
        ['baseline_unchanged', 'retain', 'uncertain_baseline_unchanged']
    )].copy()
gp = pd.read_csv(D / 'gender_probs_formal_2026-09-07.csv', encoding='utf-8-sig')
rs = pd.read_csv(D / 'reference_standard_sheet_2026-09-07.csv', encoding='utf-8-sig')
eimap = pickle.load(open(D / 'entry_institution_map_2026-09-07.pkl', 'rb'))

df = co.rename(columns={'journal': 'index_journal'})[
    ['fa_id', 'fa_name', 'index_journal', 'entry_year', 'entry_date',
     'primary_cohort', 'institution']].copy()
age_source = 'earliest_verified_publication_year' if ANALYSIS_MODE == 'revised' else 'first_pub_year'
df['academic_age'] = df['entry_year'].astype(int) - co[age_source].astype(int)

# Institution country from the entry-institution map.
cmap = {f: (v[2] if len(v) > 2 and v[2] else '未知') for f, v in eimap.items()}
df['inst_country'] = df['fa_id'].map(cmap).fillna('未知')

# Name scores, locked sampling metadata, and definite reference labels.
df = df.merge(
    gp[['fa_id', 'p_female_raw']].rename(columns={'p_female_raw': 'p_female'}),
    on='fa_id', how='left', validate='1:1',
)
reference_fields = [
    'fa_id', 'core_selected', 'special_selected', 'design_weight',
    'final_name_origin', 'final_status',
]
missing_fields = [field for field in reference_fields if field not in rs.columns]
if missing_fields:
    raise ValueError(
        'The locked reference-standard sheet is missing required fields: '
        + ', '.join(missing_fields)
    )
df = df.merge(rs[reference_fields], on='fa_id', how='left', validate='1:1')
for field in ('core_selected', 'special_selected'):
    df[field] = df[field].fillna(0).astype(int)
gold = rs.loc[rs['gender_final'].isin(['F', 'M']), ['fa_id', 'gender_final']]
df = df.merge(
    gold.rename(columns={'gender_final': 'gender_gold'}),
    on='fa_id', how='left', validate='1:1',
)
df['female_gold'] = df['gender_gold'].map({'F': 1.0, 'M': 0.0})
df['in_rs'] = df['gender_gold'].notna().astype(int)

# Outcomes, pre-entry covariates, and institution tier from linked records.
keep = ['fa_id', 'event_type', 'time_years', 'grad_event', 'grad_time_years',
        'grad_tier', 'grad_year', 'cess_event', 'n_inactivity_spells',
        'inactive_years', 'n_post_works', 'preentry_pubs', 'preentry_citations',
        'preentry_hindex', 'preentry_coauthors', 'intl_collab_share',
        'inst_id', 'inst_tier']
df = df.merge(link[keep], on='fa_id', how='left', validate='1:1')
df['inst_id'] = df['inst_id'].fillna('')

cols = ['fa_id', 'fa_name', 'index_journal', 'entry_year', 'entry_date',
        'primary_cohort', 'academic_age', 'institution', 'inst_country',
        'inst_tier', 'preentry_pubs', 'preentry_citations', 'preentry_hindex',
        'preentry_coauthors', 'intl_collab_share', 'p_female', 'gender_gold',
        'female_gold', 'in_rs', 'core_selected', 'special_selected',
        'design_weight', 'final_name_origin', 'final_status',
        'event_type', 'time_years', 'grad_event', 'grad_time_years',
        'grad_tier', 'grad_year', 'cess_event', 'n_inactivity_spells',
        'inactive_years', 'n_post_works', 'inst_id']
df = df[cols].sort_values('fa_id').reset_index(drop=True)
out = sys.argv[1] if len(sys.argv) > 1 else str(D / 'analysis_cohort_base.csv')
df.to_csv(out, index=False, encoding='utf-8-sig')
print(f"Base analysis cohort assembled: {df.shape}")
