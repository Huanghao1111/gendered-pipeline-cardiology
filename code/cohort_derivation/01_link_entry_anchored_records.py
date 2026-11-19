#!/usr/bin/env python3
"""Link entry-anchored records and derive outcomes from frozen pulls.

Per preregistration (config.yaml#outcome):
- graduation = first post-index ORIGINAL ARTICLE as LAST author,
  PubMed-indexed (ids.pmid present), 2–30 authors, no consortium byline;
- inactivity = gaps of >=5 consecutive years without any indexed work
  (reversible; re-entry allowed);
- censoring at 2024-12-31;
- tiers: tier1 = entry journals (EHJ/Circulation/JACC) + NEJM/JAMA/Lancet/
  BMJ/Nature Medicine; tier3 = all other PubMed-indexed;
  tier2 uses the locked SJR-2024 cardiology Q1 source list;
  tier3 covers the remaining PubMed-indexed journals.
Pre-entry covariates use a locked rolling 5-year window:
[entry_date - 5 years, entry_date). The entry work is excluded explicitly.
Apart from the explicit two-author minimum, the original outcome derivation is
unchanged. Formal execution is gated by the signed entry table, pre-entry facts
and their SHA-256 authorization manifest.
"""
import gzip, json, os, pickle, re, sys
from pathlib import Path
from collections import defaultdict

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from entry_anchored import (
    academic_age_audit, augment_author_id_verification, compute_preentry_covariates,
    covariates_from_original_pre, is_graduation_author_count_eligible,
    partition_original_for_entry_anchored, prepare_revised_cohort, read_table,
    require_run_authorization, sha256,
)
from outcome_eligibility import (
    load_adjudications, resolve_eligibility, write_review_log,
)

ROOT = Path(os.environ.get('GENDERED_PIPELINE_RUN_DIR', Path.cwd())).resolve()
PULL = ROOT / 'data' / 'outcome_pull_2026-09-07.jsonl.gz'
PULL2 = ROOT / 'data' / 'preentry_pull_2026-09-07.jsonl.gz'
assert PULL2.exists(), "run code/06_preentry_pull.py first"
ANALYSIS_MODE = os.environ.get("ANALYSIS_MODE", "baseline_algorithm")
if ANALYSIS_MODE not in {"baseline_algorithm", "revised"}:
    raise ValueError("ANALYSIS_MODE must be baseline_algorithm or revised")
signed, facts, authorization = require_run_authorization(
    os.environ["SIGNED_ENTRY_TABLE"], os.environ["PREENTRY_FACTS"],
    os.environ["RUN_AUTHORIZATION_MANIFEST"],
    test_mode=os.environ.get("ENTRY_ANCHORED_TEST_MODE") == "YES",
)
if authorization["mode"] != ANALYSIS_MODE:
    raise ValueError("Authorization mode does not match ANALYSIS_MODE")

co = pickle.load(open(ROOT / 'data' / 'formal_cohort_2026-09-07.pkl', 'rb'))
co['entry_work_id'] = co['work_id']
entry_merges = None
if ANALYSIS_MODE == "revised":
    merge_path = os.environ.get("ENTRY_AUTHOR_ID_MERGES")
    if not merge_path:
        raise ValueError("ENTRY_AUTHOR_ID_MERGES is required for revised analysis")
    expected_merge_hash = authorization.get("author_id_merge_changes_sha256")
    if not expected_merge_hash or sha256(merge_path) != expected_merge_hash:
        raise ValueError("Entry author-ID merge file SHA-256 mismatch")
    entry_merges = read_table(merge_path)
    co = prepare_revised_cohort(co, signed)
    academic_age_audit(signed).to_csv(ROOT / 'data' / 'academic_age_audit.csv', index=False, encoding='utf-8-sig')
co = co.loc[co['final_decision'].isin(
    ['baseline_unchanged', 'retain', 'uncertain_baseline_unchanged']
)].copy() if ANALYSIS_MODE == "revised" else co

v = pd.read_csv(ROOT / 'data' / 'author_id_verification_sheet.csv',
                encoding='utf-8-sig', dtype=str).fillna('')
if ANALYSIS_MODE == "revised":
    v = augment_author_id_verification(v, co, entry_merges)

# ---- merge map (shadow -> canonical), union-find over verified sheet
parent = {}
def find(x):
    parent.setdefault(x, x)
    while parent[x] != x:
        parent[x] = parent[parent[x]]; x = parent[x]
    return x
for pk, grp in v.groupby('pair_key'):
    ids = grp['fa_id'].tolist()
    for b in ids[1:]:
        parent[find(ids[0])] = find(b)
GROUPS = defaultdict(set)
for f in parent:
    GROUPS[find(f)].add(f)
canon = set(co.fa_id)
merge_map = {}
for fids in GROUPS.values():
    c = [f for f in fids if f in canon]
    t = c[0] if c else sorted(fids)[0]
    for f in fids:
        merge_map[f] = t
if ANALYSIS_MODE == "revised":
    retained_by_root = defaultdict(list)
    for retained_id in canon:
        retained_by_root[find(retained_id)].append(retained_id)
    collisions = {root: sorted(ids) for root, ids in retained_by_root.items() if len(ids) > 1}
    if collisions:
        raise ValueError(f"Verified aliases merge multiple retained cohort rows: {collisions}")

TIER1_SRC = {"S181568219", "S116251202", "S60865174"}
TIER1_NAMES = {"new england journal of medicine", "the new england journal of medicine",
               "jama", "the lancet", "lancet", "bmj", "the bmj", "nature medicine"}
# tier2 = SJR-2024 Cardiology & Cardiovascular Medicine Q1,
# stored as OpenAlex source IDs (n=101)
TIER2_SRC = set(pickle.load(open(ROOT / 'data' / 'tier2_source_ids_2026-09-07.pkl', 'rb')))
TIER2_SRC = {s.split('/')[-1] for s in TIER2_SRC}
CONSORTIUM = re.compile(r"(?i)\b(consortium|investigators?|collaborat\w*|committee|trialists?|"
                        r"task force|working group|study group|coalition|alliance|registry|"
                        r"group\b|on behalf of)\b")

CENSOR = "2024-12-31"

# Attribute works to canonical researchers using compact records.
person_works = defaultdict(dict)   # canon -> work_id -> slim dict
for PULLX in (PULL, PULL2):        # disjoint year windows (2005-2025 / 2000-2004)
  with gzip.open(PULLX, 'rt', encoding='utf-8') as f:
      for line in f:
          w = json.loads(line)
          hit = [(merge_map.get(a.get('id'), a.get('id')), a['pos'])
                 for a in w['authorships']
                 if a.get('id') and (a.get('id') in merge_map or a.get('id') in canon)]
          hit = [(c, p) for c, p in hit if c in canon]
          if not hit:
              continue
          slim = {
              'work_id': w['work_id'], 'title': w.get('title') or '',
              'date': w.get('date') or '', 'year': w.get('year'), 'type': w.get('type') or '',
              'pmid': w.get('pmid') or '', 'source_id': w.get('source_id') or '',
              'source': w.get('source') or '', 'cited': w.get('cited'),
              'n_auth': len(w['authorships']),
              'consortium': any(a.get('name') and CONSORTIUM.search(a['name'])
                                for a in w['authorships']),
              'intl': len({c for a in w['authorships'] for c in (a.get('ctry') or [])}) > 1,
              'coauth': {a.get('id') or a.get('name') for a in w['authorships']},
          }
          for c, p in hit:
              slot = person_works[c].setdefault(w['work_id'], [slim, []])
              slot[1].append(p)

decisions = load_adjudications(ROOT / 'data' / 'outcome_eligibility_adjudicated.csv')
review_rows, unresolved = [], set()
rows = []
preentry_audit_rows = []
for _, r in co.iterrows():
    fid = r['fa_id']
    entry_date = r['entry_date']
    entry_year = int(r['entry_year'])
    works = person_works.get(fid, {})

    post, pre, covar_pre = partition_original_for_entry_anchored(r, works)

    # ---- graduation event: programmatic screening, then locked manual eligibility
    cands = []
    for w, poss, d in post:
        if 'last' not in poss:
            continue
        if w['type'] != 'article':
            continue
        if not w['pmid']:
            continue
        if not is_graduation_author_count_eligible(w['n_auth']):
            continue
        if w['consortium']:
            continue
        eligible, required, decision = resolve_eligibility(w['work_id'], w['title'], decisions)
        if required or decision:
            review_rows.append({'fa_id': fid, 'work_id': w['work_id'], 'date': d,
                                'title': w['title'], 'review_required': required,
                                'final_decision': decision})
        if required and not decision:
            unresolved.add(w['work_id'])
        if not eligible:
            continue
        cands.append((d, w))
    cands.sort(key=lambda x: x[0])
    if cands:
        gd, gw = cands[0]
        grad_event, grad_date = 1, gd
        sid = gw['source_id'].split('/')[-1]
        sname = (gw['source'] or '').lower()
        if sid in TIER1_SRC or sname in TIER1_NAMES:
            tier = 'tier1'
        elif sid in TIER2_SRC:
            tier = 'tier2'
        else:
            tier = 'tier3'   # PubMed-indexed (pmid required above), non-Q1-cardiology
        grad_year = gw['year']
    else:
        grad_event, grad_date, tier, grad_year = 0, CENSOR, '', None

    t_years = (pd.Timestamp(grad_date) - pd.Timestamp(entry_date)).days / 365.25

    # ---- inactivity spells (all works, any position)
    years = sorted({w['year'] for w, _, _ in post + pre if w['year']})
    spells, cur = [], 0
    cess_event, cess_time_years = 0, None
    if years:
        seq = [entry_year] + [y for y in years if y > entry_year] + [2025]
        gaps = [(a, b) for a, b in zip(seq, seq[1:])]
        spells = [b - a for a, b in gaps if b - a >= 5]
        big = [(a, b) for a, b in gaps if b - a >= 5]
        if big:
            cess_event = 1
            a0 = big[0][0]   # last active year before the first >=5y gap
            cess_time_years = round((pd.Timestamp(f"{a0}-12-31") - pd.Timestamp(entry_date)).days / 365.25, 3)
    n_spells = len(spells)
    inactive_years = sum(g - 1 for g in spells)   # full calendar years without work

    # ---- pre-entry covariates: locked rolling 5-year window
    if (ANALYSIS_MODE == "revised" and r["final_decision"] == "retain"
            and r.get("covariate_algorithm", "signed_preentry_facts") == "signed_preentry_facts"):
        pre_metrics, fact_audit = compute_preentry_covariates(r, facts)
        preentry_audit_rows.extend(fact_audit)
        n_pre = pre_metrics['preentry_pubs']
        pre_cit = pre_metrics['preentry_citations']
        h = pre_metrics['preentry_hindex']
        coauthor_n = pre_metrics['preentry_coauthors']
        intl_share = pre_metrics['intl_collab_share']
    else:
        n_pre = len(covar_pre)
        cites = [w['cited'] or 0 for w, _, _ in covar_pre]
        pre_cit = sum(cites)
        srt = sorted(cites, reverse=True)
        h = sum(1 for i, c in enumerate(srt, 1) if c >= i)
        coauth = {x for w, _, _ in covar_pre for x in w['coauth']}
        coauth = {x for x in coauth if merge_map.get(x, x) != fid}
        coauthor_n = len(coauth)
        intl = [1 if w['intl'] else 0 for w, _, _ in covar_pre]
        intl_share = (sum(intl) / len(intl)) if intl else None

    # ---- competing-event resolution (A0 cessation before graduation)
    t_cens = (pd.Timestamp(CENSOR) - pd.Timestamp(entry_date)).days / 365.25
    t_grad = t_years if grad_event else float('inf')
    t_cess = cess_time_years if cess_event else float('inf')
    if t_grad <= t_cess and t_grad <= t_cens:
        event_type, time_years = 1, t_grad            # graduation first
    elif t_cess < t_cens:
        event_type, time_years = 2, t_cess            # cessation first
    else:
        event_type, time_years = 0, t_cens            # administrative censoring
    time_years = round(max(time_years, 1/365.25), 3)

    rows.append({'fa_id': fid, 'grad_event': grad_event, 'grad_time_years': round(t_years, 3),
                 'grad_tier': tier, 'grad_year': grad_year,
                 'cess_event': cess_event,
                 'cess_time_years': (cess_time_years if cess_event else None),
                 'event_type': event_type, 'time_years': time_years,
                 'n_post_works': len(post), 'n_inactivity_spells': n_spells,
                 'inactive_years': inactive_years,
                 'preentry_pubs': n_pre, 'preentry_citations': pre_cit,
                 'preentry_hindex': h, 'preentry_coauthors': coauthor_n,
                 'intl_collab_share': (round(intl_share, 4) if intl_share is not None else None)})

write_review_log(review_rows, ROOT / 'data' / 'outcome_eligibility_review.csv')
if unresolved:
    raise ValueError(f"{len(unresolved)} title-flagged outcomes require locked include/exclude "
                     "decisions in data/outcome_eligibility_adjudicated.csv")

if ANALYSIS_MODE == "revised":
    pd.DataFrame(preentry_audit_rows).to_csv(
        ROOT / 'data' / 'preentry_fact_audit.csv', index=False, encoding='utf-8-sig')

link = pd.DataFrame(rows)

# Institution tier based on citations per work
eimap = pickle.load(open(ROOT / 'data' / 'entry_institution_map_2026-09-07.pkl', 'rb'))
ist = pd.read_csv(ROOT / 'data' / 'institution_stats_2026-09-07.csv', encoding='utf-8-sig')
ist['per_work_cited'] = ist['cited'] / ist['works']
ist['pct'] = ist['per_work_cited'].rank(pct=True)
ist['inst_tier'] = ist['pct'].map(lambda p: '高' if p >= 0.75 else ('低' if p < 0.25 else '中'))
link['inst_id'] = link.fa_id.map(lambda f: eimap.get(f, (None,))[0])
link = link.merge(ist[['inst_id', 'inst_name', 'per_work_cited', 'pct', 'inst_tier']],
                  on='inst_id', how='left')
link['inst_tier'] = link['inst_tier'].fillna('缺失')

link.to_csv(ROOT / 'data' / 'cohort_linked_2026-09-07.csv', index=False, encoding='utf-8-sig')
link.to_pickle(ROOT / 'data' / 'cohort_linked_2026-09-07.pkl')
print(f"linked: {len(link)} persons | graduates: {link.grad_event.sum()} "
      f"({link.grad_event.mean():.1%}) | tier1-of-grads: "
      f"{(link.grad_tier=='tier1').sum()}")
print(link[['grad_time_years','n_post_works','preentry_pubs']].describe().round(2).to_string())
