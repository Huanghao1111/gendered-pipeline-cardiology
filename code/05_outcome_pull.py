#!/usr/bin/env python3
"""Download post-entry works and checkpoint each author group.

Verified linked author IDs are included. Completed groups are skipped on rerun.
"""
import json, os, pickle, time, sys
from pathlib import Path
from collections import defaultdict

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent.parent
OUTF = ROOT / 'data' / 'outcome_pull_2026-09-07.jsonl'
PROG = ROOT / 'data' / 'outcome_pull_progress.txt'
HDR = {'User-Agent': 'gendered-pipeline-cardiology'}
MAILTO = os.environ.get('OPENALEX_MAILTO')

co = pickle.load(open(ROOT / 'data' / 'formal_cohort_2026-09-07.pkl', 'rb'))
v = pd.read_csv(ROOT / 'data' / 'author_id_verification_sheet.csv',
                encoding='utf-8-sig', dtype=str).fillna('')

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
U = defaultdict(set)
for f in parent:
    U[find(f)].add(f)
canon = set(co.fa_id)
shadows = sorted({f for fids in U.values() for f in fids} - canon)
pull_ids = sorted(set(co.fa_id) | set(shadows))
GROUPS = [pull_ids[i:i+10] for i in range(0, len(pull_ids), 10)]
print(f"authors to pull: {len(pull_ids)} in {len(GROUPS)} groups", flush=True)

done = set()
if PROG.exists():
    done = set(PROG.read_text().split())

def trim(w):
    src = ((w.get('primary_location') or {}).get('source') or {})
    auths = []
    for a in (w.get('authorships') or []):
        au = a.get('author') or {}
        ctry = sorted({i.get('country_code') for i in (a.get('institutions') or []) if i.get('country_code')})
        auths.append({'id': au.get('id'), 'name': a.get('raw_author_name') or au.get('display_name'),
                      'pos': a.get('author_position'), 'ctry': ctry})
    return {'work_id': w['id'], 'doi': w.get('doi') or '', 'title': (w.get('title') or '')[:300],
            'year': w.get('publication_year'), 'date': w.get('publication_date') or '',
            'type': w.get('type') or '', 'cited': w.get('cited_by_count'),
            'pmid': (w.get('ids') or {}).get('pmid') or '',
            'source_id': src.get('id') or '', 'source': src.get('display_name') or '',
            'authorships': auths}

def pull_group(ids, max_pages=60):
    flt = f"author.id:{'|'.join(i.split('/')[-1] for i in ids)},publication_year:2005-2025"
    sel = "id,doi,title,publication_year,publication_date,type,cited_by_count,ids,primary_location,authorships"
    cursor, pages = '*', 0
    while True:
        j = None
        for att in range(200):                      # patient: ride out shared-IP rate-limit windows
            try:
                r = requests.get("https://api.openalex.org/works",
                                 params={"filter": flt, "select": sel, "per-page": 200, "cursor": cursor,
                                         "mailto": MAILTO},
                                 headers=HDR, timeout=45)
                if r.status_code == 429:
                    time.sleep(min(120, 20 * (att + 1))); continue
                r.raise_for_status(); j = r.json(); break
            except Exception:
                time.sleep(min(120, 8 * (att + 1)))
        if j is None:
            raise RuntimeError("group failed permanently")
        for w in j.get('results', []):
            yield w
        pages += 1
        cursor = (j.get('meta') or {}).get('next_cursor')
        if not cursor or not j.get('results') or pages >= max_pages:
            return

t0 = time.time()
n_works = 0
with open(OUTF, 'a', encoding='utf-8') as fo:
    for gi, grp in enumerate(GROUPS):
        if str(gi) in done:
            continue
        works = list(pull_group(grp))
        for w in works:
            fo.write(json.dumps(trim(w), ensure_ascii=False) + '\n')
        fo.flush()
        with open(PROG, 'a') as pf:
            pf.write(f"{gi}\n")
        n_works += len(works)
        if gi % 10 == 0:
            print(f"group {gi}/{len(GROUPS)} | works so far {n_works} | {time.time()-t0:.0f}s", flush=True)
        time.sleep(0.3)
print(f"DONE groups={len(GROUPS)} works_added={n_works} total={time.time()-t0:.0f}s", flush=True)
