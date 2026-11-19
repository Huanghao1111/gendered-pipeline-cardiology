#!/usr/bin/env python3
"""Download pre-entry works using the outcome-pull ID groups and schema."""
import json, os, pickle, time, sys, threading
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent.parent
WORK_TMP = ROOT / 'work_tmp'
WORK_TMP.mkdir(parents=True, exist_ok=True)
OUTF = WORK_TMP / 'preentry_pull.jsonl'          # local disk during run; copied to data/ at the end
PROG = WORK_TMP / 'preentry_progress.txt'
# seed from any earlier progress stored on the mount
if not PROG.exists():
    import shutil, gzip
    if (ROOT/'data/preentry_pull_progress.txt').exists():
        shutil.copy(ROOT/'data/preentry_pull_progress.txt', PROG)
    if (ROOT/'data/preentry_pull_ckpt.jsonl.gz').exists():
        with gzip.open(ROOT/'data/preentry_pull_ckpt.jsonl.gz', 'rt', encoding='utf-8') as fi, \
             open(OUTF, 'w', encoding='utf-8') as fo:
            for ln in fi:
                fo.write(ln)
HDR = {'User-Agent': 'gendered-pipeline-cardiology'}
MAILTO = os.environ.get('OPENALEX_MAILTO')
lock = threading.Lock()

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
# deviation #18: only authors with first_pub_year < 2005 can have pre-2005 works
co['first_pub_year'] = co['first_pub_year'].astype(int)
need = set(co.loc[co.first_pub_year < 2005, 'fa_id'])
canon_of = {}
for fids in U.values():
    c = [f for f in fids if f in canon]
    t = c[0] if c else sorted(fids)[0]
    for f in fids:
        canon_of[f] = t
shadows = sorted({f for fids in U.values() for f in fids} - canon)
shadows = [f for f in shadows if canon_of.get(f) in need]
pull_ids = sorted(need | set(shadows))
print(f"pre-2005-history authors: {len(need)} + shadows {len(shadows)} -> pulling {len(pull_ids)} IDs", flush=True)
GROUPS = [pull_ids[i:i+10] for i in range(0, len(pull_ids), 10)]

done = set(PROG.read_text().split()) if PROG.exists() else set()
todo = [(gi, grp) for gi, grp in enumerate(GROUPS) if str(gi) not in done]
print(f"total {len(GROUPS)} groups | done {len(done)} | todo {len(todo)}", flush=True)

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
    flt = f"author.id:{'|'.join(i.split('/')[-1] for i in ids)},publication_year:2000-2004"
    sel = "id,doi,title,publication_year,publication_date,type,cited_by_count,ids,primary_location,authorships"
    cursor, pages = '*', 0
    while True:
        j = None
        for att in range(60):
            try:
                prm = {"filter": flt, "select": sel, "per-page": 200,
                       "cursor": cursor, "mailto": MAILTO}
                _key = os.environ.get("OPENALEX_API_KEY", "").strip()
                if _key:
                    prm["api_key"] = _key
                r = requests.get("https://api.openalex.org/works",
                                 params=prm, headers=HDR, timeout=45)
                if r.status_code == 429:
                    time.sleep(min(90, 5 * (att + 1))); continue
                r.raise_for_status(); j = r.json(); break
            except Exception:
                time.sleep(min(60, 4 * (att + 1)))
        if j is None:
            raise RuntimeError("group failed permanently")
        for w in j.get('results', []):
            yield w
        pages += 1
        cursor = (j.get('meta') or {}).get('next_cursor')
        if not cursor or not j.get('results') or pages >= max_pages:
            return

def work(item):
    gi, grp = item
    return gi, [trim(w) for w in pull_group(grp)]

def persist(gi, works):
    for att in range(8):
        try:
            with lock:
                with open(OUTF, 'a', encoding='utf-8') as fo:
                    for w in works:
                        fo.write(json.dumps(w, ensure_ascii=False) + '\n')
                with open(PROG, 'a') as pf:
                    pf.write(f"{gi}\n")
            return True
        except OSError as e:
            print(f"write retry {att}: {e}", flush=True)
            time.sleep(3 * (att + 1))
    return False

def checkpoint():
    """periodic checkpoint from work_tmp to the data directory.
    JSONL is gzipped for compact checkpoints."""
    import shutil, gzip
    for att in range(8):
        try:
            shutil.copy(PROG, ROOT / 'data/preentry_pull_progress.txt')
            gz_tmp = WORK_TMP / 'preentry_pull_ckpt.jsonl.gz'
            with open(OUTF, 'rb') as fi, gzip.open(gz_tmp, 'wb', compresslevel=6) as go:
                shutil.copyfileobj(fi, go, 1 << 20)
            shutil.copy(gz_tmp, ROOT / 'data/outcome_pull_ckpt.jsonl.gz')
            return True
        except OSError as e:
            print(f"ckpt retry {att}: {e}", flush=True); time.sleep(3 * (att + 1))
    return False

from concurrent.futures import ThreadPoolExecutor, as_completed
t0 = time.time(); done_now = 0
with ThreadPoolExecutor(max_workers=2) as ex:
    futs = [ex.submit(work, item) for item in todo]
    for fut in as_completed(futs):
        try:
            gi, works = fut.result()
        except Exception as e:
            print(f"group failed: {e}", flush=True); continue
        if persist(gi, works):
            done_now += 1
        if done_now % 10 == 0:
            print(f"+{done_now} groups this run ({time.time()-t0:.0f}s)", flush=True)
        if done_now and done_now % 20 == 0:
            checkpoint()
print(f"RUN DONE: +{done_now} groups in {time.time()-t0:.0f}s", flush=True)
checkpoint()
left = len(GROUPS) - len(set(PROG.read_text().split()))
print("REMAINING:", left, flush=True)
