"""Runs the cases of sources over seeds, across processes, and records
each case as a line of JSON."""

import json
import multiprocessing
import os
import shutil
import sys
import time
import traceback
from typing import Dict, Iterable, List, Optional, Tuple

from .cases import run_case
from .sources import SOURCES, Env


def _work(item) -> List[dict]:
    name, seed, root, vsearch, bindir, keep = item
    src = SOURCES[name]
    T = os.path.join(root, f'{name}-{seed}')
    shutil.rmtree(T, ignore_errors=True)
    os.makedirs(T)
    base = {'source': name, 'seed': seed}
    try:
        cases = src.make(T, seed, Env(bindir, vsearch))
    except Exception:
        return [dict(base, case='(generator)', outcome='error',
                     reason=traceback.format_exc(limit=3))]
    out = []
    for k, case in enumerate(cases):
        t = time.time()
        try:
            rec = run_case(case, os.path.join(T, f'c{k}'), vsearch, bindir, keep)
        except Exception:
            rec = {'case': case.name, 'outcome': 'error',
                   'reason': traceback.format_exc(limit=3)}
        rec = dict(base, index=k, **rec)
        rec['wall'] = round(time.time() - t, 3)
        out.append(rec)
    if not keep and all(r['outcome'] in ('pass', 'skip') for r in out):
        shutil.rmtree(T, ignore_errors=True)
    return out


def run(names: List[str], seeds: Optional[List[int]], jobs: int, root: str,
        vsearch: str, bindir: str, out_path: Optional[str], keep: bool,
        quiet: bool) -> Dict[str, Dict[str, int]]:
    def default(n):
        if n == 'datasets':
            from .sources.datasets import default_seeds
            return default_seeds()
        return SOURCES[n].seeds
    items = [(n, s, root, vsearch, bindir, keep)
             for n in names for s in (seeds or default(n))]
    os.makedirs(root, exist_ok=True)
    tally: Dict[str, Dict[str, int]] = {}
    fout = open(out_path, 'a') if out_path else None
    shown = 0
    try:
        with multiprocessing.Pool(jobs) as pool:
            for recs in pool.imap_unordered(_work, items):
                for r in recs:
                    t = tally.setdefault(r['source'], {})
                    t[r['outcome']] = t.get(r['outcome'], 0) + 1
                    if fout:
                        fout.write(json.dumps(r, sort_keys=True) + '\n')
                    if r['outcome'] in ('fail', 'error') and shown < 200:
                        shown += 1
                        what = r.get('diffs') or [r.get('reason', '')]
                        print(f"DIFF {r['source']} seed {r['seed']} "
                              f"[{r['case']}]: {'; '.join(what)[:400]}",
                              flush=True)
                if fout:
                    fout.flush()
    finally:
        if fout:
            fout.close()
    return tally


def print_tally(tally: Dict[str, Dict[str, int]]):
    total = {}
    for name in sorted(tally):
        t = tally[name]
        for k, v in t.items():
            total[k] = total.get(k, 0) + v
        print(f"{name:12} " + ' '.join(f'{k}={t[k]}' for k in
                                       ('pass', 'fail', 'skip', 'error') if k in t))
    if len(tally) > 1:
        print(f"{'total':12} " + ' '.join(f'{k}={total[k]}' for k in
                                          ('pass', 'fail', 'skip', 'error')
                                          if k in total))


def summary(paths: Iterable[str]):
    """Outcomes per source, and the coverage tallies of the cases."""
    tally: Dict[str, Dict[str, int]] = {}
    cover: Dict[str, Dict[str, Dict[str, int]]] = {}
    caveats: Dict[str, int] = {}
    seconds = {}
    for p in paths:
        with open(p) as f:
            for line in f:
                r = json.loads(line)
                s = r['source']
                t = tally.setdefault(s, {})
                t[r['outcome']] = t.get(r['outcome'], 0) + 1
                c = cover.setdefault(s, {})
                for k, v in (r.get('coverage') or {}).items():
                    if isinstance(v, dict):
                        v = ','.join(f'{a}{b}' for a, b in sorted(v.items()))
                    if isinstance(v, (int, float)) and k not in ('exit',):
                        v = '0' if not v else ('1' if v == 1 else '>1')
                    c.setdefault(k, {})
                    c[k][str(v)] = c[k].get(str(v), 0) + 1
                for cv in r.get('caveats') or []:
                    caveats[cv] = caveats.get(cv, 0) + 1
                sec = r.get('seconds')
                if sec:
                    a = seconds.setdefault(s, [0.0, 0.0])
                    a[0] += sec['vsearch']
                    a[1] += sec['rsearch']
    print_tally(tally)
    print()
    for s in sorted(cover):
        print(f'{s}:')
        for k, vals in sorted(cover[s].items()):
            top = sorted(vals.items(), key=lambda x: -x[1])[:8]
            print(f'  {k}: ' + ', '.join(f'{v}={n}' for v, n in top))
        if s in seconds:
            print(f'  time: vsearch {seconds[s][0]:.1f}s, rsearch {seconds[s][1]:.1f}s')
    if caveats:
        print('\ncaveats:')
        for cv, n in sorted(caveats.items(), key=lambda x: -x[1]):
            print(f'  {n:7}  {cv}')
