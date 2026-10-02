"""RECONSTRUCTION AID, NOT THE AUTHORS' CODE.

Number by number comparison of the regenerated outputs of code/primeside.py (run with the
shim polya.py) and code/plots.py with the archived files in data/, in the manner of
misc/repro/compare.py (relative difference |a-b|/max(|a|,|b|), 60 digit arithmetic).
Usage: python compare_primeside.py BOOK_DIR RUN_DIR   (writes compare_primeside.json)
"""
import json, os, sys
from collections import defaultdict
import mpmath as mp
import numpy as np

mp.mp.dps = 60
TOL = mp.mpf('1e-9')
book_dir, run_dir = sys.argv[1], sys.argv[2]


def num(x):
    if isinstance(x, bool):
        return None
    if isinstance(x, (int, float)):
        return mp.mpf(x)
    return None


def walk(a, b, path, rows, other):
    if isinstance(a, dict) and isinstance(b, dict):
        for k in a:
            if k not in b:
                other.append((path + '/' + str(k), 'missing in rerun'))
                continue
            walk(a[k], b[k], path + '/' + str(k), rows, other)
        for k in b:
            if k not in a:
                other.append((path + '/' + str(k), 'extra in rerun'))
        return
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            other.append((path, f'length {len(a)} vs {len(b)}'))
        for i, (x, y) in enumerate(zip(a, b)):
            walk(x, y, f'{path}[{i}]', rows, other)
        return
    na, nb = num(a), num(b)
    if na is not None and nb is not None:
        d = mp.mpf(0) if na == nb else abs(na - nb) / max(abs(na), abs(nb), mp.mpf(10) ** -300)
        u = 0.0
        if isinstance(a, float) and isinstance(b, float) and a != 0:
            u = float((b - a) / np.spacing(abs(a)))
        rows.append((path, a, b, d, u))
        return
    if a != b:
        other.append((path, f'{str(a)[:40]} vs {str(b)[:40]}'))


summary = {}
for name in ['results_prime.json', 'extra.json', 'sens2.json', 'results_main.json']:
    fb, fr = os.path.join(book_dir, name), os.path.join(run_dir, name)
    print(f'== {name}')
    if not os.path.exists(fr):
        print('   not produced by this run (no script in the run writes it)')
        summary[name] = 'not produced'
        continue
    a, b = json.load(open(fb)), json.load(open(fr))
    rows, other = [], []
    walk(a, b, '', rows, other)
    n = len(rows)
    nbit = sum(1 for r in rows if r[3] == 0)
    nok = sum(1 for r in rows if r[3] <= TOL)
    worst = sorted(rows, key=lambda r: r[3], reverse=True)
    maxrel = worst[0][3] if rows else mp.mpf(0)
    print(f'   {n} numbers: {nbit} bitwise equal, {nok} within rel 1e-9, {n - nok} beyond; max rel diff {mp.nstr(maxrel, 3)}; structural differences: {len(other)}')
    for p, msg in other[:10]:
        print(f'   structural: {p}: {msg}')
    for p, x, y, d, u in worst[:8]:
        if d == 0:
            break
        print(f'   {p}: book {x!r} | rerun {y!r} (rel {mp.nstr(d, 3)}, {u:+.0f} ulp)')
    bykey = defaultdict(list)
    for p, x, y, d, u in rows:
        key = p.split('/')[-1].split('[')[0] if '/' in p else p
        bykey[key].append((d, u))
    keysum = {}
    for k, v in bykey.items():
        keysum[k] = dict(n=len(v), bitwise_equal=sum(1 for d, u in v if d == 0),
                         max_rel=float(max(d for d, u in v)), max_ulp=max(abs(u) for d, u in v))
        print(f'   field {k:12s}: {len(v):4d} numbers, {keysum[k]["bitwise_equal"]:4d} bitwise equal, max rel {keysum[k]["max_rel"]:.3e}, max {keysum[k]["max_ulp"]:.0f} ulp')
    summary[name] = dict(n=n, bitwise_equal=nbit, within_1e9=nok, max_rel=float(maxrel), structural=other,
                         by_field=keysum,
                         worst=[dict(path=p, book=x, rerun=y, rel=float(d), ulp=u) for p, x, y, d, u in worst[:20]])
json.dump(summary, open('compare_primeside.json', 'w'), indent=1)
