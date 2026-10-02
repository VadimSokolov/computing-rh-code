# Compare the regenerated JSON outputs with the book's data/*.json.
import json, os, glob, math
import mpmath as mp
mp.mp.dps = 60
BASE = os.environ.get('RH_REPRO_BASE', '/scratch/vsokolov/rh_book_repro')
def num(x):
    if isinstance(x, bool): return None
    if isinstance(x, (int, float)): return mp.mpf(x)
    if isinstance(x, str):
        try: return mp.mpf(x.strip().replace('[', '').split('+/-')[0])
        except Exception: return None
    return None
def walk(a, b, path, diffs, stats):
    if isinstance(a, dict) and isinstance(b, dict):
        for k in a:
            if k not in b: diffs.append((path + '/' + str(k), 'missing in rerun', '')); continue
            walk(a[k], b[k], path + '/' + str(k), diffs, stats)
        for k in b:
            if k not in a: diffs.append((path + '/' + str(k), 'extra in rerun', ''))
        return
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b): diffs.append((path, f'length {len(a)} vs {len(b)}', ''))
        for i, (x, y) in enumerate(zip(a, b)): walk(x, y, f'{path}[{i}]', diffs, stats)
        return
    na, nb = num(a), num(b)
    if na is not None and nb is not None:
        stats['n'] += 1
        if na == nb: return
        d = abs(na - nb) / max(abs(na), abs(nb), mp.mpf(10) ** -300)
        stats['maxrel'] = max(stats['maxrel'], d)
        if d > mp.mpf('1e-9'): diffs.append((path, mp.nstr(na, 15), mp.nstr(nb, 15) + f'  (rel {mp.nstr(d, 3)})'))
        return
    if a != b: diffs.append((path, str(a)[:60], str(b)[:60]))
out = ['# Reproducibility comparison: book data/*.json against a fresh run on Hopper', '']
st = json.load(open(f'{BASE}/status.json'))
out.append('## Script runs'); out.append('')
for k, v in sorted(st.items()): out.append(f'- {k}: {v}')
out.append(''); out.append('## Output comparison (relative tolerance 1e-9)'); out.append('')
# the data files of data/ against run/, then the JSON outputs of verification_reciprocal_zeta/ against verif/
for f, new, name in [(f, f'{BASE}/run/' + os.path.basename(f), os.path.basename(f)) for f in sorted(glob.glob(f'{BASE}/book_data/*.json'))] + \
        [(f, f'{BASE}/verif/' + os.path.basename(f), 'verif/' + os.path.basename(f)) for f in sorted(glob.glob(f'{BASE}/book_verif/*.json'))]:
    if not os.path.exists(new): out.append(f'- {name}: NOT REGENERATED'); continue
    a = json.load(open(f)); b = json.load(open(new)); diffs = []; stats = {'n': 0, 'maxrel': mp.mpf(0)}
    walk(a, b, '', diffs, stats)
    if not diffs: out.append(f'- {name}: agrees ({stats["n"]} numbers, max rel diff {mp.nstr(stats["maxrel"], 3)})')
    else:
        out.append(f'- {name}: {len(diffs)} differences beyond 1e-9 among {stats["n"]} numbers (max rel {mp.nstr(stats["maxrel"], 3)}); first ones:')
        for p, x, y in diffs[:12]: out.append(f'    - {p}: book {x} | rerun {y}')
for pair in [('verif/enclosures_report.md', 'book_verif/enclosures_report.md'), ('logs/spot_checks.out', 'book_verif/spot_checks_report.md'), ('logs/bosch_v7.out', 'book_bosch/verify_bosch_xi_v7_output.txt')]:
    n, o = f'{BASE}/{pair[0]}', f'{BASE}/{pair[1]}'
    if os.path.exists(n) and os.path.exists(o):
        same = open(n).read() == open(o).read()
        out.append(f'- {pair[0]}: ' + ('identical to the archived report' if same else 'DIFFERS from the archived report (see diff)'))
    else: out.append(f'- {pair[0]}: not produced')
open(f'{BASE}/compare_report.md', 'w').write('\n'.join(out) + '\n'); print('\n'.join(out))
