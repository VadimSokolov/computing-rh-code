# Compare archived outputs (ref/) with a fresh run (run/) and write compare6_report.md in BASE.
# BASE is the first argument (default: the Part VI run).
import json, glob, os, math, sys
import numpy as np
BASE = sys.argv[1] if len(sys.argv) > 1 else '/scratch/vsokolov/rh_book_repro/part6'
REF, RUN = f'{BASE}/ref', f'{BASE}/run'

def num(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool)

def walk(a, b, path, out):
    if isinstance(a, dict):
        if not isinstance(b, dict):
            out.append((path, 'type differs', None)); return
        for k in a:
            if k not in b:
                out.append((f'{path}/{k}', 'missing in run', None))
            else:
                walk(a[k], b[k], f'{path}/{k}', out)
        for k in b:
            if k not in a:
                out.append((f'{path}/{k}', 'new in run', None))
    elif isinstance(a, list):
        if not isinstance(b, list) or len(a) != len(b):
            out.append((path, f'length {len(a)} vs {len(b) if isinstance(b, list) else type(b).__name__}', None)); return
        try:
            A = np.array(a, dtype=float); Bv = np.array(b, dtype=float)
        except (ValueError, TypeError):
            A = None
        if A is not None and A.shape == Bv.shape and A.size:
            both_nan = np.isnan(A) & np.isnan(Bv)
            d = np.where(both_nan, 0.0, np.abs(A - Bv))
            scale = np.maximum(np.abs(A), np.abs(Bv))
            rel = np.where(d == 0, 0.0, d / np.where(scale == 0, 1.0, scale))
            out.append((path, 'num', (float(np.nanmax(d)), float(np.nanmax(rel)), int(A.size)))); return
        for i, (x, y) in enumerate(zip(a, b)):
            walk(x, y, f'{path}[{i}]', out)
    elif num(a) and num(b):
        if (isinstance(a, float) and math.isnan(a)) and (isinstance(b, float) and math.isnan(b)):
            d = 0.0
        else:
            d = abs(a - b)
        rel = 0.0 if d == 0 else d / max(abs(a), abs(b))
        out.append((path, 'num', (d, rel, 1)))
    else:
        out.append((path, 'equal' if a == b else f'differs: {str(a)[:50]} vs {str(b)[:50]}', None))

lines = [f'# Reproducibility: archived outputs against a fresh run on Hopper ({os.path.basename(BASE)})', '']
status = json.load(open(f'{BASE}/status.json'))
lines += ['## Script runs', ''] + [f'- {k}: {v}' for k, v in sorted(status.items())] + ['']
lines += ['## Outputs', '', 'For each archived file: the number of values compared, the largest absolute and relative differences, and the worst entries (relative difference above 1e-9).', '']
for ref in sorted(glob.glob(f'{REF}/**/*.json', recursive=True) + glob.glob(f'{REF}/**/*.npy', recursive=True)):
    rel_path = os.path.relpath(ref, REF)
    run = f'{RUN}/{rel_path}'
    if not os.path.exists(run):
        lines.append(f'### {rel_path}: NOT PRODUCED by the run'); lines.append(''); continue
    out = []
    if ref.endswith('.npy'):
        a, b = np.load(ref), np.load(run)
        walk(a.tolist(), b.tolist(), '', out)
    else:
        walk(json.load(open(ref)), json.load(open(run)), '', out)
    nums = [o for o in out if o[1] == 'num']
    other = [o for o in out if o[1] not in ('num', 'equal')]
    n = sum(o[2][2] for o in nums)
    mabs = max([o[2][0] for o in nums], default=0.0); mrel = max([o[2][1] for o in nums], default=0.0)
    verdict = 'identical' if mabs == 0 and not other else ('agrees to 1e-9' if mrel <= 1e-9 and not other else 'DIFFERS')
    lines.append(f'### {rel_path}: {verdict}')
    lines.append(f'{n} values; max abs diff {mabs:.3e}; max rel diff {mrel:.3e}')
    worst = sorted([o for o in nums if o[2][1] > 1e-9], key=lambda o: -o[2][1])[:12]
    for p, _, (d, r, k) in worst:
        lines.append(f'- {p or "/"}: abs {d:.3e}, rel {r:.3e} ({k} values)')
    for p, what, _ in other[:12]:
        lines.append(f'- {p or "/"}: {what}')
    lines.append('')
open(f'{BASE}/compare6_report.md', 'w').write('\n'.join(lines) + '\n')
print('\n'.join(lines))
