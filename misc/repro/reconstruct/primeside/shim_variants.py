"""RECONSTRUCTION AID, NOT THE AUTHORS' CODE.

Forensic check for the missing module `polya` of code/primeside.py: which evaluation
convention reproduces the last bits of the 'polya' column of data/results_prime.json?
The column is Re(xi'/xi)(alpha + i b) = (F1[0]/F[0]).real for the 18 rows (alpha, b).
Candidates:
  flint_<prec>_<div>: code/polya_flint.py, setup(prec, 200, 3), then the ratio formed
      in ball arithmetic ('ball', as code/plots.py does for data/extra.json), by Python
      complex division ('pyc') or by numpy complex128 division ('np', as the shim does);
  mp_<dps>_<div>: the same trapezoid rule (h = 1/200, |u| <= 3) summed in mpmath at
      dps digits (primeside.py sets mp.mp.dps = 30), ratio in mpmath ('mp') or numpy ('np').
Inputs: book_results_prime.json and book_extra.json (copies of the archived data files).
Output: shim_variants.json and a table on stdout.
"""
import json, time
import numpy as np, mpmath as mp
from flint import acb, ctx
from polya_flint import setup, xi_pair

book = json.load(open('book_results_prime.json'))
extra = json.load(open('book_extra.json'))['prime']
refs = {'book_polya': [r['polya'] for r in book], 'extra_polya': [r['polya'] for r in extra],
        'book_exact': [r['exact'] for r in book]}
pts = [(r['alpha'], r['b']) for r in book]


def ulps(v, r):
    return float((v - r) / np.spacing(abs(r)))


cands = {}
for prec in [240, 420]:
    t0 = time.time()
    S = setup(prec, 200, 3)
    ball, pyc, npc = [], [], []
    for a, b in pts:
        ctx.prec = prec
        A, B = xi_pair(S, acb(a - 0.5, b))
        ball.append(float((B / A).real))
        pyc.append((complex(B) / complex(A)).real)
        npc.append(float((np.array([complex(B)]) / np.array([complex(A)]))[0].real))
    cands[f'flint_{prec}_ball'] = ball
    cands[f'flint_{prec}_pyc'] = pyc
    cands[f'flint_{prec}_np'] = npc
    print(f'flint {prec} bits: {time.time() - t0:.1f}s', flush=True)


def mp_nodes(h_inv=200, umax=3):
    h = mp.mpf(1) / h_inv
    J = umax * h_inv
    Phi = []
    for j in range(J + 1):
        u = j * h
        e2 = mp.exp(2 * u)
        s = mp.mpf(0)
        n = 1
        while True:
            t = (2 * mp.pi ** 2 * n ** 4 * mp.exp(9 * u / 2) - 3 * mp.pi * n ** 2 * mp.exp(5 * u / 2)) * mp.exp(-mp.pi * n * n * e2)
            s += t
            if n > 2 and abs(t) < mp.mpf(2) ** (-mp.mp.prec - 20):
                break
            n += 1
        Phi.append(2 * s)
    return h, Phi


for dps in [30, 50]:
    t0 = time.time()
    mp.mp.dps = dps
    h, Phi = mp_nodes()
    mpd, npd = [], []
    for a, b in pts:
        s = mp.mpc(a - 0.5, b)
        A = h * Phi[0]
        B = mp.mpf(0)
        for j in range(1, len(Phi)):
            u = j * h
            ep = mp.exp(s * u)
            em = 1 / ep
            A += h * Phi[j] * (ep + em)
            B += h * u * Phi[j] * (ep - em)
        mpd.append(float(mp.re(B / A)))
        npd.append(float((np.array([complex(B)]) / np.array([complex(A)]))[0].real))
    cands[f'mp_{dps}_mp'] = mpd
    cands[f'mp_{dps}_np'] = npd
    print(f'mpmath dps {dps}: {time.time() - t0:.1f}s', flush=True)

out = {}
print(f"{'candidate':18s} " + ' '.join(f'{k:>28s}' for k in refs))
for name, vals in cands.items():
    row = {}
    cells = []
    for rk, rv in refs.items():
        u = [ulps(v, r) for v, r in zip(vals, rv)]
        nexact = sum(1 for x in u if x == 0)
        maxrel = max(abs(v - r) / abs(r) for v, r in zip(vals, rv))
        row[rk] = dict(bitwise_equal=nexact, max_rel=maxrel, ulps=u)
        cells.append(f'{nexact:2d}/18 eq, max rel {maxrel:8.2e}')
    out[name] = dict(values=vals, vs=row)
    print(f'{name:18s} ' + ' '.join(f'{c:>28s}' for c in cells), flush=True)
print('\nulps of each candidate against the archived polya column (rows in file order):')
for name in cands:
    print(f'{name:18s}', ' '.join(f'{x:+.0f}' for x in out[name]['vs']['book_polya']['ulps']))
print(f"{'extra.json polya':18s}", ' '.join(f'{ulps(v, r):+.0f}' for v, r in zip(refs['extra_polya'], refs['book_polya'])))
print(f"{'book exact':18s}", ' '.join(f'{ulps(v, r):+.0f}' for v, r in zip(refs['book_exact'], refs['book_polya'])))
json.dump(dict(points=pts, candidates=out), open('shim_variants.json', 'w'), indent=1)
