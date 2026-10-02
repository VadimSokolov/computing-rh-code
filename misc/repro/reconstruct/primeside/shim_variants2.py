"""RECONSTRUCTION AID, NOT THE AUTHORS' CODE.

Second forensic pass for the missing module `polya` of code/primeside.py. The first pass
(shim_variants.py) showed that every accurate evaluation of Re(xi'/xi)(alpha + i b) gives
the correctly rounded value (equal, bit for bit, to the 'exact' column), while the archived
'polya' column differs from it by -5 to +2 ulp. So the lost evaluator carried an error of
about 1e-15. This pass looks for a working precision and summation order that reproduces
those last bits:
  flint_<nodes>_<eval>: code/polya_flint.py with nodes at <nodes> bits, polynomial
      evaluation (acb_poly, as in the archive) at <eval> bits;
  mp_<nodes>_<eval>_<how>: the same trapezoid rule in mpmath, nodes at <nodes> digits,
      sums at <eval> digits, by Horner's rule in E = e^{sh} ('horner'), by powers of E from
      repeated multiplication ('powers') or with one exponential per node ('direct').
The ratio is formed in ball or mpmath arithmetic ('*_hp') and by numpy complex128
division of the two values rounded to double ('*_np').
Inputs: book_results_prime.json. Output: shim_variants2.json and the best candidates on stdout.
"""
import json, time
import numpy as np, mpmath as mp
from flint import acb, ctx
from polya_flint import setup, xi_pair

book = json.load(open('book_results_prime.json'))
ref = [r['polya'] for r in book]
pts = [(r['alpha'], r['b']) for r in book]
res = {}


def record(name, vals):
    u = [float((v - r) / np.spacing(abs(r))) for v, r in zip(vals, ref)]
    res[name] = dict(eq=sum(1 for x in u if x == 0), maxulp=max(abs(x) for x in u), ulps=u, values=vals)


t0 = time.time()
for pn in [240, 420]:
    S = setup(pn, 200, 3)
    for pe in [53, 64, 72, 80, 88, 96, 100, 104, 108, 112, 116, 120, 128, 160]:
        ctx.prec = pe
        hp, npd = [], []
        for a, b in pts:
            A, B = xi_pair(S, acb(a - 0.5, b))
            hp.append(float((B / A).real))
            npd.append(float((np.array([complex(B)]) / np.array([complex(A)]))[0].real))
        record(f'flint_{pn}_{pe}_hp', hp)
        record(f'flint_{pn}_{pe}_np', npd)
for p in [96, 100, 104, 112, 128]:
    S = setup(p, 200, 3)
    hp, npd = [], []
    for a, b in pts:
        ctx.prec = p
        A, B = xi_pair(S, acb(a - 0.5, b))
        hp.append(float((B / A).real))
        npd.append(float((np.array([complex(B)]) / np.array([complex(A)]))[0].real))
    record(f'flint_{p}_{p}_hp', hp)
    record(f'flint_{p}_{p}_np', npd)
print(f'flint variants {time.time() - t0:.1f}s', flush=True)


def mp_nodes(h_inv=200, umax=3):
    h = mp.mpf(1) / h_inv
    Phi = []
    for j in range(umax * h_inv + 1):
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


t0 = time.time()
for dn in [30, 50]:
    mp.mp.dps = dn
    h0, Phi0 = mp_nodes()
    for de in [16, 20, 24, 26, 28, 29, 30, 31, 32, 34]:
        mp.mp.dps = de
        h = +h0
        Phi = [+x for x in Phi0]
        cA = [h * Phi[0] / 2] + [h * Phi[j] for j in range(1, len(Phi))]
        cB = [mp.mpf(0)] + [h * (j * h) * Phi[j] for j in range(1, len(Phi))]
        for how in ['horner', 'powers', 'direct']:
            hp, npd = [], []
            for a, b in pts:
                s = mp.mpc(a - 0.5, b)
                if how == 'horner':
                    E = mp.exp(s * h); Ei = 1 / E
                    A = mp.polyval(cA[::-1], E) + mp.polyval(cA[::-1], Ei)
                    B = mp.polyval(cB[::-1], E) - mp.polyval(cB[::-1], Ei)
                elif how == 'powers':
                    E = mp.exp(s * h); Ei = 1 / E
                    A = cA[0] * 2; B = mp.mpf(0); ep = mp.mpc(1); em = mp.mpc(1)
                    for j in range(1, len(cA)):
                        ep *= E; em *= Ei
                        A += cA[j] * (ep + em); B += cB[j] * (ep - em)
                else:
                    A = cA[0] * 2; B = mp.mpf(0)
                    for j in range(1, len(cA)):
                        ep = mp.exp(s * (j * h)); em = 1 / ep
                        A += cA[j] * (ep + em); B += cB[j] * (ep - em)
                hp.append(float(mp.re(B / A)))
                npd.append(float((np.array([complex(B)]) / np.array([complex(A)]))[0].real))
            record(f'mp_{dn}_{de}_{how}_hp', hp)
            record(f'mp_{dn}_{de}_{how}_np', npd)
print(f'mpmath variants {time.time() - t0:.1f}s', flush=True)

best = sorted(res.items(), key=lambda kv: (-kv[1]['eq'], kv[1]['maxulp']))
print(f'{len(res)} candidates; best 25 by number of bitwise matches with the archived polya column:')
for name, r in best[:25]:
    print(f"{name:28s} {r['eq']:2d}/18 equal, max {r['maxulp']:.0f} ulp  ", ' '.join(f'{x:+.0f}' for x in r['ulps']))
json.dump(dict(points=pts, archived_polya=ref, candidates=res), open('shim_variants2.json', 'w'), indent=1)
