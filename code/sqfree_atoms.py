# sqfree_atoms.py: Figure fig:sqfree:atoms of Section sec:sqfree (ch/sqfree.tex): the normalised squarefree error
# Delta(x)/x^{1/4}, Delta(x) = Q(x) - 6x/pi^2, against its truncated explicit formula. Added for the book in October 2026, from
# Nick Polson's note on the squarefree series (misc/reviews/nick-squarefree/, figure squarefree_atoms.pdf there), recomputed
# independently. The truncated explicit formula is
#   1 + 2 Re sum_{k<=K} R_k x^{rho_k/2} + sum_{k=1,3,5} r_k x^{-k},  R_k = zeta(rho_k/2)/(rho_k zeta'(rho_k)),  r_k = zeta(-k)/(-2k zeta'(-2k)),
# the residues of D(s) x^s/s at s = 0, at the poles rho_k/2 of 1/zeta(2s) and at the trivial poles -1, -3, -5; the first 200 zeros are
# from Arb (acb.zeta_zeros) and the residues from mpmath at 30 digits. Q(x) is sieved exactly, and both curves are sampled on the grid
# x = 100 + j/20, j = 0, ..., 58000. The script records the Pearson correlation of the two samples and the root mean square of their
# difference for K = 10, 20, 60, 100, 200, and draws the figure with K = 60.
# Run from any folder: python3 sqfree_atoms.py; writes sqfree_atoms.json and sqfree_atoms.pdf there (the book's copy is
# fig/sqfree_atoms.pdf).
import json
import math
import time

import matplotlib
import mpmath as mp
import numpy as np
from flint import acb, ctx

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

t0 = time.time()
mp.mp.dps = 30
ctx.prec = 150
KMAX = 200
GAM = [mp.mpf(z.imag.mid().str(40, radius=False)) for z in acb.zeta_zeros(1, KMAX)]
R = []
for g in GAM:
    rho = mp.mpc(mp.mpf(1) / 2, g)
    R.append(complex(mp.zeta(rho / 2) / (rho * mp.zeta(rho, derivative=1))))
Rv = np.array(R)
half = np.array([complex(0.25, float(g) / 2) for g in GAM])
triv = [(k, float(mp.zeta(-k) / (-2 * k * mp.zeta(-2 * k, derivative=1)))) for k in (1, 3, 5)]

NMAX = 3000
sf = np.ones(NMAX + 1, dtype=np.bool_)
sf[0] = False
for p in range(2, math.isqrt(NMAX) + 1):
    sf[p * p::p * p] = False
Q = np.cumsum(sf)
x = 100 + np.arange(58001) / 20.0
c = 6 / math.pi ** 2
exact = (Q[np.floor(x).astype(int)] - c * x) / x ** 0.25
lx = np.log(x)


def approx(K):
    tot = np.ones_like(x)
    for k in range(K):
        tot += 2 * np.real(Rv[k] * np.exp(half[k] * lx))
    for k, v in triv:
        tot += v * x ** (-k)
    return tot / x ** 0.25


out = {'R_1, R_2, R_3': [[z.real, z.imag] for z in R[:3]], 'trivial residues r_1, r_3, r_5': triv, 'by K': {}}
for K in (10, 20, 60, 100, 200):
    a = approx(K)
    out['by K'][K] = {'correlation': float(np.corrcoef(exact, a)[0, 1]), 'rms of the difference': float(np.sqrt(np.mean((exact - a) ** 2))),
                      'rms of the exact curve': float(np.sqrt(np.mean(exact ** 2)))}
    print(K, out['by K'][K], flush=True)
plt.rcParams.update({'font.size': 9, 'axes.grid': True, 'grid.alpha': 0.3})
fig, ax = plt.subplots(figsize=(5.0, 2.8))
ax.plot(x, exact, color='0.55', lw=0.6, label=r'$\Delta(x)/x^{1/4}$')
ax.plot(x, approx(60), color='navy', lw=0.9, label='explicit formula, 60 zeros')
ax.set_xlabel(r'$x$')
ax.set_xlim(100, 3000)
ax.legend(fontsize=7.5, loc='lower right')
fig.tight_layout()
fig.savefig('sqfree_atoms.pdf')
out['seconds'] = time.time() - t0
json.dump(out, open('sqfree_atoms.json', 'w'), indent=1)
print('done', round(time.time() - t0, 1))
