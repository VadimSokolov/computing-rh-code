# sqfree_smooth.py: the smoothed squarefree counts of Section sec:sqfree (ch/sqfree.tex), Tables tab:sqfree:smooth and
# tab:sqfree:gamma2. Added for the book in October 2026, from Nick Polson's note on the squarefree series
# (misc/reviews/nick-squarefree/), recomputed independently.
# The exponentially smoothed count S(x) = sum_n mu(n)^2 e^{-n/x} and the Gamma(2) smoothed count E Q(G_X) = sum_n mu(n)^2 (1 + 2n/X) e^{-2n/X}
# (Q_1(x) and Q_2(X) of the section)
# are computed at 45 digits from the rapidly convergent sums over d,
#   S(x) = sum_d mu(d)/(e^{d^2/x} - 1),   E Q(G_X) = sum_d mu(d) (1/(e^t - 1) + t/(4 sinh^2(t/2))),  t = 2 d^2/X,
# (mu sieved to 2.5 x 10^5; terms with d^2/x > 140, resp. t > 140, are below 10^{-55} and omitted), checked against the direct
# sums over the squarefree n <= 60x, resp. n <= 40X, sieved to 4 x 10^7 (float64), and against the truncated sum
# -sum_{d <= 2 x 10^5} mu(d) g(2d^2/X), g(t) = 2/t - 1/(e^t - 1) - t/(4 sinh^2(t/2)), of the note. Both are compared with their explicit
# formulas: 6x/pi^2 + 1 + 2 Re sum_rho A_rho x^{rho/2} + trivial terms, A_rho = Gamma(rho/2) zeta(rho/2)/(2 zeta'(rho)), and
# 6X/pi^2 + 1 + 2 Re sum_rho B_rho X^{rho/2} + trivial terms, B_rho = (1 + rho/2) Gamma(rho/2) zeta(rho/2) 2^{-rho/2}/(2 zeta'(rho)),
# over the first 40 and 100 zeros (Arb, acb.zeta_zeros), with the residues at s = -1, ..., -6 computed as Laurent coefficients
# by the trapezoid rule on circles of radius 0.3 (128 points), and the double pole of Gamma(s) D(s) at s = -1 also in closed form.
# Run from any folder: python3 sqfree_smooth.py; writes sqfree_smooth.json there (a few minutes).
import json
import math
import time

import mpmath as mp
import numpy as np
from flint import acb, ctx

t0 = time.time()
out = {}


def rec(k, v):
    out[k] = v
    print(k, v, round(time.time() - t0, 1), flush=True)


def s(x, n=12):
    return mp.nstr(x, n)


mp.mp.dps = 45
six = 6 / mp.pi ** 2

# Moebius function to 2.5 x 10^5
Dmax = 250000
isp = np.ones(Dmax + 1, dtype=np.bool_)
isp[:2] = False
for p in range(2, math.isqrt(Dmax) + 1):
    if isp[p]:
        isp[p * p::p] = False
mu = np.ones(Dmax + 1, dtype=np.int64)
mu[0] = 0
for p in np.nonzero(isp)[0]:
    p = int(p)
    mu[p::p] *= -1
    mu[p * p::p * p] = 0
DS = [(int(d), int(mu[d])) for d in np.nonzero(mu)[0]]


def S_exp(x):
    x = mp.mpf(x)
    tot = []
    for d, m in DS:
        t = mp.mpf(d) ** 2 / x
        if t > 140:
            break
        tot.append(m / mp.expm1(t))
    return mp.fsum(tot)


def EQ_gamma2(X):
    X = mp.mpf(X)
    tot = []
    for d, m in DS:
        t = 2 * mp.mpf(d) ** 2 / X
        if t > 140:
            break
        tot.append(m * (1 / mp.expm1(t) + t / (4 * mp.sinh(t / 2) ** 2)))
    return mp.fsum(tot)


def g2(t):
    return 2 / t - 1 / mp.expm1(t) - t / (4 * mp.sinh(t / 2) ** 2)


def trunc_note(X, D=200000):
    with mp.workdps(30):
        X = mp.mpf(X)
        return -mp.fsum(m * g2(2 * mp.mpf(d) ** 2 / X) for d, m in DS if d <= D)


# direct sums over the squarefree integers
NS = 4 * 10 ** 7
sf = np.ones(NS + 1, dtype=np.bool_)
sf[0] = False
for p in np.nonzero(isp[:math.isqrt(NS) + 1])[0]:
    p = int(p)
    sf[p * p::p * p] = False
nsf = np.nonzero(sf)[0].astype(np.float64)


def S_direct(x):
    v = nsf[nsf <= 60 * x]
    return math.fsum(np.add.reduceat(np.exp(-v / x), np.arange(0, len(v), 10 ** 6)).tolist())


def EQ_direct(X):
    v = nsf[nsf <= 40 * X]
    w = (1 + 2 * v / X) * np.exp(-2 * v / X)
    return math.fsum(np.add.reduceat(w, np.arange(0, len(w), 10 ** 6)).tolist())


# zeros and residues
ctx.prec = 200
GAM = [mp.mpf(z.imag.mid().str(50, radius=False)) for z in acb.zeta_zeros(1, 100)]
A, Bc = [], []
for gm in GAM:
    rho = mp.mpc(mp.mpf(1) / 2, gm)
    zh, dz, gh = mp.zeta(rho / 2), mp.zeta(rho, derivative=1), mp.gamma(rho / 2)
    A.append(gh * zh / (2 * dz))
    Bc.append((1 + rho / 2) * gh * zh * 2 ** (-rho / 2) / (2 * dz))
RH = [mp.mpc(mp.mpf(1) / 2, gm) / 2 for gm in GAM]
rec('|A_rho| and |B_rho| for the first five zeros, and at the 40th and 100th', [[s(abs(a), 4) for a in A[:5]], [s(abs(b), 4) for b in Bc[:5]],
                                                                                 s(abs(A[39]), 3), s(abs(A[99]), 3)])


def atoms(coef, x, K):
    x = mp.mpf(x)
    return 2 * mp.re(mp.fsum(coef[k] * mp.exp(RH[k] * mp.log(x)) for k in range(K)))


def laurent(f, k, r=mp.mpf('0.3'), M=128):
    # coefficients a_{-1}, a_{-2}, a_{-3} of f at s = -k from (1/M) sum f(s_j) (s_j + k)^m, s_j = -k + r e^{2 pi i j/M}
    pts = [-k + r * mp.expjpi(2 * mp.mpf(j) / M) for j in range(M)]
    vals = [f(z) for z in pts]
    return [mp.fsum(v * (z + k) ** m for v, z in zip(vals, pts)) / M for m in (1, 2, 3)]


def fS(z):
    return mp.gamma(z) * mp.zeta(z) / mp.zeta(2 * z)


def fG(z):
    return (1 + z) * mp.gamma(z) * mp.zeta(z) * 2 ** (-z) / mp.zeta(2 * z)


LS = {k: laurent(fS, k) for k in range(1, 7)}
LG = {k: laurent(fG, k) for k in range(1, 7)}
c1 = -mp.zeta(-1) / (2 * mp.zeta(-2, derivative=1))
c2 = c1 * (mp.digamma(2) + mp.zeta(-1, derivative=1) / mp.zeta(-1) - mp.zeta(-2, derivative=2) / mp.zeta(-2, derivative=1))
rec('Gamma(s) D(s) at s = -k: [a_-1, a_-2, a_-3]', {k: [s(mp.re(v), 15) for v in LS[k]] for k in LS})
rec('(1+s) Gamma(s) D(s) 2^-s at s = -k: [a_-1, a_-2, a_-3]', {k: [s(mp.re(v), 15) for v in LG[k]] for k in LG})
rec('closed form c_1, c_2 against the contour values a_-2, a_-1 at s = -1', [s(c1, 20), s(c2, 20), s(c1 - mp.re(LS[1][1]), 3), s(c2 - mp.re(LS[1][0]), 3)])


def trivial(L, x):
    x = mp.mpf(x)
    return mp.fsum(mp.re(L[k][0] + L[k][1] * mp.log(x)) * x ** (-k) for k in L)


tab = []
for x in (10 ** 3, 10 ** 4, 3 * 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7, 10 ** 8):
    S = S_exp(x)
    err = S - six * x - 1
    a40, a100, tr = atoms(A, x, 40), atoms(A, x, 100), trivial(LS, x)
    row = {'x': x, 'S(x) - 6x/pi^2 - 1': s(err, 12), 'atoms (40)': s(a40, 6), 'atoms (100) - atoms (40)': s(a100 - a40, 3),
           'trivial terms': s(tr, 6), 'c_1 log x/x': s(c1 * mp.log(x) / x, 6), 'c_2/x': s(c2 / x, 6),
           'residual S - 6x/pi^2 - 1 - atoms - trivial': s(err - a100 - tr, 3),
           'residual without c_2/x and later terms (the note\'s comparison)': s(err - a40 - c1 * mp.log(x) / x, 3),
           '(S(x) - 6x/pi^2 - 1)/x^(1/4)': s(err / mp.mpf(x) ** 0.25, 6)}
    if x <= 10 ** 5:
        row['direct sum minus d sum'] = s(S_direct(x) - S, 3)
    tab.append(row)
    print(row, flush=True)
rec('exponential smoothing', tab)
tab = []
for X in (10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7, 10 ** 8):
    E = EQ_gamma2(X)
    err = E - six * X - 1
    a100, tr = atoms(Bc, X, 100), trivial(LG, X)
    row = {'X': X, 'E Q(G_X) - 6X/pi^2': s(E - six * X, 12), '(E Q(G_X) - 6X/pi^2 - 1)/X^(1/4)': s(err / mp.mpf(X) ** 0.25, 6),
           'atoms (100)': s(a100, 6), 'trivial terms': s(tr, 6), '2 c_1/X': s(2 * c1 / X, 6),
           'residual': s(err - a100 - tr, 3)}
    if X <= 10 ** 6:
        row['direct sum - 6X/pi^2'] = s(EQ_direct(X) - six * X, 12)
        row['-sum_{d <= 2 10^5} mu(d) g(2d^2/X)'] = s(trunc_note(X), 12)
    tab.append(row)
    print(row, flush=True)
rec('Gamma(2) smoothing', tab)
out['seconds'] = time.time() - t0
json.dump(out, open('sqfree_smooth.json', 'w'), indent=1)
print('done', round(time.time() - t0, 1))
