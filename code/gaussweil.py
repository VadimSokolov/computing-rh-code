# code/gaussweil.py (September 2026, from misc/reviews/nick-audit/B-ray-explicit/weil_table.py): Table tab:ch12:gweil,
# the Gaussian smoothing G(u,t) = sum exp(-2t(gamma~ - u)^2) of the zeros (one zero of each pair {rho, 1 - rho}, gamma > 0,
# gamma~ = gamma - i delta) from the zeros and from the primes, by the explicit formula eq:ch12:Garith:
#   G(u,t) + G(-u,t) = 2 e^{t/2-2tu^2} cos 2tu + A(u,t) - (2 pi t)^{-1/2} sum_{n>=2} Lambda(n) n^{-1/2} e^{-(log n)^2/8t} cos(u log n),
#   A(u,t) = (1/4 pi) int_R (e^{-2t(r-u)^2} + e^{-2t(r+u)^2}) m(r) dr,   m(r) = Re psi(1/4 + ir/2) - log pi,
# which is Weil's formula of Section ch:clock for h(r) = e^{-2t(r-u)^2} + e^{-2t(r+u)^2}, halved.
# Zero side: the first 1700 zeros (zeros_1700.txt, all zeros below 2197.26), with a bound for the zeros above.
# Prime side: Lambda(n) for n <= 1e9 by a segmented sieve, sums in long double. The prime sum beyond N is bounded without
# its phases by psi(x) < 1.03883 x (Rosser and Schoenfeld 1962, Theorem 12): with f(x) = x^{-1/2} e^{-(log x)^2/8t},
# decreasing for x > 1, sum_{n>N} Lambda(n) f(n) <= f(N)(1.03883 N - psi(N)) + 1.03883 e^{t/2} sqrt(2 pi t) erfc((log N - 2t)/sqrt(8t)).
# Also reported: the main term of that tail with its phase, int_N^inf x^{-1/2} e^{-(log x)^2/8t} cos(u log x) dx, which is
# what the prime number theorem predicts for it, in closed form through the Faddeeva function w.
# Reads zeros_1700.txt; writes gaussweil.json. About a minute on 16 cores (Hopper).
import json
import os
from fractions import Fraction
from multiprocessing import Pool

import mpmath as mp
import numpy as np
from scipy.special import wofz

PARAMS = [('100', '1'), ('1977.22', '1.75'), ('1000', '2'), ('100', '3')]
NMAX = 10 ** 9
CHECK = [10 ** 4, 10 ** 6, 10 ** 8, 10 ** 9]
BOUNDS = sorted(set([1] + CHECK + list(range(2 * 10 ** 7, NMAX + 1, 10 ** 7))))
SQ = int(NMAX ** 0.5) + 1
isp = np.ones(SQ + 1, bool); isp[:2] = False
for p in range(2, int(SQ ** 0.5) + 1):
    if isp[p]:
        isp[p * p::p] = False
BASEP = np.nonzero(isp)[0]
LD = lambda x: np.longdouble(Fraction(x).numerator) / np.longdouble(Fraction(x).denominator)   # decimal to long double
UT = [(LD(u), LD(t)) for u, t in PARAMS]


def sums(n, lam):   # sum of Lambda(n) f(n) cos(u log n) for each (u, t), and psi over the block
    ln = np.log(n.astype(np.longdouble))
    amp = lam / np.sqrt(n.astype(np.longdouble))
    return [np.sum(amp * np.exp(-ln * ln / (8 * t)) * np.cos(u * ln)) for u, t in UT], np.sum(lam)


def segment(k):     # the primes in (BOUNDS[k], BOUNDS[k+1]]
    lo, hi = BOUNDS[k] + 1, BOUNDS[k + 1]
    arr = np.ones(hi - lo + 1, bool)
    for p in BASEP:
        p = int(p)
        if p * p > hi:
            break
        arr[max(p * p, -(-lo // p) * p) - lo::p] = False
    if lo <= 1:
        arr[:2 - lo] = False
    pr = np.nonzero(arr)[0].astype(np.int64) + lo
    return (k,) + tuple(sums(pr, np.log(pr.astype(np.longdouble))))


if __name__ == '__main__':
    nseg = len(BOUNDS) - 1
    with Pool(int(os.environ.get('SLURM_CPUS_PER_TASK', '8'))) as pool:
        res = sorted(pool.map(segment, range(nseg)))
    S = np.zeros((nseg, len(PARAMS)), np.longdouble); PSI = np.zeros(nseg, np.longdouble)
    for k, s, psi in res:
        S[k] += s; PSI[k] += psi
    pp, pl = [], []   # prime powers p^j, j >= 2
    for p in BASEP:
        q = int(p) ** 2
        while q <= NMAX:
            pp.append(q); pl.append(int(p)); q *= int(p)
    pp, pl = np.array(pp, np.int64), np.array(pl, np.int64)
    seg = np.searchsorted(np.array(BOUNDS[1:], np.int64), pp, side='left')
    for k in set(seg.tolist()):
        s, psi = sums(pp[seg == k], np.log(pl[seg == k].astype(np.longdouble)))
        S[k] += s; PSI[k] += psi
    S, PSI = np.cumsum(S, axis=0), np.cumsum(PSI)

    mp.mp.dps = 30
    Z = [mp.mpf(x) for x in open('zeros_1700.txt').read().split()]
    T1 = Z[-1]
    m = lambda r: mp.re(mp.digamma(mp.mpf(1) / 4 + 0.5j * r)) - mp.log(mp.pi)
    out = []
    for j, (us, ts) in enumerate(PARAMS):
        u, t = mp.mpf(us), mp.mpf(ts)
        Gp = mp.fsum(mp.exp(-2 * t * (g - u) ** 2) for g in Z)
        Gm = mp.fsum(mp.exp(-2 * t * (g + u) ** 2) for g in Z)
        # zeros above T1 = 2197.26: each term of G(u,t) or G(-u,t) is at most e^{t/2} e^{-2t(gamma - u)^2} in modulus; with
        # N(r) <= r^2 and partial summation, those of G(u,t) + G(-u,t) together at most 2 e^{t/2} e^{-2t(T1-u)^2}
        # (T1^2 + T1/(2t(T1 - u))), for u < T1
        ztail = 2 * mp.exp(t / 2) * mp.exp(-2 * t * (T1 - u) ** 2) * (T1 ** 2 + T1 / (2 * t * (T1 - u)))
        polar = 2 * mp.exp(t / 2 - 2 * t * u * u) * mp.cos(2 * t * u)
        L = 12 / mp.sqrt(t)   # A(u,t) = (1/2 pi) int_R e^{-2t(r-u)^2} m(r) dr, since m is even; |r - u| > L adds < 1e-120
        A = mp.quad(lambda r: mp.exp(-2 * t * (r - u) ** 2) * m(r), [u - L + 2 * L * i / 24 for i in range(25)]) / (2 * mp.pi)
        c = 1 / mp.sqrt(2 * mp.pi * t)
        rows = []
        for N in CHECK:
            k = BOUNDS.index(N) - 1
            prime = c * mp.mpf(str(S[k, j]))
            formula = polar + A - prime
            y = mp.log(N)
            fN = mp.exp(-y / 2 - y * y / (8 * t))
            tail = c * (fN * (mp.mpf('1.03883') * N - mp.mpf(str(PSI[k]))) + mp.mpf('1.03883') * mp.exp(t / 2) * mp.sqrt(2 * mp.pi * t) * mp.erfc((y - 2 * t) / mp.sqrt(8 * t)))
            # int_y^inf e^{(1/2+iu)x - x^2/8t} dx = sqrt(2 pi t) e^{(1/2+iu)y - y^2/8t} w(z), z = (4tu + i(y - 2t))/sqrt(8t), w Faddeeva
            z = (4 * float(t) * float(u) + 1j * (float(y) - 2 * float(t))) / np.sqrt(8 * float(t))
            main = c * mp.re(mp.sqrt(2 * mp.pi * t) * mp.exp(mp.mpc(0.5, u) * y - y * y / (8 * t)) * mp.mpc(complex(wofz(z))))
            rows.append(dict(N=N, prime=mp.nstr(prime, 16), formula=mp.nstr(formula, 16), zeros_minus_formula=mp.nstr(Gp + Gm - formula, 4),
                             tail_bound=mp.nstr(tail, 3), zeros_minus_formula_with_tail_main_term=mp.nstr(Gp + Gm - formula + main, 4)))
        rec = dict(u=us, t=ts, G=mp.nstr(Gp, 16), G_minus_u=mp.nstr(Gm, 4), zeros_above_T1_at_most=mp.nstr(ztail, 3),
                   polar=mp.nstr(polar, 4), archimedean=mp.nstr(A, 16), rows=rows)
        out.append(rec)
        print(json.dumps(rec), flush=True)
    json.dump(dict(psi_1e9=mp.nstr(mp.mpf(str(PSI[-1])), 15), cases=out), open('gaussweil.json', 'w'), indent=1)
