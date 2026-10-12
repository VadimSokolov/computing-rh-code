# sqfree_sieve.py: the sieve computations of Section sec:sqfree (ch/sqfree.tex): Tables tab:sqfree:sharp and tab:sqfree:rms
# and the mean square paragraph. Added for the book in October 2026, from Nick Polson's note on the squarefree series
# (misc/reviews/nick-squarefree/), recomputed independently. The squarefree integers are sieved to N = 2 x 10^8 (numpy, a
# boolean array with the multiples of p^2 struck out for p <= sqrt N), and with c = 6/pi^2, Q(x) = #{n <= x squarefree} and
# Delta(x) = Q(x) - c x the script records
#  (1) at x = 10^3, ..., 10^8: E_1(x) = sum_{n<=x} mu(n)^2/n - c (log x + C), with C = 24 log A - gamma - 2 log 2pi, the normalised
#      error Delta(x)/x^{1/4}, and E_a(x) = sum_{n<=x} mu(n)^2 n^{-a} - c x^{1-a}/(1-a) - zeta(a)/zeta(2a) at a = 3/4, divided by x^{1/2-a};
#  (2) for the dyadic ranges [X, 2X], X = 10^4 2^j, j = 0, ..., 13, and X = 10^4, ..., 10^8, the root mean square (int_X^{2X} Delta^2 dx / X)^{1/2} / X^{1/4},
#      computed exactly from int_n^{n+1} Delta^2 = Delta_n^2 - c Delta_n + c^2/3, Delta_n = Q(n) - c n, and the extremes of
#      Delta(x)/X^{1/4} and of Delta(x)/x^{1/4} on the range (on [n, n+1) Delta decreases from Delta_n to Delta_n - c);
#  (3) the largest and least values of Delta(x)/x^{1/4} on [10^3, N] and where they occur.
# Sums of reciprocals are taken in float64 by numpy's pairwise summation in blocks of 10^6 and the blocks added by math.fsum;
# a long double recomputation checks them. Run from any folder: python3 sqfree_sieve.py; writes sqfree_sieve.json there
# (about 8 GB of memory and a few minutes).
import json
import math
import time

import mpmath as mp
import numpy as np

t0 = time.time()
out = {}


def rec(k, v):
    out[k] = v
    print(k, v, round(time.time() - t0, 1), flush=True)


mp.mp.dps = 30
N = 2 * 10 ** 8
c_mp = 6 / mp.pi ** 2
C_mp = 24 * mp.log(mp.glaisher) - mp.euler - 2 * mp.log(2 * mp.pi)
Dq_mp = mp.zeta(mp.mpf('0.75')) / mp.zeta(mp.mpf('1.5'))
c = float(c_mp)
cC = float(c_mp * C_mp)
Dq = float(Dq_mp)
r = math.isqrt(N)
small = np.ones(r + 1, dtype=np.bool_)
small[:2] = False
for p in range(2, math.isqrt(r) + 1):
    if small[p]:
        small[p * p::p] = False
primes = np.nonzero(small)[0]
sf = np.ones(N + 1, dtype=np.bool_)
sf[0] = False
for p in primes:
    sf[int(p) * int(p)::int(p) * int(p)] = False
rec('sieve done; Q(N)', int(np.count_nonzero(sf)))

checkpoints = [10 ** k for k in range(3, 9)]
ranges = sorted(set([10 ** 4 * 2 ** j for j in range(14)] + [10 ** k for k in range(4, 9)]))   # [X, 2X] with 2X <= N
acc = {X: {'int': [], 'maxabs': (0.0, 0), 'max': (-1e300, 0), 'min': (1e300, 0), 'maxx': (-1e300, 0), 'minx': (1e300, 0)} for X in ranges}
glob = {'max': (-1e300, 0), 'min': (1e300, 0)}
s1_blocks, sq_blocks, s1_ld = [], [], np.longdouble(0)
cp = {}
B = 10 ** 6
Qprev = 0
for start in range(1, N + 1, B):
    end = min(N, start + B - 1)
    blk = sf[start:end + 1]
    n = np.arange(start, end + 1, dtype=np.int64)
    q = Qprev + np.cumsum(blk, dtype=np.int64)
    Qprev = int(q[-1])
    nf = n.astype(np.float64)
    dn = q.astype(np.float64) - c * nf          # Delta_n = Delta(x) for x in [n, n+1) at its left end
    inv = np.where(blk, 1.0 / nf, 0.0)
    invq = np.where(blk, nf ** -0.75, 0.0)
    for x in checkpoints:
        if start <= x <= end:
            i = x - start
            s1 = math.fsum(s1_blocks + [float(np.sum(inv[:i + 1]))])
            sq = math.fsum(sq_blocks + [float(np.sum(invq[:i + 1]))])
            s1l = s1_ld + np.sum(inv[:i + 1].astype(np.longdouble))
            cp[x] = {'Q(x)': int(q[i]), 'Delta(x)': float(dn[i]), 'Delta(x)/x^(1/4)': float(dn[i]) / x ** 0.25,
                     'E_1(x)': s1 - (cC + c * math.log(x)), 'E_1(x), long double sum': float(s1l - np.longdouble(cC) - np.longdouble(c) * np.log(np.longdouble(x))),
                     'E_1(x) x^(3/4)': (s1 - (cC + c * math.log(x))) * x ** 0.75,
                     'E_{3/4}(x)/x^(-1/4)': (sq - 4 * c * x ** 0.25 - Dq) * x ** 0.25}
    s1_blocks.append(float(np.sum(inv)))
    sq_blocks.append(float(np.sum(invq)))
    s1_ld += np.sum(inv.astype(np.longdouble))
    # Delta(x)/x^{1/4} on [n, n+1): at most Delta_n/n^{1/4} when Delta_n > 0, at least about (Delta_n - c)/(n+1)^{1/4}
    if end >= 1000:
        lo = max(0, 1000 - start)
        hi_v = dn[lo:] / nf[lo:] ** 0.25
        lo_v = (dn[lo:] - c) / (nf[lo:] + 1) ** 0.25
        j = int(np.argmax(hi_v))
        if hi_v[j] > glob['max'][0]:
            glob['max'] = (float(hi_v[j]), int(n[lo + j]))
        j = int(np.argmin(lo_v))
        if lo_v[j] < glob['min'][0]:
            glob['min'] = (float(lo_v[j]), int(n[lo + j]) + 1)
    for X in ranges:
        a, b = max(start, X), min(end, 2 * X - 1)
        if a > b:
            continue
        sl = slice(a - start, b - start + 1)
        d = dn[sl]
        A = acc[X]
        A['int'].append(float(np.sum(d * d - c * d + c * c / 3)))
        vals = np.maximum(np.abs(d), np.abs(d - c))
        j = int(np.argmax(vals))
        if vals[j] > A['maxabs'][0]:
            A['maxabs'] = (float(vals[j]), int(n[sl][j]))
        j = int(np.argmax(d))
        if d[j] > A['max'][0]:
            A['max'] = (float(d[j]), int(n[sl][j]))
        j = int(np.argmin(d - c))
        if d[j] - c < A['min'][0]:
            A['min'] = (float(d[j] - c), int(n[sl][j]) + 1)
        xv = nf[sl]
        hv = d / xv ** 0.25
        lv = (d - c) / (xv + 1) ** 0.25
        j = int(np.argmax(hv))
        if hv[j] > A['maxx'][0]:
            A['maxx'] = (float(hv[j]), int(n[sl][j]))
        j = int(np.argmin(lv))
        if lv[j] < A['minx'][0]:
            A['minx'] = (float(lv[j]), int(n[sl][j]) + 1)
rec('constants c, c C, zeta(3/4)/zeta(3/2)', [mp.nstr(c_mp, 20), mp.nstr(c_mp * C_mp, 20), mp.nstr(Dq_mp, 20)])
rec('checkpoints', {str(x): cp[x] for x in checkpoints})
dy = []
for X in ranges:
    A = acc[X]
    I = math.fsum(A['int'])
    X4 = X ** 0.25
    dy.append({'X': X, 'rms/X^(1/4)': math.sqrt(I / X) / X4, 'max|Delta|/X^(1/4)': A['maxabs'][0] / X4, 'at': A['maxabs'][1],
               'max Delta/X^(1/4)': A['max'][0] / X4, 'min Delta/X^(1/4)': A['min'][0] / X4,
               'max Delta/x^(1/4)': A['maxx'][0], 'at max': A['maxx'][1], 'min Delta/x^(1/4)': A['minx'][0], 'at min': A['minx'][1]})
rec('dyadic ranges [X, 2X]', dy)
rec('largest and least Delta(x)/x^(1/4) on [10^3, 2 10^8], with the integer near which they occur', glob)
out['seconds'] = time.time() - t0
json.dump(out, open('sqfree_sieve.json', 'w'), indent=1)
print('done', round(time.time() - t0, 1))
