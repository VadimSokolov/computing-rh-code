# sqfree_quadratic.py: Polya kernels of the Dedekind zeta functions of imaginary quadratic fields, Section sec:sqfree:quad
# (merge of Nick Polson's note on the squarefree series, misc/reviews/nick-squarefree/, Section 8 there; recomputed independently,
# October 2026), Tables tab:sqfree:named and tab:sqfree:quad and Theorem thm:sqfree:inert. For a fundamental discriminant d < 0 let chi_d be the Kronecker symbol, a_n = sum_{m|n} chi_d(m)
# the number of ideals of norm n, theta_K(x) = sum_n a_n e^{-2 pi n x/sqrt|d|}, H_K(x) = (x^2 theta_K'(x))' and Phi_K(u) = 2 e^{u/2} H_K(e^u),
# so that xi_K(1/2 + s) = int Phi_K(u) e^{su} du with xi_K(s) = s(s-1) 2 (sqrt|d|/2pi)^s Gamma(s) zeta(s) L(s, chi_d). The script
#  (1) checks, at 30 digits, for d = -4 the Fourier identity at z = 0 and 2.5 (trapezoid rule, step 1/200 on [0, 4.5]), the evenness of
#      Phi_K at u = 0.3, 0.6 by direct summation on both sides, theta_K(1/x) = x theta_K(x) + (h/w)(x - 1) at x = 1.5, 2, and for d = -23 and
#      d = -163 the symmetry H_K(1/x) = x H_K(x) at x = 1.25;
#  (2) computes zeta_K(1/2) = zeta(1/2) L(1/2, chi_d) for d = -3, -4, -7, -8 and the sixteen discriminants named in the note;
#  (3) evaluates the bound F(y) + sum_{n>=4} a*_n F(ny) of Theorem thm:sqfree:inert, F(k) = k(k - 2) e^{-k}, a*_n = 0 unless 2 and 3 divide n
#      to even powers and then the number of divisors of n with the 2s and 3s removed, at y = 1.42 (sum to n = 2000, tail bounded by
#      tau(n) <= 2 sqrt n), and the interval of y on which it is negative;
#  (4) scans all fundamental discriminants -20000 <= d < 0: Phi_K > 0 exactly when H_K > 0 on [1, sqrt|d|/pi] (by the symmetry, and
#      since every term is positive beyond), i.e. when sum_n a_n F(n y) > 0 for 2 pi/sqrt|d| <= y <= 2; this is evaluated in float64
#      with the terms n <= 62/y_min + 1 (the rest is below 10^{-23}) on 3000 geometrically spaced y, the least value refined by golden
#      section search; least values within 10^{-8} of zero are recomputed at 30 digits. With y = 2 pi x/sqrt|d| the sum is exactly H_K(x),
#      and Phi_K(log x) = 2 sqrt x H_K(x) is proportional to sqrt y times the sum. It records the counts by chi_d(2), chi_d(3), by blocks of
#      2000 in |d|, the least value of H_K on [1, sqrt|d|/pi] and where, and how many negative kernels have the least value of H_K on
#      [1, sqrt|d|/pi] at the centre x = 1, and the least value of Phi_K on u >= 0 at u = 0.
# Run from any folder: python3 sqfree_quadratic.py; writes sqfree_quadratic.json there (about ten minutes).
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


def s(x, n=10):
    return mp.nstr(x, n)


def is_squarefree(m):
    m = abs(m)
    p = 2
    while p * p <= m:
        if m % (p * p) == 0:
            return False
        p += 1
    return True


def is_fundamental(d):
    if d % 4 == 1:
        return is_squarefree(d)
    if d % 4 == 0:
        m = d // 4
        return m % 4 in (2, 3) and is_squarefree(m)
    return False


NSPF = 5000
spf = list(range(NSPF + 1))
for p in range(2, math.isqrt(NSPF) + 1):
    if spf[p] == p:
        for q in range(p * p, NSPF + 1, p):
            if spf[q] == q:
                spf[q] = p


def chi_table(d, N):
    # Kronecker symbol (d/n), n = 0, ..., N, for a fundamental discriminant d < 0
    def chip(p):
        if p == 2:
            if d % 2 == 0:
                return 0
            return 1 if d % 8 == 1 else -1
        r = d % p
        if r == 0:
            return 0
        return 1 if pow(r, (p - 1) // 2, p) == 1 else -1
    ch = [0] * (N + 1)
    ch[1] = 1
    cache = {}
    for n in range(2, N + 1):
        p = spf[n]
        if p not in cache:
            cache[p] = chip(p)
        ch[n] = cache[p] * ch[n // p]
    return ch


def a_table(d, N):
    ch = chi_table(d, N)
    a = np.zeros(N + 1)
    for m in range(1, N + 1):
        if ch[m]:
            a[m::m] += ch[m]
    return a, ch


# (1) d = -4 at 30 digits
mp.mp.dps = 30


def kernel_mp(d, NA=4000):
    a, ch = a_table(d, NA)
    ad = abs(d)
    kap = 2 * mp.pi / mp.sqrt(ad)

    def theta(x):
        tot, n = mp.mpf(0), 1
        while n <= NA:
            t = kap * n * x
            if a[n]:
                tot += int(a[n]) * mp.exp(-t)
            if t > 120:
                break
            n += 1
        return tot

    def H(x):
        tot, n = mp.mpf(0), 1
        while n <= NA:
            t = kap * n * x
            if a[n]:
                tot += int(a[n]) * t * (t - 2) * mp.exp(-t)
            if t > 120:
                break
            n += 1
        return tot
    return theta, H, ch


theta4, H4, ch4 = kernel_mp(-4)


def Phi4(u):
    return 2 * mp.exp(u / 2) * H4(mp.exp(u))


def xiK(z, d, ch):
    ad = abs(d)
    chi = [ch[n % ad] if n % ad else 0 for n in range(ad)]
    return z * (z - 1) * 2 * (mp.sqrt(ad) / (2 * mp.pi)) ** z * mp.gamma(z) * mp.zeta(z) * mp.dirichlet(z, chi)


hstep = mp.mpf(1) / 200
nodes = [(j * hstep, Phi4(j * hstep)) for j in range(0, 901)]


def fourier4(z):
    return 2 * hstep * (mp.fsum(p * mp.cos(z * u) for u, p in nodes) - nodes[0][1] / 2)


chk = {}
for z in (mp.mpf(0), mp.mpf('2.5')):
    chk['z = ' + s(z, 3)] = [s(fourier4(z), 20), s(xiK(mp.mpf(1) / 2 + 1j * z, -4, ch4), 20)]
chk['Phi_K(-u) - Phi_K(u), u = 0.3, 0.6'] = [s(Phi4(-mp.mpf('0.3')) - Phi4(mp.mpf('0.3')), 3), s(Phi4(-mp.mpf('0.6')) - Phi4(mp.mpf('0.6')), 3)]
chk['theta(1/x) - x theta(x) - (1/4)(x - 1), x = 1.5, 2'] = [s(theta4(1 / mp.mpf(x)) - x * theta4(mp.mpf(x)) - (mp.mpf(x) - 1) / 4, 3) for x in (mp.mpf('1.5'), mp.mpf(2))]
for d in (-23, -163):
    th, H, ch = kernel_mp(d)
    chk['d = %d: H(1/x) - x H(x) at x = 1.25, and H(1.25)' % d] = [s(H(1 / mp.mpf('1.25')) - mp.mpf('1.25') * H(mp.mpf('1.25')), 3), s(H(mp.mpf('1.25')), 8)]
rec('d = -4 kernel checks (integral, xi_K)', chk)

# (2) zeta_K(1/2)
named_pos = [-11, -15, -20, -23, -24, -31, -39, -47, -71, -191]
named_neg = [-19, -43, -67, -163, -427, -1003]
zh = {}
for d in [-3, -4, -7, -8] + named_pos + named_neg:
    ad = abs(d)
    ch = chi_table(d, ad)
    zh[d] = s(mp.zeta(mp.mpf(1) / 2) * mp.dirichlet(mp.mpf(1) / 2, [ch[n] if n else 0 for n in range(ad)]), 8)
rec('zeta_K(1/2)', zh)

# (3) the inert bound


def tau(n):
    t, p = 1, 2
    while p * p <= n:
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        t *= e + 1
        p += 1
    return t * (2 if n > 1 else 1)


def astar(n):
    v2 = v3 = 0
    while n % 2 == 0:
        n //= 2
        v2 += 1
    while n % 3 == 0:
        n //= 3
        v3 += 1
    return 0 if (v2 % 2 or v3 % 2) else tau(n)


AS = [0] + [astar(n) for n in range(1, 2001)]


def F(k):
    return k * (k - 2) * mp.exp(-k)


def bound(y):
    return F(y) + mp.fsum(AS[n] * F(n * y) for n in range(4, 2001))


y0 = mp.mpf('1.42')
tailb = mp.nsum(lambda n: 2 * mp.sqrt(n) * (n * y0) ** 2 * mp.exp(-n * y0), [2001, mp.inf])
ys = [mp.mpf('0.5') + mp.mpf(j) / 1000 for j in range(1501)]
bv = [bound(y) for y in ys]
negs = [y for y, b in zip(ys, bv) if b < 0]
jm = min(range(len(ys)), key=lambda j: bv[j])
rec('inert bound at y = 1.42: F(1.42), sum_{n>=4} a*_n F(1.42 n), total, tail bound beyond n = 2000',
    [s(F(y0), 8), s(bound(y0) - F(y0), 8), s(bound(y0), 8), s(tailb, 3)])
rec('inert bound: negative for y in [a, b] on the grid of step 0.001 in [0.5, 2]; least value and where',
    [s(negs[0], 4), s(negs[-1], 4), s(bv[jm], 6), s(ys[jm], 4)])

# (4) the scan
M = 3000


def Hgrid(a, y):
    n = np.arange(1, len(a), dtype=np.float64)
    av = a[1:]
    res = np.empty(len(y))
    for i0 in range(0, len(y), 500):
        ny = np.outer(y[i0:i0 + 500], n)
        res[i0:i0 + 500] = (ny * (ny - 2) * np.exp(-ny)) @ av
    return res


def golden_np(f, a, b, it=60):
    gr = (math.sqrt(5) - 1) / 2
    c, d_ = b - gr * (b - a), a + gr * (b - a)
    fc, fd = f(c), f(d_)
    for _ in range(it):
        if fc < fd:
            b, d_, fd = d_, c, fc
            c = b - gr * (b - a)
            fc = f(c)
        else:
            a, c, fc = c, d_, fd
            d_ = a + gr * (b - a)
            fd = f(d_)
    return (a + b) / 2


discs = [d for d in range(-1, -20001, -1) if is_fundamental(d)]
table = {}
blocks = {}
rows = []
centreH = centrePhi = borderline = 0
for d in discs:
    ad = -d
    ch2 = 0 if d % 2 == 0 else (1 if d % 8 == 1 else -1)
    ch3 = 0 if d % 3 == 0 else (1 if pow(d % 3, 1, 3) == 1 else -1)
    ymin = 2 * math.pi / math.sqrt(ad)
    blk = (ad - 1) // 2000
    if ymin >= 2:
        pos, mH, xs, cH, cP = 1, None, None, 0, 0
    else:
        N = int(62 / ymin) + 2
        a, _ = a_table(d, N)
        y = np.geomspace(ymin, 2, M)
        hv = Hgrid(a, y)
        f1 = lambda t: float(Hgrid(a, np.array([t]))[0])
        j = int(np.argmin(hv))
        ys_ = golden_np(f1, y[max(j - 1, 0)], y[min(j + 1, M - 1)]) if 0 < j < M - 1 else y[j]
        mH = min(f1(ys_), float(hv[j]))
        if abs(mH) < 1e-8:
            borderline += 1
            with mp.workdps(30):
                kap = 2 * mp.pi / mp.sqrt(ad)
                mH = float(mp.fsum(int(a[n]) * (n * mp.mpf(ys_)) * (n * mp.mpf(ys_) - 2) * mp.exp(-n * mp.mpf(ys_)) for n in range(1, N + 1) if a[n]))
        pos = 1 if mH > 0 else 0
        xs = ys_ * math.sqrt(ad) / (2 * math.pi)
        cH = 1 if j == 0 else 0
        cP = 1 if int(np.argmin(np.sqrt(y) * hv)) == 0 else 0
        if not pos:
            centreH += cH
            centrePhi += cP
    key = '%d,%d' % (ch2, ch3)
    t = table.setdefault(key, [0, 0])
    t[0] += pos
    t[1] += 1
    b = blocks.setdefault(blk, [0, 0])
    b[0] += pos
    b[1] += 1
    rows.append([d, ch2, ch3, pos, None if mH is None else float('%.6g' % mH), None if xs is None else float('%.5g' % xs), cH, cP])
npos = sum(r[3] for r in rows)
rec('fundamental discriminants -20000 <= d < 0, positive kernels', [len(discs), npos])
rec('positive/total by chi(2), chi(3)', table)
rec('positive/total by blocks of 2000 in |d| (block k is 2000k < |d| <= 2000(k+1))', {k: blocks[k] for k in sorted(blocks)})
rec('negative kernels with the least value of H_K at x = 1, and of Phi_K at u = 0; least values recomputed at 30 digits',
    [centreH, centrePhi, borderline])
rec('d = 5 mod 24 (2 and 3 inert): count and positive', [sum(1 for r in rows if r[1] == -1 and r[2] == -1), sum(r[3] for r in rows if r[1] == -1 and r[2] == -1)])
nm = {r[0]: r for r in rows if r[0] in named_pos + named_neg}
COLS = '[d, chi(2), chi(3), positive, least H_K on [1, sqrt|d|/pi], at x, H_K least at x = 1, Phi_K least at u = 0]'
rec('named discriminants: ' + COLS, [nm[d] for d in named_pos + named_neg])
out['rows: ' + COLS] = rows
out['seconds'] = time.time() - t0
json.dump(out, open('sqfree_quadratic.json', 'w'), indent=1)
print('done', round(time.time() - t0, 1))
