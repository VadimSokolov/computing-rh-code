#!/usr/bin/env python3
"""Independent re-derivation of the manuscript's stated constants.

The package as received carries no code, so every number in Section sec:num and
in the prose is unreproducible from the package alone.  This script recomputes
the ones that can be settled from mpmath and the first few hundred zeta zeros,
and prints each beside the value the manuscript prints.  It is not a substitute
for the author's scripts: the torus optimisations of Theorem thm:torus, the
argument principle counts of Corollary cor:Storus, the Talbot inversions of the
density f, and the Meixner Pollaczek expansions are all out of its reach.  (The
torus bounds of Theorem thm:torus have since been proved in ball arithmetic by
the programs in torus/; the argument principle counts still have no script.)

Requires mpmath 1.3.0 and sympy (for the Moebius function).  Runs in about two
minutes at the default precision.

The script writes no file; its standard output is the report:
    python3 spot_checks.py > spot_checks_report.md
"""
import random

from sympy import mobius
from mpmath import (mp, mpf, mpc, zeta, psi, pi, log, sqrt, exp, quad, sech, j,
                    euler, gamma, findroot, diff, mangoldt, siegeltheta, siegelz,
                    zetazero, im, re, erfc, besselk, inf)

mp.dps = 25
FAILURES = []


def check(name, got, claimed, tol=None):
    """Compare a computed value against the manuscript's printed digits."""
    g = mp.nstr(got, 12)
    ok = True
    if claimed is not None and tol is not None:
        ok = abs(mpf(got) - mpf(claimed)) <= mpf(tol)
    print(f"{name:52s} {g:>22s}   paper: {claimed}")
    if not ok:
        FAILURES.append(name)


def section(title):
    print(f"\n--- {title} ---")


section("constants, Section sec:thorin and the constants table")
c = zeta(mpf(1) / 2, derivative=1) / zeta(mpf(1) / 2)
check("c = zeta'(1/2)/zeta(1/2)", c, '2.6860917096', '1e-9')
check("c = (log pi - psi(1/4))/2", (log(pi) - psi(0, mpf(1) / 4)) / 2, '2.6860917096', '1e-9')
check("c = -theta'(0)", -diff(siegeltheta, mpf(0)), '2.6860917096', '1e-9')
check("b = (euler + log pi)/2", (euler + log(pi)) / 2, '0.860973', '1e-6')

tau0 = findroot(lambda t: diff(siegeltheta, t), mpf('6.29'))
check("tau_0, the minimum of theta", tau0, '6.2898360', '1e-6')
check("theta(tau_0)", siegeltheta(tau0), '-3.5309728', '1e-6')
check("max M = -theta(tau_0)/pi", -siegeltheta(tau0) / pi, '1.1239436', '1e-6')
check("tau_pi with theta(tau_pi) = -pi",
      findroot(lambda t: siegeltheta(t) + pi, mpf('3.44')), '3.4362182', '1e-6')

g1 = im(zetazero(1))
check("M(gamma_1^2 -) = -theta(gamma_1)/pi", -siegeltheta(g1) / pi, '0.5502528', '1e-6')
check("f(0+) = -zeta(1/2)/4", -zeta(mpf(1) / 2) / 4, '0.36509', '1e-5')
check("f'(0+) = -zeta(1/2)/16", -zeta(mpf(1) / 2) / 16, '0.09127', '1e-5')
check("tail constant c/(2 sqrt pi)", c / (2 * sqrt(pi)), '0.75773', '1e-5')


section("Lemma lem:W0 and the sign change of W_0")
TH0 = diff(siegeltheta, mpf(0))


def F(x):
    x = mpf(x)
    return 2 * TH0 + quad(lambda u: exp(-u / 4) * (1 - exp(-u ** 2 / (16 * x))) / (1 - exp(-u)),
                          [0, 1, mp.inf])


for xv, cl in [('0.003', '0.0366850'), ('0.0035', '-0.0442285'), ('0.004', '-0.1145624')]:
    check(f"F({xv})", F(xv), cl, '1e-6')

lo, hi = mpf('0.003'), mpf('0.0035')
for _ in range(60):
    mid = (lo + hi) / 2
    if F(mid) > 0:
        lo = mid
    else:
        hi = mid
check("x_*, the zero of W_0", (lo + hi) / 2, '0.0032173906', '1e-9')


section("Theorem thm:w and Theorem thm:id")


def P(x, N=20000):
    x = mpf(x)
    s = mpf(0)
    for n in range(2, N):
        L = mangoldt(n)
        if L:
            s += L / sqrt(n) * exp(-log(n) ** 2 / (4 * x))
    return s / sqrt(4 * pi * x)


check("P(0.004), the small-x prime bound", P('0.004'), '1.99e-13', '1e-14')

ZEROS = [im(zetazero(k)) for k in range(1, 301)]


def W(x):
    x = mpf(x)
    return sum(exp(-x * g ** 2) for g in ZEROS)


def W0(x):
    x = mpf(x)
    return quad(lambda t: exp(-x * t ** 2) * diff(siegeltheta, t), [0, 10, 60, mp.inf]) / pi


check("W_0(0.1), closed form vs direct", W0('0.1'), '-0.86250960', '1e-7')
for xv, cl in [('0.001', '1.00025'), ('0.01', '1.00249'), ('0.1', '0.86251'),
               ('1', '0.49140'), ('10', '0.21524')]:
    x = mpf(xv)
    val = exp(x / 4) - P(xv) if x <= mpf('0.01') else W(x) - W0(x)
    check(f"w({xv})", val, cl, '1e-5')


section("Theorem thm:class: the excursions of S")


def S_at(k, side):
    g = ZEROS[k - 1]
    n_below = k - 1 if side == 'minus' else k
    return n_below - siegeltheta(g) / pi - 1


vals = [(S_at(k, s), k, s) for k in range(1, 301) for s in ('minus', 'plus')]
mn, mx = min(vals), max(vals)
check(f"min S over first 300, at gamma_{mn[1]}", mn[0], '-1.1454808', '1e-6')
check(f"max S over first 300, at gamma_{mx[1]}", mx[0], '1.0975638', '1e-6')


section("Theorem thm:necklace and Proposition prop:signs")


def M_k(k, q):
    return sum(mobius(k // d) * q ** d for d in range(1, k + 1) if k % d == 0) / k


def b_n(n, q):
    s, i = mpf(0), 0
    while n % (2 ** i) == 0:
        s += M_k(n // (2 ** i), q)
        i += 1
    return s


q, y = mpf(2) ** mpf('-0.5'), mpf('0.5')
prod = mpf(1)
for n in range(1, 400):
    prod *= (1 + y ** n) ** b_n(n, q)
print(f"{'binary cyclotomic identity, 1/(1-qy)':52s} {mp.nstr(1/(1-q*y),25)}")
print(f"{'                          product form':52s} {mp.nstr(prod,25)}")
if abs(prod - 1 / (1 - q * y)) > mpf('1e-24'):
    FAILURES.append('binary cyclotomic identity')

for n, expr, lbl in [(2, (q**2 + q) / 2, 'b_2 = (q^2+q)/2'),
                     (3, (q**3 - q) / 3, 'b_3 = (q^3-q)/3'),
                     (4, q**4 / 4 + q**2 / 4 + q / 2, 'b_4'),
                     (6, (q**6 + q**3 - q**2 - q) / 6, 'b_6')]:
    check(lbl, b_n(n, q), mp.nstr(expr, 12), '1e-12')

for p in (2, 3, 5, 7):
    qq = mpf(p) ** mpf('-0.5')
    neg = sum(1 for n in range(2, 61) if b_n(n, qq) < 0)
    check(f"negative b_n for 2<=n<=60, p={p}", mpf(neg), '42', '0')


section("Theorem thm:kappaexact and Corollary cor:Tfalse")


def xi(s):
    s = mpc(s)
    return mpf(1) / 2 * s * (s - 1) * pi ** (-s / 2) * gamma(s / 2) * zeta(s)


def kappa_xi(rho):
    """kappa from the functional-equation form of Proposition prop:firstorder."""
    return rho * xi(3 - rho) / (2 * pi * diff(xi, rho))


def kappa_closed(g):
    """kappa from the closed form of Theorem thm:kappaexact."""
    A = exp(j * siegeltheta(g)) * zeta(mpc(mpf(5) / 2, g))
    Zp = diff(siegelz, g)
    Rk = ((g ** 2 - mpf(15) / 4) * im(A) - 4 * g * re(A)) / (4 * pi ** 2 * Zp)
    Ik = (4 * g * im(A) + (g ** 2 - mpf(15) / 4) * re(A)) / (4 * pi ** 2 * Zp)
    return mpc(Rk, Ik)


for k, cl in [(1, '-4.762-1.944i'), (2, '-8.673-0.627i'),
              (3, '-10.149-3.872i'), (4, '-15.165+3.705i')]:
    rho = mpc(mpf(1) / 2, ZEROS[k - 1])
    a, b = kappa_xi(rho), kappa_closed(ZEROS[k - 1])
    print(f"kappa(rho_{k}): xi form {mp.nstr(a,7):>26s}  closed form {mp.nstr(b,7):>26s}   paper: {cl}")
    if abs(a - b) > mpf('1e-6') * abs(a):
        FAILURES.append(f'kappa forms disagree at rho_{k}')

print("\nscanning the first 460 zeros for Re kappa > 0 (Corollary cor:Tfalse)")
mp.dps = 20
pos = []
for k in range(1, 461):
    g = im(zetazero(k))
    if re(kappa_closed(g)) > 0:
        pos.append((k, g))
print(f"  sign flips found: {len(pos)}   paper: exactly three")
for k, g in pos:
    print(f"    gamma_{k} = {mp.nstr(g,10)}   kappa = {mp.nstr(kappa_closed(g),7)}")
if len(pos) != 3 or [k for k, _ in pos] != [213, 289, 379]:
    FAILURES.append('kappa sign scan')
mp.dps = 25

section("Proposition prop:hmfalse (assumption (H) is false)")
mp.dps = 25
_f = lambda h: quad(lambda u: (1 + h / (2 * u)) * exp(-u - h ** 2 / (4 * u)), [0, h, inf]) * exp(-h)
_sf = lambda h: exp(-h) * quad(lambda u: exp(-u - h ** 2 / (4 * u)), [0, h, inf])

tot = quad(_f, [0, 1, 10, inf])
print(f"int f over (0,inf)            {mp.nstr(tot,12):>22s}   must be 1")
if abs(tot - 1) > mpf('1e-10'):
    FAILURES.append('harmonic mean density does not integrate to 1')

# closed forms, as a cross check on the quadrature: sf = h e^-h K_1(h), f = h e^-h (K_0+K_1)
for h in ['0.25', '1', '3']:
    h = mpf(h)
    a, b = _sf(h), h * exp(-h) * besselk(1, h)
    c, d = _f(h), h * exp(-h) * (besselk(0, h) + besselk(1, h))
    print(f"h={mp.nstr(h,4):>6s}  sf {mp.nstr(a,10):>16s} vs {mp.nstr(b,10):>16s}   f {mp.nstr(c,10):>16s} vs {mp.nstr(d,10):>16s}")
    if abs(a - b) > mpf('1e-12') or abs(c - d) > mpf('1e-12'):
        FAILURES.append(f'harmonic mean quadrature vs Bessel closed form at h={h}')

# Monte Carlo, because the density was derived by hand and one transcription error would hide
random.seed(20260920)
N = 200000
draws = sorted(2 * x * y / (x + y) for x, y in
               ((random.expovariate(1.0), random.expovariate(1.0)) for _ in range(N)))
for h in ['0.1', '1', '2']:
    h = mpf(h)
    mc = mpf(sum(1 for v in draws if v > h)) / N
    ex = _sf(h)
    print(f"h={mp.nstr(h,4):>6s}  P(H>h) Monte Carlo {mp.nstr(mc,6):>10s}   closed {mp.nstr(ex,6):>10s}")
    if abs(mc - ex) > mpf('0.004'):
        FAILURES.append(f'harmonic mean Monte Carlo disagrees at h={h}')

# the contradiction itself: a GGC with zero drift and Thorin mass <= 1 has a nonincreasing
# density, and this one rises above f(0+) = 1
xm = findroot(lambda h: diff(_f, h), mpf('0.1'))
fm = _f(xm)
print(f"f(0+)  1   argmax {mp.nstr(xm,6):>10s}   f(argmax) {mp.nstr(fm,7):>10s}   paper: 1.11125 at 0.09765")
if not (fm > 1) or abs(xm - mpf('0.09765')) > mpf('1e-4') or abs(fm - mpf('1.11125')) > mpf('1e-4'):
    FAILURES.append('harmonic mean density maximum')

# and beta <= 1, which is what lets the lemma apply
for k in (4, 6, 8):
    s = mpf(10) ** k
    phi = quad(lambda h: exp(-s * h) * _f(h), [0, 1 / s, 1, 10, inf])
    print(f"s=1e{k}   s*phi(s) {mp.nstr(s*phi,10):>14s}   psi(s)/log s {mp.nstr(-log(phi)/log(s),10):>14s}   both -> 1")
    if abs(-log(phi) / log(s) - 1) > mpf('0.01'):
        FAILURES.append(f'harmonic mean beta estimate at s=1e{k}')


section("Table tab:symwin (zeros of xi in the doubled windows)")
mp.dps = 20
_g = []
_k = 1
while True:
    _v = im(zetazero(_k))
    _g.append(_v)
    if _v > 440:
        break
    _k += 1
for a, b, claimed in [('0.5', '40', 21), ('40', '70', 27), ('70', '140', 78), ('195', '215', 26)]:
    lo, hi = 2 * mpf(a), 2 * mpf(b)
    idx = [i + 1 for i, v in enumerate(_g) if lo < v < hi]
    print(f"Im w in ({a},{b}) -> gamma in ({mp.nstr(lo,6)},{mp.nstr(hi,6)}): "
          f"{len(idx)} zeros, indices {idx[0]}..{idx[-1]}   table: {claimed}")
    if len(idx) != claimed:
        FAILURES.append(f'xi zero count for window ({a},{b})')
# the last window is the one the paper's K_8 count disagrees with, so record the margin
print(f"  margin at the endpoints: gamma_195 = {mp.nstr(_g[194],8)} < 390, "
      f"gamma_222 = {mp.nstr(_g[221],8)} > 430")
mp.dps = 25


section("Round 2 corrections")

# prop:zetaright: the omitted terms do not enter with weight one, so the truncation
# explains the whole of the row's gap and round 1's "unexplained remainder" was wrong.
mp.dps = 25
_beta, _s = mpf(2), mpf(3)
_u = sqrt(_s)
_lhs = log(zeta(_beta)) - log(zeta(_beta + _u)) + _u * zeta(_beta, derivative=1) / zeta(_beta)
_trunc = mp.fsum(mangoldt(n) * mpf(n) ** (-_beta) * ((1 - mpf(n) ** (-_u)) / log(n) - _u)
                 for n in range(2, 3000))
_wt = (1 - mpf(3000) ** (-_u)) / log(3000) - _u
print(f"log LHS closed form           {mp.nstr(_lhs,12):>18s}   paper: -0.587879")
print(f"series truncated at n<3000    {mp.nstr(_trunc,12):>18s}   paper: -0.587379")
print(f"truncation error              {mp.nstr(_lhs-_trunc,8):>18s}   row's printed gap: -5.00e-4")
print(f"weight carried at n=3000      {mp.nstr(_wt,8):>18s}   round 1 assumed 1")
if abs(_lhs - _trunc) < mpf('5.0e-4'):
    FAILURES.append('prop:zetaright truncation no longer exceeds the printed gap')
if abs(_lhs - mpf('-0.587879')) > mpf('1e-6'):
    FAILURES.append('prop:zetaright left side')

# thm:kappaexact discussion: seven of the first 460 zeros have S(gamma^-) < -1, not two.
mp.dps = 22
_low = []
for k in range(1, 461):
    g = im(zetazero(k))
    Sm = mpf(k - 2) - siegeltheta(g) / pi
    if Sm < -1:
        _low.append((k, Sm + mpf(1) / 2))
print(f"zeros with S(gamma^-) < -1 among the first 460: {len(_low)}   paper now: seven")
print("   " + ", ".join(f"k={k} mid={mp.nstr(m,5)}" for k, m in _low))
_negk = [m for k, m in _low if k not in (213, 289, 379)]
print(f"   largest midpoint in modulus with Re kappa < 0: {mp.nstr(min(_negk),6)}   paper: -0.556")
if len(_low) != 7 or [k for k, _ in _low] != [127, 196, 233, 289, 368, 380, 401]:
    FAILURES.append('S(gamma^-) < -1 enumeration')
if abs(min(_negk) - mpf('-0.5564929')) > mpf('1e-5'):
    FAILURES.append('kappa threshold lower bracket')
mp.dps = 25

# prop:onlyhalf: the five printed minima, and the grid that the proof now discloses.
def _minRe(sig, hi=60, steps=3000):
    best = mpf(10)
    for i in range(1, steps + 1):
        z = mpc(sig, mpf(hi) * i / steps)
        v = re(zeta(z, derivative=1) / zeta(z))
        if v < best:
            best = v
    return best
mp.dps = 18
for sig, claimed in [('0.55', '-1.049'), ('0.6', '-0.995'), ('0.75', '-0.867'),
                     ('0.9', '-0.769'), ('0.99', '-0.721')]:
    got = _minRe(mpf(sig))
    print(f"min Re zeta'/zeta on (0,60], sigma={sig:>5s}   {mp.nstr(got,7):>11s}   paper: {claimed}")
    if abs(got - mpf(claimed)) > mpf('1e-3'):
        FAILURES.append(f'prop:onlyhalf minimum at sigma={sig}')
_worst = max(_minRe(mpf(i) / 40, steps=1500) for i in range(21, 40))
print(f"least negative minimum on the grid step 0.025   {mp.nstr(_worst,6):>11s}   paper: below -0.729")
if _worst > mpf('-0.729'):
    FAILURES.append('prop:onlyhalf grid scan')
mp.dps = 25


section("summary")
print(f"  checks failed: {len(FAILURES)}")
for f in FAILURES:
    print(f"    {f}")
