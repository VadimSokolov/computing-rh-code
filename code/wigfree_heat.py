# wigfree_heat.py: the heat flow as finite free convolution, for the subsection on the heat flow and finite free probability in the
# free probability section (misc/reviews/nick-wigner/merge/wigfree.tex, label sec:wigfree:heat). Added for the book in October 2026,
# from the first author's note on the semicircle law and the Thorin framework (misc/reviews/nick-wigner/source-v7/, version 7,
# Section 11), recomputed independently.
# 1. The identity e^{-(t/2) D^2} p = p boxplus_d q_t, q_t(x) = t^{d/2} He_d(x/sqrt t) = e^{-(t/2) D^2} x^d, with the finite free
#    additive convolution of Marcus, Spielman and Srivastava in coefficient form, and the semigroup q_s boxplus_d q_t = q_{s+t}
#    (sympy, d = 5, generic coefficients).
# 2. The example of the note: p with roots -3, -1, 1/2, 2, 1 +- 0.8 i. The discriminant of e^{-(t/2) D^2} p in x is a polynomial
#    in t (sympy, exact); its least positive root t* is the time at which the complex pair meets the real axis (50 digits), and the
#    largest imaginary part of the roots is recorded at t = 0, 0.1, 0.2, 0.3, t*, 0.35, and the paths of the six roots for
#    0 <= t <= 0.4 are kept for the figure fig/wigfree_heat.pdf (drawn by wigfree_plots.py).
# 3. The finite semicircle: the roots of He_d(sqrt d x), d = 20, 40, 80, 160, 320 (eigenvalues of the Jacobi matrix with
#    off diagonal entries sqrt k, divided by sqrt d; checked against sympy's He_20), against the quantiles q_k of the standard
#    semicircle law, F(q_k) = (k - 1/2)/d; the largest deviation over all roots and over the roots other than the two extreme
#    ones. The roots and quantiles for d = 20 and 160 are kept for the figure.
# Run from any folder: python3 wigfree_heat.py; writes wigfree_heat.json there (under a minute).
import json
import time

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.optimize import brentq

t0 = time.time()
out = {}


def rec(k, v):
    out[k] = v
    print(k, v, round(time.time() - t0, 1), flush=True)


x, t, s = sp.symbols('x t s')


def heat(p, tt, d):
    # e^{-(tt/2) D^2} p = sum_j (-tt/2)^j p^{(2j)}/j!
    return sp.expand(sum((-tt / 2) ** j * sp.diff(p, x, 2 * j) / sp.factorial(j) for j in range(d // 2 + 1)))


def boxplus(p, q, d):
    # MSS: for p = sum_i x^{d-i} (-1)^i a_i, q = sum_j x^{d-j} (-1)^j b_j,
    # p boxplus_d q = sum_k x^{d-k} (-1)^k sum_{i+j=k} (d-i)! (d-j)!/(d! (d-k)!) a_i b_j
    P = sp.Poly(p, x)
    Q = sp.Poly(q, x)
    a = [(-1) ** i * P.coeff_monomial(x ** (d - i)) for i in range(d + 1)]
    b = [(-1) ** j * Q.coeff_monomial(x ** (d - j)) for j in range(d + 1)]
    r = 0
    for k in range(d + 1):
        ck = sum(sp.factorial(d - i) * sp.factorial(d - (k - i)) / (sp.factorial(d) * sp.factorial(d - k)) * a[i] * b[k - i] for i in range(k + 1))
        r += (-1) ** k * ck * x ** (d - k)
    return sp.expand(r)


d = 5
cs = sp.symbols('c0:6')
pg = sum(cs[i] * x ** i for i in range(d)) + x ** d
qt = heat(x ** d, t, d)
qs = heat(x ** d, s, d)
He5 = sp.expand(sp.hermite_prob(5, x)) if hasattr(sp, 'hermite_prob') else None
rec('q_t(x) - t^{5/2} He_5(x/sqrt t) (symbolic)', str(sp.simplify(qt - sp.expand(t ** sp.Rational(5, 2) * sp.hermite_prob(5, x / sp.sqrt(t))))) if He5 is not None else 'hermite_prob unavailable')
rec('p boxplus_5 q_t - e^{-(t/2) D^2} p (generic monic p, symbolic)', str(sp.simplify(boxplus(pg, qt, d) - heat(pg, t, d))))
rec('q_s boxplus_5 q_t - q_{s+t} (symbolic)', str(sp.simplify(boxplus(qs, qt, d) - heat(x ** d, s + t, d))))

# 2. the example
p = sp.expand((x + 3) * (x + 1) * (x - sp.Rational(1, 2)) * (x - 2) * ((x - 1) ** 2 + sp.Rational(16, 25)))
pt = heat(p, t, 6)
D = sp.Poly(sp.discriminant(pt, x), t)
rec('degree of the discriminant in t', D.degree())
mp.mp.dps = 50
roots_t = [r for r in sp.Poly(D, t).nroots(n=50, maxsteps=200) if abs(sp.im(r)) < 1e-40]
pos = sorted(float(sp.re(r)) for r in roots_t if sp.re(r) > 0)
tstar = min((sp.re(r) for r in roots_t if sp.re(r) > 0))
rec('least positive root t* of the discriminant (50 digits)', str(sp.N(tstar, 40)))
rec('real roots of the discriminant in (0, 2)', ['%.12g' % r for r in pos if r < 2])
tstar_q = sp.Rational(str(sp.N(tstar, 45)))
imax = {}
for lab, tv in (('0', 0), ('0.1', sp.Rational(1, 10)), ('0.2', sp.Rational(1, 5)), ('0.3', sp.Rational(3, 10)), ('t*', tstar_q), ('0.35', sp.Rational(7, 20))):
    rts = sp.Poly(pt.subs(t, tv), x).nroots(n=30, maxsteps=200)
    imax[lab] = ['%.6g' % max(abs(float(sp.im(r))) for r in rts), sorted(['%.6g' % float(sp.re(r)) for r in rts], key=float)]
rec('largest |Im| of the roots of e^{-(t/2) D^2} p, and the real parts, at t = 0, 0.1, 0.2, 0.3, t*, 0.35', imax)
rec('isolated pair: (x - 1)^2 + 0.64 reaches the axis at t = 0.64', str(sp.solve(sp.discriminant(heat(sp.expand((x - 1) ** 2 + sp.Rational(16, 25)), t, 2), x), t)))
# the paths of the six roots for the figure, t = 0, 0.0025, ..., 0.4 (numpy roots of the coefficients, sorted by real part)
cpoly = [sp.lambdify(t, c) for c in sp.Poly(pt, x).all_coeffs()]
paths = {'t': [], 're': [], 'im': []}
for tv in np.linspace(0, 0.4, 161):
    r = np.roots([float(c(tv)) for c in cpoly])
    r = r[np.lexsort((r.imag, r.real))]
    paths['t'].append(round(float(tv), 4))
    paths['re'].append([float('%.8g' % v) for v in r.real])
    paths['im'].append([float('%.8g' % v) for v in r.imag])
out['paths of the roots of e^{-(t/2) D^2} p'] = paths

# 3. the finite semicircle


def F(q):
    return 0.5 + (q * np.sqrt(4 - q * q)) / (4 * np.pi) + np.arcsin(q / 2) / np.pi


sc = {}
for dd in (20, 40, 80, 160, 320):
    off = np.sqrt(np.arange(1, dd))
    J = np.diag(off, 1) + np.diag(off, -1)
    r = np.sort(np.linalg.eigvalsh(J)) / np.sqrt(dd)
    qk = np.array([brentq(lambda q: F(q) - (k - 0.5) / dd, -2, 2, xtol=1e-15) for k in range(1, dd + 1)])
    dev = np.abs(r - qk)
    sc[dd] = ['%.6f' % dev.max(), '%.6f' % dev[1:-1].max(), '%.6f' % r[-1], '%.4f' % (dev.max() * dd ** (2 / 3))]
    if dd in (20, 160):
        out['roots and semicircle quantiles, d = %d' % dd] = {'roots': [float('%.8g' % v) for v in r], 'quantiles': [float('%.8g' % v) for v in qk]}
rec('roots of He_d(sqrt d x) against semicircle quantiles: largest deviation, without the two extremes, largest root, deviation times d^{2/3}', sc)
# check of the Jacobi matrix route against sympy's Hermite polynomial at d = 20
r20 = np.sort(np.linalg.eigvalsh(np.diag(np.sqrt(np.arange(1, 20)), 1) + np.diag(np.sqrt(np.arange(1, 20)), -1)))
rs = np.sort(np.array([float(v) for v in sp.Poly(sp.hermite_prob(20, x), x).nroots(n=30)]))
rec('roots of He_20: Jacobi matrix against sympy, largest difference', '%.2g' % np.abs(r20 - rs).max())
out['seconds'] = round(time.time() - t0, 1)
json.dump(out, open('wigfree_heat.json', 'w'), indent=1)
print('done', out['seconds'])
