# sqfree_gammanu.py: Gamma(nu) smoothing of the squarefree count, Section sec:sqfree (ch/sqfree.tex), Table tab:sqfree:nu.
# Added for the book in October 2026, from Nick Polson's note on the squarefree series (misc/reviews/nick-squarefree/),
# recomputed independently. With R_nu(h) = h^{nu+1} sum_{m>=1} m^nu e^{-mh}/Gamma(nu+1), the right Riemann sum with step h of the
# Gamma(nu+1) density, the mixing density of the Mertens representation is f_nu(t) = 2 nu (1 - R_nu(t))/t^2, so it is nonnegative
# exactly when R_nu <= 1. For h <= 2.5 the script uses the expansion
#   1 - R_nu(h) = -(h^{nu+1}/Gamma(nu+1)) sum_{k>=0} zeta(-nu-k) (-h)^k/k!     (|h| < 2 pi),
# which follows from Li_{-nu}(e^{-h}) = Gamma(nu+1) h^{-nu-1} + sum_k zeta(-nu-k)(-h)^k/k!, and for h > 2.5 the direct sum; the two are
# compared at h = 1, 2 and 2.5. For twenty shapes nu it records the least value of q_nu(h) = Gamma(nu+1)(1 - R_nu(h))/h^{nu+1}, whose
# limit at h = 0 is -zeta(-nu), and of 1 - R_nu(h) itself, on a grid of 1199 points (400 geometrically spaced in [10^{-4}, 10^{0.4}]
# and 799 equally spaced in (2.5, 60]) refined by golden section search, and the maximal runs of grid points where
# 1 - R_nu < 0. It also checks, at 40 digits: f_nu = -2 g_nu' with
# g_nu(t) = nu/t - sum_m Gamma(nu, mt)/Gamma(nu) (numerical derivative, nu = 0.5, 2, 3); the closed forms of R_1 and R_2; the expansion
# g(t) = 1/2 - t^3/360 + O(t^5) of the Gamma(2) case; and the Mellin transform int_0^oo v^s f_2(v) dv = -2 Gamma(s+2) zeta(s) at
# s = 1/2 and s -> 0 (mass one), by mpmath's quadrature. Run from any folder: python3 sqfree_gammanu.py; writes sqfree_gammanu.json.
import json
import time

import mpmath as mp

t0 = time.time()
out = {}


def rec(k, v):
    out[k] = v
    print(k, v, round(time.time() - t0, 1), flush=True)


def s(x, n=10):
    return mp.nstr(x, n)


mp.mp.dps = 40


ZT = {}


def zetas(nu, K=400):
    # zeta(-nu-k), k = 0, ..., K-1, computed once per shape
    if nu not in ZT:
        ZT[nu] = [mp.zeta(-nu - k) for k in range(K)]
    return ZT[nu]


def omr_series(nu, h):
    zt = zetas(nu)
    tot, term_fac, prev = mp.mpf(0), mp.mpf(1), mp.inf
    for k in range(len(zt)):
        if k > 0:
            term_fac *= -h / k
        term = zt[k] * term_fac
        tot += term
        # for integer nu every other zeta(-nu-k) vanishes, so stop on two consecutive small terms
        if k > 10 and abs(term) + abs(prev) < mp.mpf(10) ** -45 * (abs(tot) + mp.mpf(10) ** -60):
            return -h ** (nu + 1) / mp.gamma(nu + 1) * tot
        prev = term
    raise RuntimeError('series did not converge')


def omr_direct(nu, h):
    tot, m = mp.mpf(0), 1
    while True:
        term = mp.mpf(m) ** nu * mp.exp(-m * h)
        tot += term
        if m * h > 40 and term < mp.mpf(10) ** -50 * tot:
            break
        m += 1
    return 1 - h ** (nu + 1) * tot / mp.gamma(nu + 1)


def omr(nu, h):
    return omr_series(nu, h) if h <= mp.mpf('2.5') else omr_direct(nu, h)


chk = {}
for nu in (mp.mpf('0.5'), mp.mpf(2), mp.mpf('2.5'), mp.mpf(8)):
    chk[s(nu, 3)] = [s(omr_series(nu, h) - omr_direct(nu, h), 3) for h in (mp.mpf(1), mp.mpf(2), mp.mpf('2.5'))]
rec('series minus direct sum at h = 1, 2, 2.5', chk)
h = mp.mpf('1.3')
R1 = (h / 2) ** 2 / mp.sinh(h / 2) ** 2
R2 = (h / 2) ** 3 * mp.cosh(h / 2) / mp.sinh(h / 2) ** 3
rec('closed forms of R_1 and R_2 at h = 1.3 against the sums', [s(1 - R1 - omr(mp.mpf(1), h), 3), s(1 - R2 - omr(mp.mpf(2), h), 3)])


def g_nu(nu, t):
    return nu / t - mp.nsum(lambda m: mp.gammainc(nu, m * t, mp.inf, regularized=True), [1, mp.inf])


der = {}
for nu in (mp.mpf('0.5'), mp.mpf(2), mp.mpf(3)):
    der[s(nu, 3)] = [s(-2 * mp.diff(lambda t: g_nu(nu, t), t) - 2 * nu * omr(nu, t) / t ** 2, 3) for t in (mp.mpf('0.7'), mp.mpf(3))]
rec('-2 g_nu\'(t) - 2 nu (1 - R_nu(t))/t^2 at t = 0.7, 3', der)


def g2(t):
    return 2 / t - 1 / mp.expm1(t) - t / (4 * mp.sinh(t / 2) ** 2)


rec('(g(t) - 1/2 + t^3/360)/t^5 at t = 0.01, 0.001', [s((g2(mp.mpf('0.01')) - mp.mpf(1) / 2 + mp.mpf('0.01') ** 3 / 360) / mp.mpf('0.01') ** 5, 8),
                                                      s((g2(mp.mpf('0.001')) - mp.mpf(1) / 2 + mp.mpf('0.001') ** 3 / 360) / mp.mpf('0.001') ** 5, 8)])


def f2(v):
    if v < mp.mpf('1e-6'):
        return v ** 2 / 60
    return 4 / v ** 2 - v * mp.cosh(v / 2) / (2 * mp.sinh(v / 2) ** 3)


mass = mp.quad(f2, [0, 1, 10, 100, mp.inf])
m12 = mp.quad(lambda v: mp.sqrt(v) * f2(v), [0, 1, 10, 100, mp.inf])
rec('Gamma(2) mixing density: mass - 1; int v^(1/2) f_2 - (-2 Gamma(5/2) zeta(1/2)); -2 Gamma(5/2) zeta(1/2)',
    [s(mass - 1, 3), s(m12 + 2 * mp.gamma(mp.mpf('2.5')) * mp.zeta(mp.mpf('0.5')), 3), s(-2 * mp.gamma(mp.mpf('2.5')) * mp.zeta(mp.mpf('0.5')), 12)])
rec('Gamma(2): min of f_2 on the grid t = 0.01 j, j = 1, ..., 4000', s(min(f2(mp.mpf(j) / 100) for j in range(1, 4001)), 6))

grid = [mp.mpf(10) ** (-4 + 4.4 * mp.mpf(j) / 399) for j in range(400)] + [mp.mpf('2.5') + (60 - mp.mpf('2.5')) * mp.mpf(j) / 799 for j in range(1, 800)]


def golden(f, a, b, it=80):
    gr = (mp.sqrt(5) - 1) / 2
    c, d = b - gr * (b - a), a + gr * (b - a)
    fc, fd = f(c), f(d)
    for _ in range(it):
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - gr * (b - a)
            fc = f(c)
        else:
            a, c, fc = c, d, fd
            d = a + gr * (b - a)
            fd = f(d)
    return (a + b) / 2


rows = []
for nu in ('0.25', '0.5', '0.75', '1', '1.25', '1.5', '1.75', '1.9', '1.99', '2', '2.01', '2.1', '2.5', '3', '3.5', '4', '5', '6', '7', '8'):
    nu = mp.mpf(nu)
    vals = [omr(nu, hh) for hh in grid]
    qv = [mp.gamma(nu + 1) * v / hh ** (nu + 1) for v, hh in zip(vals, grid)]
    jq = min(range(len(grid)), key=lambda j: qv[j])
    jr = min(range(len(grid)), key=lambda j: vals[j])
    hq = golden(lambda t: mp.gamma(nu + 1) * omr(nu, t) / t ** (nu + 1), grid[max(jq - 1, 0)], grid[min(jq + 1, len(grid) - 1)])
    hr = golden(lambda t: omr(nu, t), grid[max(jr - 1, 0)], grid[min(jr + 1, len(grid) - 1)])
    neg = [hh for v, hh in zip(vals, grid) if v < 0]
    runs, cur = [], None
    for v, hh in zip(vals, grid):
        if v < 0:
            cur = [hh, hh] if cur is None else [cur[0], hh]
        elif cur is not None:
            runs.append(cur)
            cur = None
    if cur is not None:
        runs.append(cur)
    note_range = [v for v, hh in zip(vals, grid) if mp.mpf('0.02') <= hh <= 15]
    rows.append({'nu': s(nu, 4), '-zeta(-nu)': s(-mp.zeta(-nu), 6), 'min q_nu': s(mp.gamma(nu + 1) * omr(nu, hq) / hq ** (nu + 1), 6), 'at h': s(hq, 5),
                 'min (1 - R_nu)': s(omr(nu, hr), 6), 'at h ': s(hr, 5),
                 'grid points with 1 - R_nu < 0': len(neg), 'negative on [h_min, h_max] of the grid': [s(neg[0], 4), s(neg[-1], 4)] if neg else None,
                 'maximal runs of negative grid points [first, last]': [[s(a, 4), s(b, 4)] for a, b in runs],
                 'min (1 - R_nu) on the note\'s range [0.02, 15]': s(min(note_range), 6)})
    print(rows[-1], flush=True)
rec('shapes', rows)
out['seconds'] = time.time() - t0
json.dump(out, open('sqfree_gammanu.json', 'w'), indent=1)
print('done', round(time.time() - t0, 1))
