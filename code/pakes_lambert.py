# pakes_lambert.py: the Lambert law of Section pk:sec:lambert (ch/pakes.tex): Proposition pk:prop:bv, Theorem pk:thm:lambert,
# Corollaries pk:cor:levy and pk:cor:nuttall, and the two bounds used in the proof of the theorem.
# Added for the book in October 2026, for the merge of the first author's chapter on Pakes, the Lambert law and the Thorin
# condition (misc/reviews/nick-pakes/); adapted from the triage checks B_lambert.py and B_levy.py of that folder.
# Notation as in the section: W is the principal branch of Lambert's function (written \mathcal W in the book), the solution of
# W(u) e^{W(u)} = u analytic off (-infinity, -1/e]; r(b) = b e^{-b cot b}/sin b for 0 < b < pi; U_W is the law of r(B) with B
# uniform on (0, pi); the theorem says W(u) = (1/pi) int_0^pi log(1 + u/r(b)) db.
# Checks:
#   1. boundary values W(-r(b) + i0) = -b cot b + i b, the identity b sin^2 b (log r)'(b) = (b - sin b cos b)^2 + sin^4 b,
#      the limits r(0+) = 1/e (with r(b) - 1/e ~ b^2/(2e)) and r(pi-) = infinity;
#   2. the Thorin integral against mpmath's W at 40 digits, at seventeen real and complex points of the cut plane, and
#      W(2.7) at 60 digits by Newton's method on w e^w = 27/10 and by the integral (the input must be the rational 27/10:
#      the binary double nearest 2.7 changes W from the 17th decimal on), with an Arb value if python-flint provides one;
#   3. Nuttall's moments (1/pi) int_0^pi r(b)^{-n} db = n^n/n! for n <= 12, and nu^nu/Gamma(nu+1) at four real nu;
#   4. the Levy sides of Corollary pk:cor:levy in double precision (scipy): W(u) = int_0^infty (1 - e^{-ux}) k(x) dx/x with
#      k(x) = (1/pi) int_0^pi e^{-x r(b)} db, and W(w^2) = int_0^infty 4 sin^2(wx/2) k_s(x) dx/x with
#      k_s(x) = (1/pi) int_0^pi e^{-x sqrt(r(b))} db; near b = pi, r(b) overflows and exp(-x r) = 0 is exact in double precision
#      (at 25 digits mpmath evaluates exp(-x r(b)) with r(b) near exp(1e8) and does not finish);
#   5. the growth facts behind the proof: W(u) - log u + log log u -> 0, the divergence of int log z U_W(dz) (partial
#      integrals up to b = pi - eps exceed log(1/eps) by a bounded amount), and the tail U_W(z > y) ~ 1/log y;
#   6. the lower bound Re W(u) >= log R - log(log R + pi) on circles |u| = R >= e (both sides of the cut included);
#   7. the density of U_W near 1/e against sqrt(e/2) pi^{-1} (z - 1/e)^{-1/2}.
# Run on Hopper (about ten seconds): python3 pakes_lambert.py; writes pakes_lambert.json in the current folder.
import json
import time
import warnings

import mpmath as mp
import numpy as np
from scipy.integrate import quad
from scipy.special import lambertw

warnings.filterwarnings('ignore')
mp.mp.dps = 40
PI = mp.pi
E = mp.e
out = {}
t_start = time.time()


def r_of_b(b):
    b = mp.mpf(b)
    return b * mp.exp(-b * mp.cot(b)) / mp.sin(b)


def t_of_b(b):
    """1/r(b), the atom of the image of U_W under z -> 1/z."""
    b = mp.mpf(b)
    return mp.sin(b) * mp.exp(b * mp.cot(b)) / b


def logr_of_b(b):
    b = mp.mpf(b)
    return mp.log(b) - b * mp.cot(b) - mp.log(mp.sin(b))


def b_of_r(y):
    """Inverse of r on (0, pi), by bisection (r is increasing)."""
    y = mp.mpf(y)
    lo, hi = mp.mpf('1e-30'), PI - mp.mpf('1e-6')
    for _ in range(200):
        mid = (lo + hi) / 2
        if r_of_b(mid) < y:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


# ---------------------------------------------------------------- 1. boundary values (Proposition pk:prop:bv)
res = []
for b in ['0.05', '0.5', '1.0', '1.3', '2.0', '2.5', '3.0', '3.1']:
    b = mp.mpf(b)
    y = r_of_b(b)
    target = -b * mp.cot(b) + 1j * b
    w_above = mp.lambertw(mp.mpc(-y, mp.mpf('1e-60')))
    w_below = mp.lambertw(mp.mpc(-y, -mp.mpf('1e-60')))
    resid = target * mp.exp(target) + y
    res.append({'b': str(b), 'r(b)': mp.nstr(y, 20),
                'abs(W(-r+i0) - (-b cot b + ib))': mp.nstr(abs(w_above - target), 3),
                'abs(W(-r-i0) - conj)': mp.nstr(abs(w_below - mp.conj(target)), 3),
                'abs(w e^w + r) at w = -b cot b + ib': mp.nstr(abs(resid), 3)})
out['boundary_values'] = res
mono = []
for b in ['0.01', '0.7', '1.5', '2.9', '3.14']:
    b = mp.mpf(b)
    dlog = mp.diff(lambda x: mp.log(r_of_b(x)), b)
    lhs = b * mp.sin(b) ** 2 * dlog
    rhs = (b - mp.sin(b) * mp.cos(b)) ** 2 + mp.sin(b) ** 4
    formula = 1 / b - 2 * mp.cot(b) + b / mp.sin(b) ** 2
    mono.append({'b': str(b), 'b sin^2 b (log r)\'': mp.nstr(lhs, 15), '(b - sin b cos b)^2 + sin^4 b': mp.nstr(rhs, 15),
                 'abs((log r)\' - (1/b - 2 cot b + b/sin^2 b))': mp.nstr(abs(dlog - formula), 3)})
out['monotonicity'] = mono
out['limits'] = {'r(1e-12) - 1/e': mp.nstr(r_of_b(mp.mpf('1e-12')) - 1 / E, 5),
                 'r(pi - 1e-2)': mp.nstr(r_of_b(PI - mp.mpf('1e-2')), 8),
                 'r(pi - 1e-3)': mp.nstr(r_of_b(PI - mp.mpf('1e-3')), 8)}
out['(r(b) - 1/e)/b^2 at b = 1e-2, 1e-4, 1e-6'] = [mp.nstr((r_of_b(mp.mpf(b)) - 1 / E) / mp.mpf(b) ** 2, 12)
                                                     for b in ['1e-2', '1e-4', '1e-6']]
out['1/(2e)'] = mp.nstr(1 / (2 * E), 12)
print('1. boundary values done', round(time.time() - t_start, 1), 's', flush=True)


# ---------------------------------------------------------------- 2. the Thorin representation (Theorem pk:thm:lambert)
def thorin_W(u):
    """(1/pi) int_0^pi log(1 + u t(b)) db, with breakpoints at the logarithmic transitions."""
    u = mp.mpc(u)
    pts = [mp.mpf(0), PI / 4, PI / 2, 3 * PI / 4]
    if abs(u) > 1 / E:  # where |u| t(b) = 1, which matters next to the cut
        b0 = b_of_r(abs(u))
        for d in ['0', '1e-3', '1e-2', '0.1']:
            for sgn in (-1, 1):
                p = b0 + sgn * mp.mpf(d)
                if 0 < p < PI:
                    pts.append(p)
    if abs(u) > 10:  # transition near pi for large |u|: pi - b ~ pi/log|u|
        L = mp.log(abs(u))
        for c in ['0.5', '1', '2', '4']:
            p = PI - mp.mpf(c) * PI / L
            if 0 < p < PI:
                pts.append(p)
    pts.append(PI)
    pts = sorted(set(pts))

    def f(b):
        x = u * t_of_b(b)
        # near b = pi, t(b) underflows to astronomically small values; there log(1 + x) = x to |x|^2/2,
        # and forming 1 + x exactly would exhaust memory
        if abs(x) < mp.mpf(10) ** (-2 * mp.mp.dps):
            return x
        return mp.log(1 + x)
    return mp.quad(f, pts) / PI


points = ['-0.3678', '-0.36', '-0.2', '0.5', '1', '2.7', '10', '1e3', '1e6', '1e12',
          '-2+1j', '3j', '-5+0.001j', '-100+1j', '1e4j', '-0.5+0.1j', '-1e6+1e3j']
rep = []
worst = mp.mpf(0)
for p in points:
    u = mp.mpc(complex(p)) if 'j' in p else mp.mpf(p)
    lhs = mp.lambertw(u)
    rhs = thorin_W(u)
    d = abs(lhs - rhs)
    worst = max(worst, d)
    rep.append({'u': p, 'W(u)': mp.nstr(lhs, 25), 'Thorin integral': mp.nstr(rhs, 25), 'abs diff': mp.nstr(d, 3)})
    print('u', p, 'diff', mp.nstr(d, 3), flush=True)
out['thorin_representation'] = rep
out['thorin_representation_max_abs_diff'] = mp.nstr(worst, 3)
with mp.workdps(60):
    u = mp.mpf(27) / 10
    w = mp.mpf(1)
    for _ in range(12):
        w = w - (w * mp.exp(w) - u) / (mp.exp(w) * (w + 1))
    out['W(27/10), Newton, 60 digits'] = mp.nstr(w, 50)
    out['W(27/10), Thorin integral, 60 digits'] = mp.nstr(thorin_W(u).real, 50)
    ud = mp.mpf(2.7)  # the binary double nearest 2.7
    out['double nearest 2.7'] = mp.nstr(ud, 55)
    out['W(double nearest 2.7)'] = mp.nstr(mp.lambertw(ud), 40)
    out['W(double) - W(27/10)'] = mp.nstr(mp.lambertw(ud) - w, 6)
    out["(double - 27/10) W'(27/10)"] = mp.nstr((ud - u) * w / (u * (1 + w)), 6)
try:
    import flint
    flint.ctx.prec = 200
    out['W(27/10), Arb'] = (flint.arb(27) / 10).lambertw().str(40, radius=True)
except Exception as exc:  # python-flint missing, or no lambertw method
    out['W(27/10), Arb'] = 'not available: %r' % (exc,)
with mp.workdps(30):
    b = mp.mpf('1.3')
    yv = r_of_b(b)
    out['abs(W(-r(1.3)+i0) - (-1.3 cot 1.3 + 1.3i)), 30 digits'] = mp.nstr(
        abs(mp.lambertw(mp.mpc(-yv, mp.mpf('1e-40'))) - (-b * mp.cot(b) + 1j * b)), 3)
print('2. Thorin representation done', round(time.time() - t_start, 1), 's', flush=True)

# ---------------------------------------------------------------- 3. Nuttall's moments (Corollary pk:cor:nuttall)
nut = []
for n in range(0, 13):
    m = mp.quad(lambda b: t_of_b(b) ** n, [0, PI / 4, PI / 2, 3 * PI / 4, PI]) / PI
    exact = mp.mpf(n) ** n / mp.factorial(n) if n > 0 else mp.mpf(1)
    nut.append({'n': n, 'moment': mp.nstr(m, 20), 'n^n/n!': mp.nstr(exact, 20), 'abs diff': mp.nstr(abs(m - exact), 3)})
out['nuttall_integer'] = nut
nutr = []
for nu in ['0.5', '1.5', '2.7', '7.25']:
    nu = mp.mpf(nu)
    m = mp.quad(lambda b: t_of_b(b) ** nu, [0, PI / 4, PI / 2, 3 * PI / 4, PI]) / PI
    exact = nu ** nu / mp.gamma(nu + 1)
    nutr.append({'nu': str(nu), 'moment': mp.nstr(m, 20), 'nu^nu/Gamma(nu+1)': mp.nstr(exact, 20),
                 'abs diff': mp.nstr(abs(m - exact), 3)})
out['nuttall_real'] = nutr
out['n^n/n!, n = 0..6'] = [mp.nstr(mp.mpf(n) ** n / mp.factorial(n) if n else 1, 6) for n in range(7)]
print('3. Nuttall done', round(time.time() - t_start, 1), 's', flush=True)


# ---------------------------------------------------------------- 4. the Levy sides (Corollary pk:cor:levy), double precision
def logr_np(b):
    return np.log(b) - b / np.tan(b) - np.log(np.sin(b))


def k_np(x, power=1.0):
    """(1/pi) int_0^pi exp(-x r(b)^power) db, with r^power = exp(power log r)."""
    f = lambda b: np.exp(-x * np.exp(power * logr_np(b))) if b > 0 else np.exp(-x * np.exp(-power))
    val, err = quad(f, 0.0, np.pi, points=[1e-3, 0.3, 1.0, 2.0, 2.8], limit=400, epsabs=1e-15, epsrel=1e-13)
    return val / np.pi


def levy_W(u):
    f = lambda x: (-np.expm1(-u * x)) * k_np(x) / x
    pieces = [0, 1e-6, 1e-3, 0.1, 1, 5, 20, 60, 150, 400]
    return sum(quad(f, a, b, limit=400, epsabs=1e-14, epsrel=1e-12)[0] for a, b in zip(pieces[:-1], pieces[1:]))


def levy_sym(w):
    f = lambda x: 4 * np.sin(w * x / 2) ** 2 * k_np(x, 0.5) / x
    edges = list(np.arange(0.0, 120.0 + 1e-9, np.pi / w))
    return sum(quad(f, a, b, limit=200, epsabs=1e-15, epsrel=1e-12)[0] for a, b in zip(edges[:-1], edges[1:]))


lev = []
for u in [1.0, 2.7, 10.0]:
    v1, v0 = float(levy_W(u)), float(lambertw(u).real)
    lev.append({'u': u, 'Levy integral': repr(v1), 'W(u)': repr(v0), 'abs diff': '%.2e' % abs(v1 - v0)})
out['levy_positive'] = lev
sym = []
for w in [1.0, 2.0]:
    v1, v0 = float(levy_sym(w)), float(lambertw(w * w).real)
    sym.append({'w': w, 'symmetric Levy integral': repr(v1), 'W(w^2)': repr(v0), 'abs diff': '%.2e' % abs(v1 - v0)})
out['levy_symmetric'] = sym
out['k_s(x) near 0 (k_s(0+) = 1, the mass of U_W)'] = [
    {'x': x, 'k_s(x)': repr(float(k_np(x, 0.5))), '1 - 1/(2 log(1/x))': repr(float(1 - 1 / (2 * np.log(1 / x))))}
    for x in [1e-3, 1e-9, 1e-30]]
out['k(x) against sqrt(e/2) e^{-x/e}/sqrt(pi x)'] = [
    {'x': x, 'k(x)': '%.10e' % k_np(x), 'asymptote': '%.10e' % (np.sqrt(np.e / 2) * np.exp(-x / np.e) / np.sqrt(np.pi * x))}
    for x in [10.0, 100.0, 1000.0]]
print('4. Levy sides done', round(time.time() - t_start, 1), 's', flush=True)

# ---------------------------------------------------------------- 5. growth facts
asym = []
for u in ['1e2', '1e4', '1e8', '1e16', '1e32', '1e64', '1e128']:
    u = mp.mpf(u)
    w = mp.lambertw(u)
    asym.append({'u': mp.nstr(u, 3), 'W - log u': mp.nstr(w - mp.log(u), 12), '-log log u': mp.nstr(-mp.log(mp.log(u)), 12),
                 'W - log u + log log u': mp.nstr(w - mp.log(u) + mp.log(mp.log(u)), 12)})
out['W_minus_log'] = asym
div = []
for eps in ['1e-2', '1e-4', '1e-8', '1e-16']:
    eps = mp.mpf(eps)
    mid = PI - 10 * eps if PI - 10 * eps > PI / 2 else PI / 2 + eps
    val = mp.quad(logr_of_b, [0, PI / 2, mid, PI - eps]) / PI
    div.append({'eps': mp.nstr(eps, 2), '(1/pi) int_0^{pi-eps} log r(b) db': mp.nstr(val, 10), 'log(1/eps)': mp.nstr(mp.log(1 / eps), 10)})
out['int_log_z_partial'] = div
tail = []
for yv in ['1e2', '1e4', '1e8', '1e16']:
    yv = mp.mpf(yv)
    tl = 1 - b_of_r(yv) / PI
    tail.append({'y': mp.nstr(yv, 3), 'U_W(z > y)': mp.nstr(tl, 10), '1/log y': mp.nstr(1 / mp.log(yv), 10)})
out['U_W_tail'] = tail
print('5. growth facts done', round(time.time() - t_start, 1), 's', flush=True)

# ---------------------------------------------------------------- 6. the lower bound used with Borel and Caratheodory
bc = []
for R in ['3', '10', '100', '1e4', '1e8']:
    R = mp.mpf(R)
    worst, worst_phi, maxratio = mp.inf, None, mp.mpf(0)
    for j in range(0, 401):
        phi = -PI + 2 * PI * j / 400
        u = R * mp.expj(phi)
        if j == 0:  # the two sides of the cut: upper side at j = 0, lower side at j = 400
            u = mp.mpc(-R, mp.mpf('1e-30'))
        elif j == 400:
            u = mp.mpc(-R, -mp.mpf('1e-30'))
        w = mp.lambertw(u)
        margin = w.real - (mp.log(R) - mp.log(mp.log(R) + PI))
        if margin < worst:
            worst, worst_phi = margin, phi
        maxratio = max(maxratio, abs(w) / (mp.log(R) + PI))
    bc.append({'R': mp.nstr(R, 3), 'min over |u| = R of Re W - (log R - log(log R + pi))': mp.nstr(worst, 8),
               'at arg u': mp.nstr(worst_phi, 6), 'max |W| / (log R + pi)': mp.nstr(maxratio, 6)})
out['lower_bound_on_circles'] = bc
print('6. lower bound done', round(time.time() - t_start, 1), 's', flush=True)

# ---------------------------------------------------------------- 7. the density of U_W near 1/e
dens = []
for b in ['1e-2', '1e-3', '1e-4']:
    b = mp.mpf(b)
    z = r_of_b(b)
    n = 1 / (PI * mp.diff(r_of_b, b))
    dens.append({'z - 1/e': mp.nstr(z - 1 / E, 6), 'density': mp.nstr(n, 10),
                 'sqrt(e/2)/(pi sqrt(z - 1/e))': mp.nstr(mp.sqrt(E / 2) / (PI * mp.sqrt(z - 1 / E)), 10)})
out['U_W_density_near_1_over_e'] = dens

out['seconds'] = round(time.time() - t_start, 1)
json.dump(out, open('pakes_lambert.json', 'w'), indent=1, default=str)
print(json.dumps(out, indent=1, default=str))
