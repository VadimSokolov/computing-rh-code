# Added for the book (not part of the authors' package): recomputes the density v_alpha/pi of almost.py (the first 1700
# ordinates on the line plus the fifteen hypothetical quadruples at real part 0.75) at the basepoints 0.6 and 0.7 of Figure
# ag:fig:almost on a fine grid over the window [284, 300] of its middle and right panels, and near every hypothetical
# ordinate in the range [0, 2100] of almost.py, and lists every local minimum below zero with its height and depth, for
# the sentence of Chapter ch:almost that quotes the dip of the density near the height 292.4.
#
# Model: exactly that of almost.py, P(a, x) = a/(pi (a^2 + x^2)) and
#   r(theta) = sum_gamma [P(alpha - 1/2, theta - gamma) + P(alpha - 1/2, theta + gamma)]
#            + sum_j [P(alpha - b, theta - t_j) + P(alpha - (1 - b), theta - t_j) + (the same at theta + t_j)],
# with b = 0.75 and t_j = 20 j^(5/3) < 2000; r0 is the first sum alone (zeros on the line only).
# Where negative values can occur: for 1/2 < alpha < b write u = b - alpha and v = alpha - (1 - b), so v > u > 0. Then
# P(v, x) - P(u, x) > 0 exactly when x^2 > u v, so outside the intervals |theta - t_j| < sqrt(u v) (0.15 at alpha = 0.7,
# 0.229 at alpha = 0.6) every group of hypothetical kernels is positive, and so is r, since r0 > 0 and theta + t_j >= 20.
# Each such interval is scanned with step 1e-5, the window of the figure with step 1e-4; every grid minimum is refined
# by a root of the analytic derivative (brentq) and every negative one is checked in 30 digit arithmetic (mpmath).
# Run from code/ after ./setup.sh (it reads the zero lists and data/almost.json); writes almost_dip.json.
import json, math
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad
import mpmath as mp

g = np.concatenate([np.load(f) for f in ['g_1_400.npy', 'g_401_1000.npy', 'g_1001_1700.npy']])
b = 0.75
tj = []
k = 1
while True:                                   # the heights of almost.py: count grows like t^0.6
    t = (k / 1.0) ** (1 / 0.6) * 20
    if t > 2000:
        break
    tj.append(t)
    k += 1
CH = 4000


def P(a, x):
    return a / (np.pi * (a * a + x * x))


def dP(a, x):
    return -2 * a * x / (np.pi * (a * a + x * x) ** 2)


def r0(th, alpha):
    """zeros on the line only, on an array of heights"""
    th = np.asarray(th, dtype=float)
    out = np.empty_like(th)
    a = alpha - 0.5
    for i in range(0, len(th), CH):
        x = th[i:i + CH, None]
        out[i:i + CH] = (P(a, x - g[None, :]) + P(a, x + g[None, :])).sum(axis=1)
    return out


def off(th, alpha):
    th = np.asarray(th, dtype=float)
    s = np.zeros_like(th)
    for t in tj:
        s += P(alpha - b, th - t) + P(alpha - (1 - b), th - t) + P(alpha - b, th + t) + P(alpha - (1 - b), th + t)
    return s


def r(th, alpha):
    return r0(th, alpha) + off(th, alpha)


def r1(x, alpha):
    return float(r(np.array([x]), alpha)[0])


def dr1(x, alpha):
    a = alpha - 0.5
    d = float(np.sum(dP(a, x - g) + dP(a, x + g)))
    for t in tj:
        d += dP(alpha - b, x - t) + dP(alpha - (1 - b), x - t) + dP(alpha - b, x + t) + dP(alpha - (1 - b), x + t)
    return d


mp.mp.dps = 30
gm = [mp.mpf(float(v)) for v in g]


def r_mp(x, alpha):
    x = mp.mpf(x); a = mp.mpf(alpha) - mp.mpf('0.5')
    Pm = lambda w, y: w / (mp.pi * (w * w + y * y))
    s = mp.fsum(Pm(a, x - v) + Pm(a, x + v) for v in gm)
    for t in tj:
        t = mp.mpf(t)
        for w in (mp.mpf(alpha) - mp.mpf(b), mp.mpf(alpha) - (1 - mp.mpf(b))):
            s += Pm(w, x - t) + Pm(w, x + t)
    return s


def negative_minima(th, val, alpha):
    """every interior grid minimum of val, refined; returns those below zero with their negative interval and mass"""
    res = []
    idx = np.where((val[1:-1] < val[:-2]) & (val[1:-1] <= val[2:]))[0] + 1
    for i in idx:
        lo, hi = th[i - 1], th[i + 1]
        if dr1(lo, alpha) < 0 < dr1(hi, alpha):
            x = brentq(dr1, lo, hi, args=(alpha,), xtol=1e-13)
        else:
            x = th[i]
        v = r1(x, alpha)
        if v >= 0:
            continue
        j = i
        while j > 0 and val[j] < 0:
            j -= 1
        left = brentq(r1, th[j], th[j + 1], args=(alpha,), xtol=1e-12) if val[j] >= 0 else None
        j = i
        while j < len(val) - 1 and val[j] < 0:
            j += 1
        right = brentq(r1, th[j - 1], th[j], args=(alpha,), xtol=1e-12) if val[j] >= 0 else None
        mass = None
        if left is not None and right is not None:
            mass = quad(lambda y: -r1(y, alpha), left, right, points=[x], limit=200, epsabs=1e-12, epsrel=1e-10)[0]
        vmp = r_mp(x, alpha)
        near = min(tj, key=lambda t: abs(t - x))
        res.append(dict(theta=x, value=v, value_mp30=float(vmp), float_error=float(abs(vmp - v)),
                        background=float(r0(np.array([x]), alpha)[0]),
                        second_derivative=(dr1(x + 1e-6, alpha) - dr1(x - 1e-6, alpha)) / 2e-6,
                        negative_interval=[left, right], negative_mass=mass,
                        nearest_hypothetical_ordinate=near, offset_from_it=x - near))
    return res


out = dict(model='almost.py: first 1700 ordinates on the line, hypothetical quadruples 0.75 +- i t_j, 0.25 +- i t_j',
           tj=tj, n_ordinates=int(len(g)), largest_ordinate=float(g[-1]))

# 1. the recomputation agrees with the archived curves of almost.py (almost.json keeps every tenth point, step 0.5)
A = json.load(open('almost.json'))
th_a = np.linspace(0, 2100, 42001)
chk = {}
for alpha in (0.6, 0.7, 0.8):
    C = A['curve_%s' % alpha]
    thj = np.array(C['th'])
    rr, rr0 = r(thj, alpha), r0(thj, alpha)
    ra = np.array(C['r']); i = int(np.argmin(ra))
    chk[str(alpha)] = dict(max_abs_diff_r=float(np.max(np.abs(rr - ra))), max_abs_diff_r0=float(np.max(np.abs(rr0 - np.array(C['r0'])))),
                           same_grid=bool(np.allclose(thj, th_a[::10], atol=1e-12, rtol=0)),
                           archived_min=dict(value=float(ra[i]), theta=float(thj[i])),
                           archived_near_292=[dict(theta=float(thj[m]), value=float(ra[m])) for m in np.where((thj > 290.9) & (thj < 294.1))[0]])
    print('almost.json alpha', alpha, chk[str(alpha)]['max_abs_diff_r'], chk[str(alpha)]['max_abs_diff_r0'], 'archived min', chk[str(alpha)]['archived_min'], flush=True)
out['check_against_almost_json'] = chk

# 2. the window of the figure (middle and right panels): [284, 300]
W0, W1 = 284.0, 300.0
thw = np.linspace(W0, W1, 160001)                # step 1e-4
thf = np.linspace(W0, W1, 3201)                  # step 0.005, the grid of the published figure (vertices at multiples of 0.005)
ths = np.linspace(W0, W1, 8001)                  # step 0.002, stored for plotting
win = dict(window=[W0, W1], grid_step=1e-4, curves_step=0.002)
for alpha in (0.6, 0.7):
    bg = r0(thw, alpha); val = bg + off(thw, alpha)
    neg = negative_minima(thw, val, alpha)
    ff = r(thf, alpha); f0 = r0(thf, alpha); i = int(np.argmin(ff)); i0 = int(np.argmin(f0))
    m = int(np.argmin(np.abs(th_a - 292.4)))
    win[str(alpha)] = dict(
        negative_minima=neg,
        grid_min=dict(value=float(val.min()), theta=float(thw[np.argmin(val)])),
        background_min=dict(value=float(bg.min()), theta=float(thw[np.argmin(bg)])),
        figure_grid_min=dict(value=float(ff[i]), theta=float(thf[i])),
        figure_grid_background_min=dict(value=float(f0[i0]), theta=float(thf[i0])),
        almost_py_grid_value_at_292_4=dict(theta=float(th_a[m]), value=r1(th_a[m], alpha)),
        curves=dict(theta0=W0, step=0.002, with_hypothetical=[round(float(v), 6) for v in r(ths, alpha)],
                    line_only=[round(float(v), 6) for v in r0(ths, alpha)]))
    for q in neg:
        print('window alpha %.1f: local minimum %.6f at theta %.6f (t_j %.6f), background %.4f, negative on [%.5f, %.5f], mass %.5f, mp30 %.3e'
              % (alpha, q['value'], q['theta'], q['nearest_hypothetical_ordinate'], q['background'], q['negative_interval'][0], q['negative_interval'][1], q['negative_mass'], q['float_error']), flush=True)
    print('  figure grid (step 0.005) min %.4f at %.3f; background min %.4f at %.3f; almost.py grid value at %.2f: %.4f'
          % (ff[i], thf[i], f0[i0], thf[i0], th_a[m], win[str(alpha)]['almost_py_grid_value_at_292_4']['value']), flush=True)
out['figure_window'] = win

# 3. the whole range of almost.py: negative values only within sqrt(u v) of a hypothetical ordinate
full = {}
for alpha in (0.6, 0.7, 0.8):
    u, v = b - alpha, alpha - (1 - b)
    if u <= 0:
        rg = r(th_a, alpha)
        full[str(alpha)] = dict(note='alpha > b: every kernel is positive', min_on_almost_py_grid=float(rg.min()),
                                theta_of_min=float(th_a[np.argmin(rg)]))
        print('alpha', alpha, 'min on the grid of almost.py', rg.min(), flush=True)
        continue
    delta = math.sqrt(u * v)
    lst, wmin = [], []
    for t in tj:
        th = np.linspace(t - 0.25, t + 0.25, 50001)                # step 1e-5
        val = r(th, alpha)
        lst += negative_minima(th, val, alpha)
        wmin.append(dict(t=t, grid_min=float(val.min()), theta=float(th[np.argmin(val)])))
    far = np.min([np.abs(th_a - t) for t in tj], axis=0) > delta
    rg = r(th_a, alpha)
    # exact negative mass on [0, T] (each negative interval once), for comparison with the Riemann sums of almost.py
    iv = {}
    for q in lst:
        iv[(round(q['negative_interval'][0], 9), round(q['negative_interval'][1], 9))] = q['negative_mass']
    masses = {str(T): sum(m for (l, h), m in iv.items() if h <= T) for T in (250, 1000, 2000)}
    full[str(alpha)] = dict(delta=delta, negative_minima=lst, window_minima=wmin, negative_mass_on_0_T=masses,
                            min_on_almost_py_grid_away_from_intervals=float(rg[far].min()))
    print('alpha', alpha, 'grid minimum near each t_j:', ' '.join('%.1f:%.4f' % (w['t'], w['grid_min']) for w in wmin), flush=True)
    print('alpha', alpha, 'exact negative mass on [0, T]:', masses, flush=True)
    for q in lst:
        print('alpha %.1f t_j %8.3f: minimum %.6f at %.6f (offset %+.2e), background %.4f, width of negative interval %.5f, mass %.5f'
              % (alpha, q['nearest_hypothetical_ordinate'], q['value'], q['theta'], q['offset_from_it'], q['background'],
                 q['negative_interval'][1] - q['negative_interval'][0], q['negative_mass']), flush=True)
    print('alpha', alpha, 'delta', delta, 'min of r on the grid of almost.py away from the intervals', full[str(alpha)]['min_on_almost_py_grid_away_from_intervals'], flush=True)
out['full_range'] = full
json.dump(out, open('almost_dip.json', 'w'))
print('ok')
