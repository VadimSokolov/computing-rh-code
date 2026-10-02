#!/usr/bin/env python3
"""Window maxima of log|hat f_{kappa*}(v)|/v, after Proposition rz:prop:harmonic.

hat f_kappa(v) = prod_{k>=1} cos(v/(2(k+kappa))) is the characteristic function of the
random harmonic series H_kappa of Proposition rz:prop:harmonic, and kappa* = 6.61845...
is the solution of A(kappa) = -(1 + log 2 pi)/2 in part (i), where
    A(kappa) = A(0) - (psi(1+kappa) + gamma_E)/2,
    A(0) = (gamma_E - 2 log 2 + int_1^inf y^-2 log(1+e^{-2y}) dy + int_0^1 y^-2 log cosh y dy)/2.
Part (ii) says that limsup_{v->inf} v^{-1} log|hat f_kappa(v)| = -pi/4.

For each window [lo, hi] the script computes max_{lo <= v <= hi} L(v)/v, where
    L(v) = log|hat f_{kappa*}(v)| = sum_{k<=K} log|cos(v/(2(k+kappa*)))| + tail(v),
with K = ceil(2 vmax) direct terms and the tail k > K summed exactly through
Hurwitz zeta values: from log cos x = -sum_{m>=1} (1-4^{-m}) zeta(2m) (2x/pi)^{2m}/m,
valid for |x| < pi/2,
    tail(v) = -sum_{m>=1} (1-4^{-m}) zeta(2m)/m * zeta(2m, K+1+kappa*) * (v/pi)^{2m}.
Working with L avoids the underflow of hat f itself, which is about e^{-pi v/4}.

The maximum is found exactly, not on a grid.  hat f vanishes at v = (2j-1) pi (k+kappa*),
j, k >= 1, and between consecutive zeros L is concave (each log|cos| is), so
g(v) = v L'(v) - L(v) is decreasing there (g' = v L'') and L(v)/v has a single
critical point, its maximum, in each zero free interval.  The script lists every
zero in the window, finds the root of g in each interval (or takes the window end
when g does not change sign) in double precision with brentq, and recomputes L/v at
the five best candidates with mpmath at 30 digits, checking that the neighbours at
distance 1e-6 are lower.

Windows: [0.95V, 1.05V] for V = 100, 300, 600, 1000 and 10000 (the values quoted in
the book) and for V = 2000 and 5000 (the trend), and, for V = 100, 300, 600, 1000, a
family of other windows, to see whether any of them reproduces the values
-0.700, -0.758, -0.771, -0.779 printed before this script existed.

Writes harmonic_decay.json in the working directory.  Requires numpy, scipy and
mpmath 1.3; about half a minute on an Intel node of Hopper.
    python3 harmonic_decay.py
"""
import json
import math
import platform
import time

import numpy as np
import scipy
from scipy.optimize import brentq
import mpmath
from mpmath import mp, mpf

mp.dps = 30
T0 = time.time()

# ---------------------------------------------------------------- A(0) and kappa*
I1 = mp.quad(lambda y: mp.log(1 + mp.exp(-2 * y)) / y ** 2, [1, 10, mp.inf])
I2 = mp.quad(lambda y: mp.log(mp.cosh(y)) / y ** 2, [0, 1])
A0 = (mp.euler - 2 * mp.log(2) + I1 + I2) / 2
TARGET = -(1 + mp.log(2 * mp.pi)) / 2
kappa = mp.findroot(lambda k: A0 - (mp.digamma(1 + k) + mp.euler) / 2 - TARGET, mpf('6.6'))
assert abs(A0 - mpf('-0.1485757')) < mpf('1e-7'), A0
assert abs(kappa - mpf('6.61845')) < mpf('1e-5'), kappa
KAP = float(kappa)
print('A(0) =', mp.nstr(A0, 15), ' kappa* =', mp.nstr(kappa, 15), flush=True)

MTAIL = 40


class Evaluator:
    """L(v), L'(v) and g(v) = v L'(v) - L(v) in double precision for 0 < v <= vmax."""

    def __init__(self, vmax, K=None):
        self.K = int(math.ceil(2 * vmax)) if K is None else K
        self.c = 2.0 * (np.arange(1, self.K + 1) + KAP)
        a = mpf(self.K + 1) + kappa
        self.a_mp = a
        self.a = float(a)
        # tail(v) = -sum_m C_m (v/(pi a))^{2m}, C_m = (1-4^-m) zeta(2m)/m a^{2m} zeta(2m, a)
        self.C_mp = [(1 - mpf(4) ** (-m)) * mp.zeta(2 * m) / m * a ** (2 * m) * mp.zeta(2 * m, a)
                     for m in range(1, MTAIL + 1)]
        self.C = np.array([float(x) for x in self.C_mp])
        self.m = np.arange(1, MTAIL + 1)
        assert vmax / (math.pi * self.a) < 0.2

    def L(self, v):
        r = (v / (math.pi * self.a)) ** 2
        return float(np.sum(np.log(np.abs(np.cos(v / self.c)))) - np.sum(self.C * r ** self.m))

    def dL(self, v):
        r = (v / (math.pi * self.a)) ** 2
        return float(-np.sum(np.tan(v / self.c) / self.c) - np.sum(2 * self.m * self.C * r ** self.m) / v)

    def g(self, v):
        return v * self.dL(v) - self.L(v)

    # the same in mpmath
    def L_mp(self, v):
        v = mpf(v)
        s = mp.fsum(mp.log(abs(mp.cos(v / (2 * (k + kappa))))) for k in range(1, self.K + 1))
        r = (v / (mp.pi * self.a_mp)) ** 2
        return s - mp.fsum(C * r ** m for m, C in enumerate(self.C_mp, 1))


def zeros_in(lo, hi):
    """The zeros (2j-1) pi (k+kappa*) of hat f in the open interval (lo, hi)."""
    zs = []
    k = 1
    while math.pi * (k + KAP) < hi:
        base = math.pi * (k + KAP)
        j = max(1, math.ceil((lo / base + 1) / 2))
        while (2 * j - 1) * base < hi:
            z = (2 * j - 1) * base
            if z > lo:
                zs.append(z)
            j += 1
        k += 1
    return sorted(zs)


def window_max(ev, lo, hi, refine=5):
    zs = zeros_in(lo, hi)
    pts = [lo] + zs + [hi]
    cands = []
    for i in range(len(pts) - 1):
        l, r = pts[i], pts[i + 1]
        if r <= l:
            continue
        eps = min(1e-9 * l, (r - l) / 4)
        a_ = l + eps if i > 0 else l
        b_ = r - eps if i < len(pts) - 2 else r
        ga, gb = ev.g(a_), ev.g(b_)
        if ga <= 0:
            vs = a_
        elif gb >= 0:
            vs = b_
        else:
            vs = brentq(ev.g, a_, b_, xtol=1e-13, maxiter=300)
        cands.append((ev.L(vs) / vs, vs, l, r))
    cands.sort(reverse=True)
    best = []
    for val, vs, l, r in cands[:refine]:
        vm = mpf(vs)
        value_mp = ev.L_mp(vm) / vm
        # a local maximum: the neighbours at distance 1e-6, inside the interval, are lower
        nb = [x for x in (vm - mpf('1e-6'), vm + mpf('1e-6')) if l < x < r and lo <= x <= hi]
        is_max = all(ev.L_mp(x) / x < value_mp for x in nb)
        best.append({'v': vs, 'value_double': val, 'value_mp': float(value_mp),
                     'mp_local_max_check': is_max, 'at_window_end': not (lo < vs < hi),
                     'interval': [l, r]})
    best.sort(key=lambda d: -d['value_mp'])
    return {'window': [lo, hi], 'zeros_in_window': len(zs), 'intervals': len(cands),
            'max': best[0]['value_mp'], 'argmax': best[0]['v'],
            'max_double': cands[0][0], 'best': best}


# ---------------------------------------------------------------- checks of the evaluator
ev_chk = Evaluator(1100)
ev_big = Evaluator(1100, K=200000)
tail_dev = max(abs(ev_chk.L(v) - ev_big.L(v)) for v in (95.3, 301.7, 612.9, 1049.1))
mp_dev = max(abs(ev_chk.L(v) - float(ev_chk.L_mp(v))) for v in (95.3, 1049.1))
# a direct product in mpmath (no logarithms, no tail formula) at a moderate v
v_t = mpf('37.25')
prod_direct = mp.fprod(mp.cos(v_t / (2 * (k + kappa))) for k in range(1, 400001))
# the factors k > 400000 multiply by exp(-v^2/8 sum 1/(k+kappa)^2), up to about 1e-13 in the logarithm
prod_direct *= mp.exp(-v_t ** 2 / 8 * mp.zeta(2, 400001 + kappa))
prod_dev = abs(mp.log(abs(prod_direct)) - ev_chk.L_mp(v_t))
# zeros of hat f are where the text puts them: first zero (k=1, j=1) at pi(1+kappa*)
z1 = math.pi * (1 + KAP)
zero_check = ev_chk.L(z1 * (1 + 1e-12))
assert tail_dev < 1e-9 and mp_dev < 1e-9 and prod_dev < mpf('1e-12') and zero_check < -20
print(f'checks: tail {tail_dev:.1e}, double vs mpmath {mp_dev:.1e}, '
      f'direct product {float(prod_dev):.1e}', flush=True)

# ---------------------------------------------------------------- the windows
BOOK = {'100': -0.700, '300': -0.758, '600': -0.771, '1000': -0.779}
main = {}
for V in (100, 300, 600, 1000, 2000, 5000, 10000):
    t1 = time.time()
    ev = Evaluator(1.05 * V)
    res = window_max(ev, 0.95 * V, 1.05 * V)
    res['seconds'] = round(time.time() - t1, 1)
    main[str(V)] = res
    print(f"V={V:>6d} [0.95V,1.05V]: max L/v = {res['max']:.6f} at v = {res['argmax']:.4f} "
          f"({res['zeros_in_window']} zeros; double {res['max_double']:.6f}); "
          f"-pi/4 + 1/V * {(res['max'] + math.pi / 4) * V:.2f}", flush=True)

ALT = [('[0.9V,1.1V]', lambda V: (0.9 * V, 1.1 * V)),
       ('[0.975V,1.025V]', lambda V: (0.975 * V, 1.025 * V)),
       ('[0.8V,1.2V]', lambda V: (0.8 * V, 1.2 * V)),
       ('[0.9V,V]', lambda V: (0.9 * V, V)), ('[V,1.1V]', lambda V: (V, 1.1 * V)),
       ('[0.95V,V]', lambda V: (0.95 * V, V)), ('[V,1.05V]', lambda V: (V, 1.05 * V)),
       ('[V-5,V+5]', lambda V: (V - 5, V + 5)), ('[V-10,V+10]', lambda V: (V - 10, V + 10)),
       ('[V-20,V+20]', lambda V: (V - 20, V + 20)), ('[V-25,V+25]', lambda V: (V - 25, V + 25)),
       ('[V-50,V+50]', lambda V: (V - 50, V + 50)), ('[V,V+10]', lambda V: (V, V + 10)),
       ('[V,V+20]', lambda V: (V, V + 20)), ('[V,V+50]', lambda V: (V, V + 50)),
       ('[V-10,V]', lambda V: (V - 10, V)), ('[V-20,V]', lambda V: (V - 20, V)),
       ('[V-50,V]', lambda V: (V - 50, V))]
alt = {}
for name, win in ALT:
    alt[name] = {}
    for V in (100, 300, 600, 1000):
        lo, hi = win(V)
        ev = Evaluator(hi)
        r = window_max(ev, lo, hi, refine=2)
        alt[name][str(V)] = {'window': [lo, hi], 'max': r['max'], 'argmax': r['argmax']}
    vals = [alt[name][str(V)]['max'] for V in (100, 300, 600, 1000)]
    match = all(abs(round(x, 3) - BOOK[str(V)]) < 1e-9 for x, V in zip(vals, (100, 300, 600, 1000)))
    alt[name]['reproduces_printed_values'] = match
    print(f"{name:>16s}: " + ', '.join(f'{x:.4f}' for x in vals) + ('   <- printed' if match else ''),
          flush=True)

out = {
    'script': 'harmonic_decay.py',
    'what': 'max of log|hat f_{kappa*}(v)|/v over windows, Chapter rz3 after Proposition rz:prop:harmonic',
    'A0': float(A0), 'kappa_star': mp.nstr(kappa, 20), 'minus_pi_over_4': -math.pi / 4,
    'checks': {'tail_formula_vs_200000_direct_terms_max_dev': tail_dev,
               'double_vs_mpmath_max_dev': mp_dev,
               'log_of_direct_mpmath_product_vs_L_at_37.25': float(prod_dev),
               'L_just_above_first_zero': zero_check},
    'main_windows_0.95V_1.05V': main,
    'book_values_before': BOOK,
    'other_windows': alt,
    'software': {'python': platform.python_version(), 'numpy': np.__version__,
                 'scipy': scipy.__version__, 'mpmath': mpmath.__version__,
                 'host': platform.node()},
    'seconds_total': round(time.time() - T0, 1),
}
with open('harmonic_decay.json', 'w') as fh:
    json.dump(out, fh, indent=1)
print('wrote harmonic_decay.json in', round(time.time() - T0, 1), 's')
