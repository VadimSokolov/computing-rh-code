#!/usr/bin/env python3
"""Truncation errors of the sawtooth expansion of S, after Proposition rz:prop:sawtrunc.

Chapter rz3 compares two truncations of the formal expansion of
S(tau) = pi^{-1} arg zeta(1/2 + i tau) on 10 <= tau <= 50:

  the sawtooth series truncated by the index N,
      S_saw(tau) = sum_{2 <= N <= N0} w_N B_N Sigma_{log N}(tau),
  and the classical series,
      S_cl(tau)  = -(1/pi) sum_{2 <= n <= N0} w_n a_n sin(tau log n),

with the Cesaro weight w_N = 1 - log N / log N0.  Here a_n = Lambda(n)/(sqrt(n) log n),
B_N = sum_{n^m = N} d_m a_n (Theorem rz:thm:saw), d_m is defined by
sum_m d_m m^{-v} = 1/eta(1+v), that is d_m = sum_{2^j m' = m} mu(m')/m'
(Proposition rz:prop:eta), and
    Sigma_lambda(tau) = #{k >= 0 : (2k+1) pi / lambda <= tau} - lambda tau / (2 pi)
is the sawtooth of Lemma rz:lem:saw, which equals floor(y + 1/2) - y with
y = lambda tau / (2 pi).

S is computed as N(tau) - 1 - theta(tau)/pi, with theta the Riemann Siegel theta
function evaluated by mpmath at 25 digits and N(tau) the exact count of the ten
zeros with ordinate below 50 (mpmath.zetazero, checked against Odlyzko's table).
It is checked against pi^{-1} times the principal argument of zeta(1/2 + i tau),
computed by mpmath at 300 random points; the two agree wherever |S| < 1, and the
script records max |S| on the grid.

Root mean squares are taken over the uniform grids
tau_i = 10 + i h, i = 0, ..., 40/h, for h = 1e-3, 5e-4 and 2.5e-4, and over the grid of
step 1e-3 shifted by 5e-4 (tau = 10.0005, ..., 49.9995), so that the agreement of the
printed digits across grids shows that they are converged.  The script also computes
the errors without the weight, and (1/2) sum_{N <= N0} |B_N|, and it includes
N0 = 100 to test whether the values printed before this script existed (0.230,
0.228, 0.223 for N0 = 10^3, 10^4, 10^5) were shifted by one power of ten.

Writes sawtooth_truncation.json in the working directory.  Requires numpy and
mpmath 1.3; about a minute on an Intel node of Hopper.
    python3 sawtooth_truncation.py
"""
import json
import math
import platform
import random
import time

import numpy as np
import mpmath
from mpmath import mp, mpf, mpc, zetazero, siegeltheta, zeta, arg

mp.dps = 25
T0 = time.time()

# ---------------------------------------------------------------- the zeros
ODLYZKO = ['14.134725142', '21.022039639', '25.010857580', '30.424876126',
           '32.935061588', '37.586178159', '40.918719012', '43.327073281',
           '48.005150881', '49.773832478', '52.970321478']
zeros = [zetazero(n).imag for n in range(1, len(ODLYZKO) + 1)]
zero_dev = max(abs(z - mpf(o)) for z, o in zip(zeros, ODLYZKO))
assert zero_dev < mpf('1e-8'), zero_dev
gam = np.array([float(z) for z in zeros if z < 50])
assert len(gam) == 10 and float(zeros[10]) > 50

# ---------------------------------------------------------------- the grid and S
A, B = 10, 50
DEN = 4000                       # finest step 1/4000 = 2.5e-4
n_f = (B - A) * DEN + 1
tau = A + np.arange(n_f) / DEN
theta = np.array([float(siegeltheta(mpf(A) + mpf(i) / DEN)) for i in range(n_f)])
S = np.searchsorted(gam, tau, side='right') - 1 - theta / math.pi
t_theta = time.time() - T0

GRIDS = [('h=1e-3', slice(0, None, 4)),
         ('h=1e-3 shifted by 5e-4', slice(2, None, 4)),
         ('h=5e-4', slice(0, None, 2)),
         ('h=2.5e-4', slice(0, None, 1))]
grid_info = {name: {'points': int(tau[sl].size), 'first': float(tau[sl][0]),
                    'last': float(tau[sl][-1])} for name, sl in GRIDS}

# check S against the principal argument of zeta
rng = random.Random(20260927)
arg_dev, n_arg = mpf(0), 0
while n_arg < 300:
    t = mpf(rng.uniform(A, B))
    if min(abs(t - z) for z in zeros) < mpf('1e-6'):
        continue
    s_count = sum(1 for z in zeros if z <= t) - 1 - siegeltheta(t) / mp.pi
    s_arg = arg(zeta(mpc(mpf(1) / 2, t))) / mp.pi
    arg_dev = max(arg_dev, abs(s_count - s_arg))
    n_arg += 1
assert arg_dev < mpf('1e-15'), arg_dev
# a few values of S for the record
S_samples = {str(t): float(sum(1 for z in zeros if z <= t) - 1 - siegeltheta(t) / mp.pi)
             for t in (10, 20, 30, 40, 50)}


# ---------------------------------------------------------------- the weights B_N
def mobius(n):
    r, k, p = 1, n, 2
    while p * p <= k:
        if k % p == 0:
            k //= p
            if k % p == 0:
                return 0
            r = -r
        p += 1
    return -r if k > 1 else r


def d(m):
    """d_m = sum_{2^j m' = m} mu(m')/m' (Proposition rz:prop:eta)."""
    s, j = 0.0, 0
    while m % (2 ** j) == 0:
        mm = m // 2 ** j
        s += mobius(mm) / mm
        j += 1
    return s


def B_pp(p, e):
    """B_{p^e} = sum_{m | e} d_m a_{p^{e/m}}, with a_{p^r} = p^{-r/2}/r."""
    return sum(d(m) * p ** (-(e // m) / 2) / (e // m) for m in range(1, e + 1) if e % m == 0)


def necklace_b(n, q):
    """b_n(q) = sum_{2^i | n} M_{n/2^i}(q), M_k(q) = (1/k) sum_{c | k} mu(k/c) q^c (Theorem rz:thm:necklace)."""
    def M(k):
        return sum(mobius(k // c) * q ** c for c in range(1, k + 1) if k % c == 0) / k
    s, i = 0.0, 0
    while n % (2 ** i) == 0:
        s += M(n // 2 ** i)
        i += 1
    return s


NMAX = 10 ** 5
sieve = np.ones(NMAX + 1, bool)
sieve[:2] = False
for p in range(2, int(NMAX ** 0.5) + 1):
    if sieve[p]:
        sieve[p * p::p] = False
pp = []                          # (N, p, e) for the prime powers N = p^e <= NMAX
for p in np.nonzero(sieve)[0]:
    p, e, N = int(p), 1, int(p)
    while N <= NMAX:
        pp.append((N, p, e))
        e += 1
        N *= p
pp.sort()
N_all = np.array([x[0] for x in pp], dtype=np.int64)
lam_all = np.log(N_all.astype(float))
B_all = np.array([B_pp(p, e) for _, p, e in pp])
a_all = np.array([p ** (-e / 2) / e for _, p, e in pp])       # a_{p^e} = Lambda/(sqrt(N) log N)
neck_dev = max(abs(B_pp(p, e) - necklace_b(e, p ** -0.5)) for _, p, e in pp)
first_B = {'B_2': B_pp(2, 1), 'B_4': B_pp(2, 2), 'B_8': B_pp(2, 3), 'B_9': B_pp(3, 2)}
first_B_text = {'B_2': 2 ** -0.5, 'B_4': 0.25 + 1 / (2 * math.sqrt(2)),
                'B_8': -1 / (6 * math.sqrt(2)), 'B_9': 1 / 6 + 1 / (2 * math.sqrt(3))}
first_B_dev = max(abs(first_B[k] - first_B_text[k]) for k in first_B)
assert neck_dev < 1e-12 and first_B_dev < 1e-15

# the sawtooth formula floor(y + 1/2) - y against the count in Lemma rz:lem:saw
rng2 = random.Random(7)
saw_dev = 0.0
for _ in range(2000):
    lam_, t_ = rng2.uniform(0.5, 12.0), rng2.uniform(A, B)
    cnt = 0
    while (2 * cnt + 1) * math.pi / lam_ <= t_:
        cnt += 1
    y = lam_ * t_ / (2 * math.pi)
    saw_dev = max(saw_dev, abs((cnt - y) - (math.floor(y + 0.5) - y)))
assert saw_dev < 1e-12


# ---------------------------------------------------------------- the truncations
def rms(x):
    return float(np.sqrt(np.mean(x * x)))


BOOK = {'saw_cesaro': {'1000': 0.230, '10000': 0.228, '100000': 0.223},
        'classical_cesaro': {'1000': 0.147, '10000': 0.128, '100000': 0.116},
        'rms_S': 0.324,
        'half_sum_abs_B': {'1000': 9.0, '10000': 18.5, '100000': 40.6}}

results = {}
CH = 64
for N0 in (10 ** 2, 10 ** 3, 10 ** 4, 10 ** 5):
    t1 = time.time()
    sel = N_all <= N0
    lam, Bv, av = lam_all[sel], B_all[sel], a_all[sel]
    w = 1.0 - lam / math.log(N0)
    saw_w = np.zeros(n_f)
    saw_1 = np.zeros(n_f)
    cl_w = np.zeros(n_f)
    cl_1 = np.zeros(n_f)
    for c0 in range(0, lam.size, CH):
        c = slice(c0, c0 + CH)
        X = np.outer(lam[c], tau)                  # lambda tau
        Y = X / (2 * math.pi)
        saw = np.floor(Y + 0.5) - Y                # Sigma_lambda(tau)
        sn = np.sin(X)
        saw_w += (w[c] * Bv[c]) @ saw
        saw_1 += Bv[c] @ saw
        cl_w -= (w[c] * av[c]) @ sn / math.pi
        cl_1 -= av[c] @ sn / math.pi
    per_grid = {}
    for name, sl in GRIDS:
        per_grid[name] = {
            'rms_S': rms(S[sl]),
            'rms_err_saw_cesaro': rms(saw_w[sl] - S[sl]),
            'rms_err_classical_cesaro': rms(cl_w[sl] - S[sl]),
            'rms_err_saw_unweighted': rms(saw_1[sl] - S[sl]),
            'rms_err_classical_unweighted': rms(cl_1[sl] - S[sl]),
        }
    results[str(N0)] = {'prime_powers': int(lam.size),
                        'half_sum_abs_B': float(0.5 * np.abs(Bv).sum()),
                        'grids': per_grid,
                        'seconds': round(time.time() - t1, 1)}
    g = per_grid['h=1e-3']
    print(f"N0={N0:>6d}  saw(Cesaro) {g['rms_err_saw_cesaro']:.6f}  classical(Cesaro) "
          f"{g['rms_err_classical_cesaro']:.6f}  saw(w=1) {g['rms_err_saw_unweighted']:.6f}  "
          f"classical(w=1) {g['rms_err_classical_unweighted']:.6f}  "
          f"(1/2)sum|B| {0.5 * np.abs(Bv).sum():.4f}", flush=True)

print('rms S on the grids:', {name: round(rms(S[sl]), 6) for name, sl in GRIDS})
print('grid spread of each quantity (max - min over the four grids):')
spread = {}
for N0, r in results.items():
    spread[N0] = {}
    for q in ('rms_S', 'rms_err_saw_cesaro', 'rms_err_classical_cesaro',
              'rms_err_saw_unweighted', 'rms_err_classical_unweighted'):
        vals = [r['grids'][name][q] for name, _ in GRIDS]
        spread[N0][q] = max(vals) - min(vals)
    print(N0, {k: f"{v:.2e}" for k, v in spread[N0].items()})

out = {
    'script': 'sawtooth_truncation.py',
    'what': 'RMS of S and RMS errors of the truncated sawtooth and classical series on '
            'tau in [10,50], Chapter rz3 after Proposition rz:prop:sawtrunc',
    'definitions': {
        'S': 'N(tau) - 1 - theta(tau)/pi, exact count of the ten zeros below 50',
        'saw': 'sum_{2<=N<=N0} w_N B_N Sigma_{log N}(tau), Sigma = floor(y+1/2) - y, y = tau log N/(2 pi)',
        'classical': '-(1/pi) sum_{2<=n<=N0} w_n Lambda(n) n^{-1/2} (log n)^{-1} sin(tau log n)',
        'cesaro_weight': 'w_N = 1 - log N/log N0 (unweighted: w = 1)',
        'rms': 'square root of the mean of squares over the grid points',
    },
    'grids': grid_info,
    'checks': {
        'zeros_vs_odlyzko_max_dev': float(zero_dev),
        'zeros_below_50': [float(z) for z in zeros[:10]],
        'S_vs_principal_arg_zeta_max_dev_300_random_points': float(arg_dev),
        'max_abs_S_on_finest_grid': float(np.abs(S).max()),
        'mean_S_on_finest_grid': float(S.mean()),
        'S_at': S_samples,
        'B_N_vs_necklace_formula_max_dev': float(neck_dev),
        'first_B_vs_text_max_dev': float(first_B_dev),
        'sawtooth_formula_vs_count_max_dev': float(saw_dev),
    },
    'results': results,
    'grid_spread': spread,
    'book_values_before': BOOK,
    'software': {'python': platform.python_version(), 'numpy': np.__version__,
                 'mpmath': mpmath.__version__, 'host': platform.node()},
    'seconds_theta': round(t_theta, 1),
    'seconds_total': round(time.time() - T0, 1),
}
with open('sawtooth_truncation.json', 'w') as fh:
    json.dump(out, fh, indent=1)
print('wrote sawtooth_truncation.json in', round(time.time() - T0, 1), 's')
