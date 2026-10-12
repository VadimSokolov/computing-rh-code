# wigbiane_ineq.py: numerical tests of the power sum inequalities behind the collision bounds of the section on Laguerre's
# quantity and the collisions of zeros along the de Bruijn Newman flow (ch/ch21.tex, after Section sec:ch21:hankel; draft
# misc/reviews/nick-wigner/merge/wigbiane.tex). Added for the book in October 2026, from the first author's note on the
# semicircle law and the Thorin framework (misc/reviews/nick-wigner/, version 7, Section 8), whose proofs the draft checks line
# by line; the computations here are independent of the note's scripts.
# For a real sequence x_1, ..., x_n write S_k = sum_j x_j^k. Part A minimises E/S_2^3, E = 3 S_2^3 - 6 S_2 S_4 + 4 S_3^2 (the
# sextic inequality says E >= 0, with equality only for a pair {a, -a}), by BFGS with the exact gradient from random starts of
# five kinds (Gaussian, Cauchy, sparse, a perturbed pair, a perturbed pair of equal signs), n = 2, ..., 20, once without and
# once with the constraint S_1 = 0 (SLSQP). Part B maximises 2 S_4/S_2^2 subject to S_3 = 0 (the quartic inequality says at
# most 1, with equality only for a pair). Part C maximises (S_4 - S_2^2/2)/(S_2^{1/2}|S_3|) (the shadow lemma (i) says at
# most 1). Part D checks the shadow lemma (ii), |Re (1 - i e)^{-k} - 1| <= k(k+1) e^2/2, for k = 1, ..., 8 on a grid of
# 4001 values of e in [10^-6, 10^4] at 40 digits. Part E draws 40000 random configurations of real zeros and conjugate pairs,
# takes a gap between two adjacent real zeros, finds a critical point u of Lambda(x) = sum (x - z)^{-2} over the gap
# (bisection safeguarded Newton from the least value on a grid of 1999 points), drops the configurations with Lambda(u) <= 0
# (possible when a pair lies close to the gap) and those with |S_3| >= 10^-9 S_2^{3/2} at the computed point, and checks
# the two inequalities used in the local and upper bounds at u: 6 S_4 - 3 S_2^2 <= eta S_2^2 with
# eta = (36 m B_4 + 108 B_4^2 + 72 (m + 6 B_4)^{1/2} B_5 + 120 B_6)/m^2, B_p = sum b^2/|u - a|^p over the pairs a +- ib, and
# 6 S_4 >= 3 P^2 - 6 W, P the share of the two zeros bounding the gap and W = sum |u - z|^{-4} over the nonreal zeros.
# The best points of Parts A to C are re-evaluated at 50 digits. Run from any folder: python3 wigbiane_ineq.py; writes
# wigbiane_ineq.json there. Uses SLURM_CPUS_PER_TASK processes (default 1).
import json
import multiprocessing as mpr
import os
import time

import mpmath as mp
import numpy as np
from scipy.optimize import minimize

t0 = time.time()
out = {}


def rec(k, v):
    out[k] = v
    print(k, v, round(time.time() - t0, 1), flush=True)
    json.dump(out, open('wigbiane_ineq.json', 'w'), indent=1)


def sums(x):
    x2 = x * x
    return x2.sum(), (x2 * x).sum(), (x2 * x2).sum()


def sextic(x):
    S2, S3, S4 = sums(x)
    E = 3 * S2 ** 3 - 6 * S2 * S4 + 4 * S3 ** 2
    dE = 18 * S2 ** 2 * x - 6 * (2 * x * S4 + 4 * S2 * x ** 3) + 24 * S3 * x ** 2
    return E / S2 ** 3, dE / S2 ** 3 - 6 * E * x / S2 ** 4


def mp_sextic(x):
    with mp.workdps(50):
        v = [mp.mpf(float(a)) for a in x]
        S2 = sum(a ** 2 for a in v)
        S3 = sum(a ** 3 for a in v)
        S4 = sum(a ** 4 for a in v)
        return float((3 * S2 ** 3 - 6 * S2 * S4 + 4 * S3 ** 2) / S2 ** 3)


def start(n, kind, rng):
    if kind == 0:
        return rng.standard_normal(n)
    if kind == 1:
        return rng.standard_cauchy(n)
    if kind == 2:
        x = np.zeros(n)
        k = rng.integers(1, n + 1)
        x[rng.choice(n, k, replace=False)] = rng.standard_normal(k)
        if not x.any():
            x[0] = 1.0
        return x
    x = rng.standard_normal(n) * 10.0 ** rng.uniform(-4, 0)
    x[0] += 1.0
    x[1] += -1.0 if kind == 3 else 1.0
    return x


def part_a(args):
    n, seed, nstart, constrained = args
    rng = np.random.default_rng(seed)
    best = (np.inf, None)
    lowest_eval = np.inf
    for i in range(nstart):
        x0 = start(n, i % 5, rng)
        if constrained:
            x0 = x0 - x0.mean()
            if not np.any(np.abs(x0) > 1e-12):
                continue
            x0 = x0 / np.sqrt((x0 ** 2).sum())
            r = minimize(lambda x: sextic(x)[0], x0, jac=lambda x: sextic(x)[1], method='SLSQP',
                         constraints=[{'type': 'eq', 'fun': lambda x: x.sum(), 'jac': lambda x: np.ones_like(x)},
                                      {'type': 'eq', 'fun': lambda x: (x * x).sum() - 1, 'jac': lambda x: 2 * x}],
                         options={'maxiter': 500, 'ftol': 1e-15})
        else:
            x0 = x0 / np.sqrt((x0 ** 2).sum())
            r = minimize(lambda x: sextic(x)[0], x0, jac=lambda x: sextic(x)[1], method='BFGS', options={'gtol': 1e-13})
        x = r.x / np.sqrt((r.x ** 2).sum())
        f = sextic(x)[0]
        lowest_eval = min(lowest_eval, f)
        if f < best[0]:
            best = (f, x)
    f, x = best
    order = np.argsort(-np.abs(x))
    xs = x[order]
    return {'n': n, 'S1 = 0 imposed': constrained, 'starts': nstart, 'least E/S2^3 found': f,
            'least at 50 digits': mp_sextic(x), 'two largest terms': [float(xs[0]), float(xs[1])],
            'largest other |term|': float(np.abs(xs[2:]).max()) if n > 2 else 0.0}


def quartic_obj(x):
    S2, S3, S4 = sums(x)
    return -2 * S4 / S2 ** 2


def part_b(args):
    n, seed, nstart = args
    rng = np.random.default_rng(seed)
    best = (-np.inf, None)
    for i in range(nstart):
        x0 = start(n, i % 5, rng)
        x0 = x0 / np.sqrt((x0 ** 2).sum())
        r = minimize(lambda x: -2 * (x ** 4).sum(), x0, jac=lambda x: -8 * x ** 3, method='SLSQP',
                     constraints=[{'type': 'eq', 'fun': lambda x: (x ** 3).sum(), 'jac': lambda x: 3 * x ** 2},
                                  {'type': 'eq', 'fun': lambda x: (x * x).sum() - 1, 'jac': lambda x: 2 * x}],
                     options={'maxiter': 500, 'ftol': 1e-15})
        x = r.x
        S2, S3, S4 = sums(x)
        if abs(S3) > 1e-10 * S2 ** 1.5:
            continue
        v = 2 * S4 / S2 ** 2
        if v > best[0]:
            best = (v, x / np.sqrt(S2))
    v, x = best
    with mp.workdps(50):
        w = [mp.mpf(float(a)) for a in x]
        S2 = sum(a ** 2 for a in w)
        S3 = sum(a ** 3 for a in w)
        S4 = sum(a ** 4 for a in w)
        v50 = float(2 * S4 / S2 ** 2)
        s3 = float(S3)
    order = np.argsort(-np.abs(x))
    return {'n': n, 'starts': nstart, 'largest 2 S4/S2^2 found with S3 = 0': v, 'at 50 digits': v50, 'S3 there': s3,
            'two largest terms': [float(x[order[0]]), float(x[order[1]])]}


def shadow_ratio(x):
    S2, S3, S4 = sums(x)
    return (S4 - S2 ** 2 / 2) / (np.sqrt(S2) * abs(S3))


def part_c(args):
    n, seed, nstart = args
    rng = np.random.default_rng(seed)
    best = (-np.inf, None)
    for i in range(nstart):
        x0 = start(n, i % 5, rng)
        x0 = x0 / np.sqrt((x0 ** 2).sum())
        r = minimize(lambda x: -shadow_ratio(x), x0, method='Nelder-Mead',
                     options={'maxiter': 4000 * n, 'xatol': 1e-12, 'fatol': 1e-14})
        v = shadow_ratio(r.x)
        if np.isfinite(v) and v > best[0]:
            best = (v, r.x / np.sqrt((r.x ** 2).sum()))
    v, x = best
    order = np.argsort(-np.abs(x))
    return {'n': n, 'starts': nstart, 'largest (S4 - S2^2/2)/(S2^(1/2)|S3|) found': v,
            'three largest terms': [float(x[order[j]]) for j in range(min(3, n))]}


def part_e(seed):
    rng = np.random.default_rng(seed)
    nr = int(rng.integers(3, 12))
    nc = int(rng.integers(1, 8))
    real = np.sort(rng.uniform(-1, 1, nr) * rng.uniform(1, 4))
    a = rng.uniform(-3, 3, nc)
    b = np.abs(rng.standard_normal(nc)) * 10.0 ** rng.uniform(-3, 0)
    zs = np.concatenate([real, a + 1j * b, a - 1j * b])
    k = int(rng.integers(0, nr - 1))
    lo, hi = real[k], real[k + 1]
    if hi - lo < 1e-6:
        return None
    xs = lo + (hi - lo) * (0.5 - 0.5 * np.cos(np.pi * np.linspace(0, 1, 2001)[1:-1]))
    lam = np.real(((xs[:, None] - zs[None, :]) ** -2.0).sum(axis=1))
    i = int(np.argmin(lam))
    u = xs[i]
    a_, b_ = xs[max(i - 1, 0)], xs[min(i + 1, len(xs) - 1)]
    for _ in range(200):
        d = u - zs
        l1 = np.real((-2 * d ** -3.0).sum())
        l2 = np.real((6 * d ** -4.0).sum())
        if l1 > 0:
            b_ = u
        else:
            a_ = u
        un = u - l1 / l2 if l2 > 0 else (a_ + b_) / 2
        if not (a_ < un < b_):
            un = (a_ + b_) / 2
        if abs(un - u) < 1e-15 * (1 + abs(u)):
            u = un
            break
        u = un
    d = u - zs
    S2 = np.real((d ** -2.0).sum())
    S3 = np.real((d ** -3.0).sum())
    S4 = np.real((d ** -4.0).sum())
    if S2 <= 0:
        return None
    da = u - a
    B = dict((p, float((b ** 2 / np.abs(da) ** p).sum())) for p in (4, 5, 6))
    m = S2
    eta = (36 * m * B[4] + 108 * B[4] ** 2 + 72 * np.sqrt(m + 6 * B[4]) * B[5] + 120 * B[6]) / m ** 2
    P = (u - lo) ** -2 + (u - hi) ** -2
    W = float((2 * np.abs(a + 1j * b - u) ** -4.0).sum())
    lhs = 6 * S4 - 3 * S2 ** 2
    return {'S3 relative': abs(S3) / S2 ** 1.5, 'local: lhs/(eta m^2)': lhs / (eta * m * m),
            'upper': (3 * P * P - 6 * W - 6 * S4) / (3 * P * P + 6 * W + 6 * abs(S4)), 'eta': eta}


if __name__ == '__main__':
    nproc = int(os.environ.get('SLURM_CPUS_PER_TASK', '1'))
    pool = mpr.get_context('fork').Pool(nproc)
    ns = list(range(2, 21))
    ra = pool.map(part_a, [(n, 1000 + n, 2000, False) for n in ns])
    rec('A: sextic, least E/S2^3 over the sphere, no constraint', ra)
    rec('A: least value over all n', min(r['least at 50 digits'] for r in ra))
    rc = pool.map(part_a, [(n, 2000 + n, 1000, True) for n in ns])
    rec('A: sextic with S1 = 0 imposed', rc)
    rec('A: least value over all n with S1 = 0', min(r['least at 50 digits'] for r in rc))
    rb = pool.map(part_b, [(n, 3000 + n, 1000) for n in ns])
    rec('B: quartic, largest 2 S4/S2^2 subject to S3 = 0', rb)
    rec('B: largest value over all n', max(r['at 50 digits'] for r in rb))
    rcc = pool.map(part_c, [(n, 4000 + n, 200) for n in ns])
    rec('C: shadow lemma (i), largest ratio', rcc)
    rec('C: largest ratio over all n', max(r['largest (S4 - S2^2/2)/(S2^(1/2)|S3|) found'] for r in rcc))
    mp.mp.dps = 40
    dres = {}
    for k in range(1, 9):
        worst = mp.mpf(0)
        at = None
        for j in range(4001):
            e = mp.mpf(10) ** (mp.mpf(-6) + mp.mpf(10) * j / 4000)
            r = abs(mp.re((1 - 1j * e) ** (-k)) - 1) / (k * (k + 1) * e * e / 2)
            if r > worst:
                worst, at = r, e
        dres[str(k)] = {'largest ratio': mp.nstr(worst, 12), 'at e': mp.nstr(at, 6)}
    rec('D: shadow lemma (ii), largest |Re(1 - ie)^-k - 1|/(k(k+1)e^2/2) on the grid', dres)
    re_ = [r for r in pool.map(part_e, range(5000, 45000)) if r is not None]
    rec('E: configurations drawn', len(re_))
    rec('E: largest |S3|/S2^(3/2) at the computed minimisers', max(r['S3 relative'] for r in re_))
    re_ = [r for r in re_ if r['S3 relative'] < 1e-9]
    rec('E: configurations with |S3| < 1e-9 S2^(3/2) at the minimiser (kept)', len(re_))
    rec('E: largest (6 S4 - 3 S2^2)/(eta m^2) (local bound needs <= 1)', max(r['local: lhs/(eta m^2)'] for r in re_))
    rec('E: largest (3 P^2 - 6 W - 6 S4)/(3 P^2 + 6 W + 6 |S4|) (upper bound needs <= 0)', max(r['upper'] for r in re_))
    rec('E: largest eta', max(r['eta'] for r in re_))
    out['seconds'] = time.time() - t0
    json.dump(out, open('wigbiane_ineq.json', 'w'), indent=1)
    print('done', round(time.time() - t0, 1))
