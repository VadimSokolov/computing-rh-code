"""Floating point multistart search for the torus extremes of Theorem rz:thm:torus.

For each n and band [a, b] it estimates
  m = min |P_n(x,w)|,  M = max -Re(P_n'/P_n)(x,w),  l = max |P_n'/P_n|(x,w)
over w in the torus and x in [a, b], by 300000 random samples followed by L-BFGS-B from the best
400 samples and from 1500 random starts. Floating point, not certified; used to locate the extremes.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import json, sys, time
import numpy as np
from multiprocessing import Pool
from scipy.optimize import minimize
from torus_common import Torus

NS = [7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 32, 36, 40]


def objective(T, kind):
    def f(z):
        x, w = z[0], z[1:]
        P, P1, dP, dP1 = T.PP_grad(x, w)
        if kind == "m":
            val = abs(P) ** 2
            g = 2 * np.real(np.conj(P) * dP)
        else:
            R = P1 / P
            dR = (dP1 * P - P1 * dP) / P ** 2
            if kind == "M":  # minimise Re(P'/P) = maximise -Re(P'/P)
                val = R.real
                g = dR.real
            else:  # "l": maximise |P'/P|^2
                val = -abs(R) ** 2
                g = -2 * np.real(np.conj(R) * dR)
        return val, g
    return f


def vec_obj(T, kind, X, W):
    P, P1 = T.PP(X, W)
    if kind == "m":
        return np.abs(P) ** 2
    R = P1 / P
    return R.real if kind == "M" else -np.abs(R) ** 2


def search(args):
    n, a, b, kind, seed = args
    T = Torus(n)
    rng = np.random.default_rng(seed)
    f = objective(T, kind)
    best = None
    # stage 1: random samples
    N = 300000
    X = rng.uniform(a, b, N)
    W = rng.uniform(0, 2 * np.pi, (N, T.K))
    vals = np.concatenate([vec_obj(T, kind, X[i:i + 50000], W[i:i + 50000]) for i in range(0, N, 50000)])
    idx = np.argsort(vals)[:400]
    starts = [np.concatenate([[X[i]], W[i]]) for i in idx]
    starts += [np.concatenate([[rng.uniform(a, b)], rng.uniform(0, 2 * np.pi, T.K)]) for _ in range(1500)]
    starts += [np.concatenate([[a], np.zeros(T.K)]), np.concatenate([[b], np.zeros(T.K)])]
    bounds = [(a, b)] + [(None, None)] * T.K
    results = []
    for z0 in starts:
        res = minimize(f, z0, jac=True, method="L-BFGS-B", bounds=bounds, options={"maxiter": 500, "ftol": 1e-15, "gtol": 1e-11})
        results.append((res.fun, res.x))
    results.sort(key=lambda t: t[0])
    val, z = results[0]
    # count distinct local optima near the best
    z = z.copy()
    z[1:] = np.mod(z[1:], 2 * np.pi)
    if kind == "m":
        out = np.sqrt(val)
    elif kind == "M":
        out = -val
    else:
        out = np.sqrt(-val)
    top = [(-v if kind == "M" else (np.sqrt(v) if kind == "m" else np.sqrt(-v))) for v, _ in results[:10]]
    return {"n": n, "band": [a, b], "kind": kind, "value": float(out), "x": float(z[0]),
            "w": [float(t) for t in z[1:]], "primes": T.ps, "top10": [float(t) for t in top]}


if __name__ == "__main__":
    tasks = []
    for n in NS:
        for kind in ("m", "M", "l"):
            if n >= 16:
                tasks.append((n, 1.2, 1.3, kind, 1000 * n + ord(kind)))
            tasks.append((n, 1.0, 1.5, kind, 2000 * n + ord(kind)))
    t0 = time.time()
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 32) as pool:
        res = pool.map(search, tasks, chunksize=1)
    json.dump(res, open("explore.json", "w"), indent=1)
    for r in res:
        print(f"n={r['n']:2d} band={r['band']} {r['kind']}: {r['value']:.6f} at x={r['x']:.5f} w={np.round(r['w'], 4).tolist()}  top={np.round(r['top10'][:4], 6).tolist()}")
    print("time", time.time() - t0)
