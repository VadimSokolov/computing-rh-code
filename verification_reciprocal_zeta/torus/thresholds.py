"""Thresholds of Theorem rz:thm:torus under several readings, in floating point.

Reads explore.json (the multistart estimates of m, M, l on the bands [1.2,1.3] and [1,1.5]) and computes:
  y1 = root of the first condition  log(y/pi) - 2/y - M - (l r + r' + r/y)/(y (m - r/y)) = 0,
       with r = r(1.2), r' = r'(1.2) for the band [1.2,1.3] and r(1), r'(1) for the single band [1,1.5];
  y2 = root of the second condition log(y/pi) - 1/y - B_n(y) = 0 (band [1.2,1.3] only), with B_n(y) the maximum
       over 0.3 <= sigma <= 1/2 and the torus, found by multistart L-BFGS-B;
under the printed (rounded) m, M of Table rz:tab:torus and under the unrounded estimates.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import json, sys, math
import numpy as np
from multiprocessing import Pool
from scipy.optimize import minimize, brentq
from torus_common import Torus

PRINTED = {16: (1.06, 0.86, 30), 17: (1.00, 0.89, 34), 18: (0.95, 0.93, 37), 19: (0.90, 0.96, 42),
           20: (0.85, 0.98, 46), 21: (0.80, 1.01, 50), 22: (0.76, 1.04, 56), 23: (0.72, 1.06, 62),
           24: (0.68, 1.09, 69), 32: (0.36, 1.26, 214), 36: (0.20, 1.33, 566), 40: (0.056, 1.40, 4898)}
SINGLE = {7: 14, 9: 23, 10: 31, 11: 42, 12: 57, 13: 85, 14: 131, 15: 229}


def cond1(y, m, M, l, r, rp):
    if m - r / y <= 0:
        return -1e9
    return math.log(y / math.pi) - 2 / y - M - (l * r + rp + r / y) / (y * (m - r / y))


def root1(m, M, l, r, rp):
    lo = r / m * (1 + 1e-12)
    hi = lo * 2 + 10
    while cond1(hi, m, M, l, r, rp) <= 0:
        hi *= 2
    lo2 = lo
    # cond1 is increasing on (r/m, oo); find the sign change
    return brentq(lambda y: cond1(y, m, M, l, r, rp), lo2 * (1 + 1e-9) + 1e-9, hi, xtol=1e-12)


def make_q(T, y):
    """q(sigma, w) and gradient, to be maximised (we minimise -q)."""
    absA = np.abs(T.alpha)

    def rr(x):
        e = np.exp(-2 * x * T.logk)
        return float(np.sum(absA * e)), float(np.sum(2 * T.logk * absA * e))

    def f(z):
        sg, w = z[0], z[1:]
        x1, x2 = 1.5 - sg, 1 + sg
        P1, D1, dP1, _ = T.PP_grad(x1, w)
        P2, D2, dP2, _ = T.PP_grad(x2, w)
        r1, rp1 = rr(x1)
        r2, rp2 = rr(x2)
        A1, A2 = abs(P1), abs(P2)
        a1 = A1 + r1 / y
        a2 = A2 - r2 / y
        if a2 <= 0:
            return 1e6, np.zeros_like(z)
        den = 2 * sg - 0.5
        num = math.log(a1) - math.log(a2)
        q = num / den
        g = np.zeros_like(z)
        # w derivatives
        dA1w = np.real(np.conj(P1) * dP1[1:]) / A1
        dA2w = np.real(np.conj(P2) * dP2[1:]) / A2
        g[1:] = (dA1w / a1 - dA2w / a2) / den
        # sigma derivative: x1 = 1.5 - sigma, x2 = 1 + sigma; d r/dx = -r'
        da1 = -np.real(np.conj(P1) * D1) / A1 + rp1 / y
        da2 = np.real(np.conj(P2) * D2) / A2 + rp2 / y
        g[0] = (da1 / a1 - da2 / a2) / den - 2 * num / den ** 2
        return -q, -g
    return f


def Bn(T, y, rng, seeds=None, nrand=300):
    f = make_q(T, y)
    starts = [] if seeds is None else list(seeds)
    starts += [np.concatenate([[rng.uniform(0.3, 0.5)], rng.uniform(0, 2 * np.pi, T.K)]) for _ in range(nrand)]
    best = []
    for z0 in starts:
        res = minimize(f, z0, jac=True, method="L-BFGS-B", bounds=[(0.3, 0.5)] + [(None, None)] * T.K,
                       options={"maxiter": 400, "ftol": 1e-14, "gtol": 1e-10})
        best.append((-res.fun, res.x))
    best.sort(key=lambda t: -t[0])
    return best[0][0], [b[1] for b in best[:40]]


def min_outer(T, y, rng):
    """min over sigma in [0.3,0.5] and w of |P(1+sigma,w)| - r(1+sigma)/y (must be > 0 for B_n finite)."""
    absA = np.abs(T.alpha)

    def f(z):
        x, w = 1 + z[0], z[1:]
        P, D, dP, _ = T.PP_grad(x, w)
        e = np.exp(-2 * x * T.logk)
        r = np.sum(absA * e)
        rp = np.sum(2 * T.logk * absA * e)
        A = abs(P)
        g = np.zeros_like(z)
        g[1:] = np.real(np.conj(P) * dP[1:]) / A
        g[0] = np.real(np.conj(P) * D) / A + rp / y
        return A - r / y, g
    vals = []
    for _ in range(300):
        z0 = np.concatenate([[rng.uniform(0.3, 0.5)], rng.uniform(0, 2 * np.pi, T.K)])
        res = minimize(f, z0, jac=True, method="L-BFGS-B", bounds=[(0.3, 0.5)] + [(None, None)] * T.K)
        vals.append(res.fun)
    return min(vals)


def work(n):
    ex = {(r["kind"], tuple(r["band"])): r for r in EX if r["n"] == n}
    T = Torus(n)
    rng = np.random.default_rng(7 * n)
    out = {"n": n}
    if n in PRINTED:
        m, M, l = (ex[(k, (1.2, 1.3))]["value"] for k in ("m", "M", "l"))
        r, rp = T.r(1.2), T.rp(1.2)
        pm, pM, pth = PRINTED[n]
        out.update({"m": m, "M": M, "l": l, "r12": r, "rp12": rp, "printed": [pm, pM, pth]})
        out["y1_unrounded"] = root1(m, M, l, r, rp)
        out["y1_printed_mM"] = root1(pm, pM, l, r, rp)
        out["y1_printed_m_unrounded_M"] = root1(pm, M, l, r, rp)
        out["y1_unrounded_l0"] = root1(m, M, 0.0, r, rp)
        # second condition
        lo = out["y1_unrounded"] * 0.5
        # find y where B finite
        def g2(y, seeds):
            if min_outer(T, y, rng) <= 0:
                return -1e9, seeds
            b, s = Bn(T, y, rng, seeds)
            return math.log(y / math.pi) - 1 / y - b, s
        seeds = None
        ys = []
        y = max(pth, 5.0)
        val, seeds = g2(y, seeds)
        ys.append((y, val))
        if val > 0:
            hi = y
            lo_ = y / 1.5
            while True:
                v, seeds = g2(lo_, seeds)
                ys.append((lo_, v))
                if v <= 0:
                    break
                hi = lo_
                lo_ /= 1.5
        else:
            lo_ = y
            hi = y * 1.5
            while True:
                v, seeds = g2(hi, seeds)
                ys.append((hi, v))
                if v > 0:
                    break
                lo_ = hi
                hi *= 1.5
        for _ in range(30):
            mid = 0.5 * (lo_ + hi)
            v, seeds = g2(mid, seeds)
            ys.append((mid, v))
            if v > 0:
                hi = mid
            else:
                lo_ = mid
            if hi - lo_ < 1e-3:
                break
        out["y2"] = hi
        out["B_at_printed_threshold"] = Bn(T, float(pth), rng, seeds)[0]
        out["rhs_at_printed_threshold"] = math.log(pth / math.pi) - 1 / pth
    if n in SINGLE:
        m, M, l = (ex[(k, (1.0, 1.5))]["value"] for k in ("m", "M", "l"))
        r, rp = T.r(1.0), T.rp(1.0)
        out["single"] = {"m": m, "M": M, "l": l, "r1": r, "rp1": rp, "printed_threshold": SINGLE[n],
                         "y1": root1(m, M, l, r, rp)}
    return out


EX = json.load(open("explore.json"))
if __name__ == "__main__":
    ns = sorted(set(PRINTED) | set(SINGLE))
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 24) as pool:
        res = pool.map(work, ns, chunksize=1)
    json.dump(res, open("thresholds.json", "w"), indent=1)
    for o in res:
        n = o["n"]
        if "printed" in o:
            print(f"n={n}: m={o['m']:.5f} M={o['M']:.5f} l={o['l']:.5f} r(1.2)={o['r12']:.4f} r'(1.2)={o['rp12']:.4f} printed={o['printed']}")
            print(f"   y1 unrounded={o['y1_unrounded']:.4f}  y1 printed m,M={o['y1_printed_mM']:.4f}  y1 printed m,unrounded M={o['y1_printed_m_unrounded_M']:.4f}  y1(l=0)={o['y1_unrounded_l0']:.3f}")
            print(f"   y2 (second condition)={o['y2']:.3f}  B(printed threshold)={o['B_at_printed_threshold']:.4f} vs log(y/pi)-1/y={o['rhs_at_printed_threshold']:.4f}")
        if "single" in o:
            s = o["single"]
            print(f"n={n} single band: m={s['m']:.5f} M={s['M']:.5f} l={s['l']:.5f} r(1)={s['r1']:.4f} r'(1)={s['rp1']:.4f} y1={s['y1']:.4f} printed={s['printed_threshold']}")
