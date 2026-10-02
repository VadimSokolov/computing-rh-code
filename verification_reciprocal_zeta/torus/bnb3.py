"""Leaf verification of the floating point branch and bound (bnb2.dfs with the enclosures of bnb.py and bnb_q.py).

One task = one sub box (usually a small box around the extreme point found in a piece of the full run, where the
float margins are smallest, or a random sub box). The float branch and bound is rerun on the sub box, keeping the 40
leaves with the smallest float margin and 60 random leaves. Every kept leaf is then checked three ways:
  mirror : the same disc enclosure in Arb at 128 bits (arb_disc.ArbDisc.certify), exact decimal target;
  naive  : an independent Arb enclosure (interval evaluation of the partial derivatives, no hand made radius
           formulas), bisecting along the largest contribution, with a box budget (arb_leafcheck.check2 for m, M, l;
           arb_disc.ArbQNaive.check2 for q);
  points : the function at 300 random points of the leaf, in floating point, against the target.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, time, json, math
import numpy as np
from torus_common import Torus
from bnb import Enc
from bnb_q import EncQ
from bnb2 import dfs, LeafKeeper
from arb_disc import ArbDisc, ArbQNaive


def point_values(T, kind, C, H, rng, y=None, npts=300):
    Z = C + H * rng.uniform(-1, 1, (npts, C.shape[0]))
    if kind == "q":
        absA = np.abs(T.alpha)
        s = Z[:, 0]
        w = Z[:, 1:]
        x1, x2 = 1.5 - s, 1 + s
        P1, _ = T.PP(x1, w)
        P2, _ = T.PP(x2, w)
        r1 = (absA * np.exp(-2 * np.outer(x1, T.logk))).sum(1)
        r2 = (absA * np.exp(-2 * np.outer(x2, T.logk))).sum(1)
        a2 = np.abs(P2) - r2 / y
        with np.errstate(invalid="ignore", divide="ignore"):
            v = (np.log(np.abs(P1) + r1 / y) - np.log(a2)) / (2 * s - 0.5)
        v = np.where(a2 > 0, v, np.inf)
        return v.max()
    P, P1 = T.PP(Z[:, 0], Z[:, 1:])
    if kind == "m":
        return np.abs(P).min()
    if kind == "M":
        return (-P1 / P).real.max()
    return np.abs(P1 / P).max()


def main():
    tid = int(os.environ.get("SLURM_ARRAY_TASK_ID", sys.argv[1]))
    tasks = json.load(open(sys.argv[2]))
    budget = int(sys.argv[3]) if len(sys.argv) > 3 else 1500
    t = tasks[tid]
    n, kind = t["n"], t["kind"]
    T = Torus(n)
    lo, hi = np.array(t["lo"], float), np.array(t["hi"], float)
    C0, H0 = (0.5 * (lo + hi))[None, :], (0.5 * (hi - lo))[None, :]
    keeper = LeafKeeper(nsmall=40, nrand=60, seed=tid + 11)
    label = f"n={n} {kind} {t.get('label', '')}"
    if kind == "q":
        E = EncQ(T, t["y"])
        res = dfs(E.bound, False, t["target"], C0, H0, time_limit=t["time_limit"], keeper=keeper, label=label)
    else:
        E = Enc(T)
        tgt = t["target"] ** 2 if kind in ("m", "l") else t["target"]
        res = dfs(lambda C, H: E.bound(kind, C, H), kind == "m", tgt, C0, H0, time_limit=t["time_limit"],
                  keeper=keeper, label=label)
    res["task"] = t
    print("float:", {k: v for k, v in res.items() if k not in ("task", "best_pt", "at")}, flush=True)
    if res["status"] in ("violated", "B infinite") or keeper.small is None:
        json.dump(res, open(f"result_{tid}.json", "w"), indent=1)
        return
    AD = ArbDisc(n)
    if kind == "q":
        NV = ArbQNaive(n)
    else:
        from arb_leafcheck import ArbTorus, check2
        AT = ArbTorus(n)
    rng = np.random.default_rng(tid)
    small = [(keeper.small[1][i], keeper.small[2][i], float(keeper.small[0][i]), "small")
             for i in range(keeper.small[1].shape[0])]
    rand = [(c, h, None, "random") for c, h in keeper.rand]
    rows = []
    ta = time.time()
    for c, h, fm, cls in small + rand:
        c = [float(v) for v in c]
        h = [float(v) for v in h]
        if kind == "q":
            fb = float(E.bound(np.array([c]), np.array([h]))[1][0])
        else:
            fb = float(E.bound(kind, np.array([c]), np.array([h]))[1][0])
        ok_m, b_m = AD.certify(kind, t["target"], c, h, y=t.get("y"))
        if kind == "q":
            ok_n, used = NV.check2(c, h, t["y"], budget=budget)
        else:
            ok_n, used = check2(AT, kind, t["target"], c, h, budget=budget)
        ext = point_values(T, kind, np.array(c), np.array(h), rng, y=t.get("y"))
        if kind == "m":
            viol = bool(ext < t["target"])
        elif kind == "q":
            viol = bool(ext >= math.log(t["y"] / math.pi) - 1 / t["y"])
        else:
            viol = bool(ext > t["target"])
        rows.append({"cls": cls, "float_margin": fm, "c": c, "h": h, "float_bound": fb, "mirror_ok": ok_m, "mirror_bound": b_m,
                     "bound_diff": None if (b_m is None or not math.isfinite(fb)) else abs(fb - b_m),
                     "naive_ok": ok_n, "naive_boxes": used, "point_extreme": float(ext), "point_violation": viol})
    summ = {}
    for cls in ("small", "random"):
        rr = [r for r in rows if r["cls"] == cls]
        summ[cls] = {"leaves": len(rr), "mirror_ok": sum(r["mirror_ok"] for r in rr),
                     "naive_ok": sum(r["naive_ok"] for r in rr), "point_violations": sum(r["point_violation"] for r in rr),
                     "min_float_margin": min((r["float_margin"] for r in rr if r["float_margin"] is not None), default=None),
                     "max_bound_diff": max((r["bound_diff"] for r in rr if r["bound_diff"] is not None), default=None)}
    res["leaf_checks"] = summ
    res["leaf_rows"] = rows
    res["check_seconds"] = time.time() - ta
    json.dump(res, open(f"result_{tid}.json", "w"), indent=1)
    print("checks:", json.dumps(summ), f"{time.time()-ta:.0f}s", flush=True)


if __name__ == "__main__":
    main()
