"""Branch and bound pieces for the ratio kinds M and l with the enclosure of bnb_m.EncM (uses a certified lower bound
m_low of |P_n| on the band, from the kind m runs). Otherwise as bnb2.py: depth first, margin SLACK = 1e-9, bisection
along the largest contribution. Keeps the 40 leaves of smallest float margin and 60 random leaves and checks them
with the Arb version of the same enclosure (arb_disc.ArbDisc.bound_ratio_mlow) and with the independent naive Arb
enclosure with the same lower bound (arb_leafcheck.check2_mlow)."""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, json, time, math
import numpy as np
from torus_common import Torus
from bnb_m import EncM
from bnb2 import dfs, LeafKeeper


def main():
    tid = int(os.environ.get("SLURM_ARRAY_TASK_ID", sys.argv[1]))
    tasks = json.load(open(sys.argv[2]))
    t = tasks[tid]
    n, kind = t["n"], t["kind"]
    T = Torus(n)
    E = EncM(T, t["m_low"])
    lo, hi = np.array(t["lo"], float), np.array(t["hi"], float)
    C0, H0 = (0.5 * (lo + hi))[None, :], (0.5 * (hi - lo))[None, :]
    keeper = LeafKeeper(nsmall=40, nrand=60, seed=tid + 101)
    tgt = t["target"] ** 2 if kind == "l" else t["target"]
    res = dfs(lambda C, H: E.bound(kind, C, H), False, tgt, C0, H0, time_limit=t["time_limit"], keeper=keeper,
              label=f"n={n} {kind} {t.get('label', '')}")
    if res.get("best") is not None and kind == "l":
        res["best"] = math.sqrt(res["best"])
    res["task"] = t
    if t.get("check", True) and keeper.small is not None and res["status"] != "violated":
        from arb_disc import ArbDisc
        from arb_leafcheck import ArbTorus, check2_mlow
        AD, AT = ArbDisc(n), ArbTorus(n)
        leafs = [(keeper.small[1][i], keeper.small[2][i], "small") for i in range(keeper.small[1].shape[0])]
        leafs += [(c, h, "random") for c, h in keeper.rand]
        ok_m = ok_n = 0
        maxdiff = 0.0
        ta = time.time()
        for c, h, cls in leafs:
            c = [float(v) for v in c]
            h = [float(v) for v in h]
            fb = float(E.bound(kind, np.array([c]), np.array([h]))[1][0])
            okm, bm = AD.certify_mlow(kind, t["target"], t["m_low"], c, h)
            ok_m += okm
            if bm is not None:
                maxdiff = max(maxdiff, abs(bm - fb))
            okn, _ = check2_mlow(AT, kind, t["target"], c, h, t["m_low"], budget=t.get("budget", 1500))
            ok_n += okn
        res["arb_check"] = {"leaves": len(leafs), "mirror_ok": ok_m, "naive_ok": ok_n, "max_bound_diff": maxdiff,
                            "seconds": time.time() - ta}
    json.dump(res, open(f"result_{tid}.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k not in ("task", "best_pt")}), flush=True)


if __name__ == "__main__":
    main()
