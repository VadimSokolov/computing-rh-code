"""Second condition of Theorem rz:thm:torus through G (bnb_g.py): floating point search and ball arithmetic proof.

A task is a q piece of tasks_full.json (sigma in [0.3, 0.5], w_2 and w_3 in quarters, the other phases in [0, 2 pi])
or a dyadic sub box of it (cuts [[dim, k, d], ...]: along dim the k-th of 2^d equal parts). The double precision
depth first search (bnb2.dfs with bnb_g.EncG, margin 1e-9) keeps every final box; with "arb" set, every final box is
then proved in Arb (arb_g.ArbG.bound_g > 0) on the exact dyadic box of the exact piece, as in bnb7.py (sigma bounds
3/10 and 1/2, phase bounds rational multiples of pi), with exact bisection in Arb for a box that fails, and the
checks that the dyadic positions are consistent and distinct and that the exact volumes add up to the sub box."""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, json, time
import numpy as np
from flint import arb, fmpq, ctx
from torus_common import Torus
from bnb_g import EncG
from bnb2 import dfs
from bnb7 import AllLeaves, exact_bound, to_arb

ctx.prec = 128


class VolLeaves(AllLeaves):
    def __init__(self, vol0, label, keep):
        super().__init__()
        self.vol0, self.vol, self.label, self.keep, self.last, self.count = vol0, 0.0, label, keep, time.time(), 0

    def add(self, C, H, margin):
        if self.keep:
            super().add(C, H, margin)
        self.count += C.shape[0]
        self.vol += float(np.prod(2 * H, axis=1).sum())
        now = time.time()
        if now - self.last > 120:
            self.last = now
            print(f"  {self.label}: leaves {self.count} discarded volume fraction {self.vol / self.vol0:.6f}", flush=True)


def main():
    tid = int(os.environ.get("SLURM_ARRAY_TASK_ID", sys.argv[1]))
    items = json.load(open(sys.argv[2]))[tid]
    if isinstance(items, dict):
        items = [items]
    full = json.load(open("tasks_full.json"))
    out = []
    for it in items:
        t0 = time.time()
        t = full[it["q_index"]]
        assert t["kind"] == "q"
        n, y = t["n"], t["y"]
        T = Torus(n)
        dim = len(t["lo"])
        sk, sd = [0] * dim, [0] * dim
        for j, k, d in it.get("cuts", []):
            sk[j], sd[j] = k, d
        flo = [t["lo"][i] + (t["hi"][i] - t["lo"][i]) * sk[i] / 2 ** sd[i] for i in range(dim)]
        fhi = [t["lo"][i] + (t["hi"][i] - t["lo"][i]) * (sk[i] + 1) / 2 ** sd[i] for i in range(dim)]
        lo, hi = np.array(flo), np.array(fhi)
        C0, H0 = (0.5 * (lo + hi))[None, :], (0.5 * (hi - lo))[None, :]
        label = f"n={n} G piece {t['piece']} cuts {it.get('cuts', [])}"
        keep = VolLeaves(float(np.prod(hi - lo)), label, it.get("arb", True))
        E = EncG(T, y)
        res = dfs(E.bound, True, 0.0, C0, H0, time_limit=it.get("time_limit"), keeper=keep, label=label)
        rec = {"q_index": it["q_index"], "n": n, "y": y, "piece": t["piece"], "cuts": it.get("cuts", []),
               "float_status": res["status"], "processed": res.get("processed"), "float_leaves": keep.count,
               "min_centre_G": res.get("best"), "float_seconds": time.time() - t0,
               "discarded_volume_fraction": keep.vol / keep.vol0}
        if res["status"] == "violated":
            rec["at"], rec["value"] = res.get("at"), res.get("value")
        if res["status"] != "certified" or not it.get("arb", True):
            rec["status"] = "float " + res["status"]
            out.append(rec)
            print(json.dumps(rec), flush=True)
            continue
        from arb_g import ArbG
        AG = ArbG(n)
        elo = [to_arb(exact_bound(v, i)) for i, v in enumerate(t["lo"])]
        ehi = [to_arb(exact_bound(v, i)) for i, v in enumerate(t["hi"])]
        eL = [b - a for a, b in zip(elo, ehi)]
        Cs, Hs = np.concatenate(keep.C), np.concatenate(keep.H)
        Ls = hi - lo
        dd = np.rint(np.log2(Ls[None, :] / (2 * Hs))).astype(np.int64)
        kk = np.rint((Cs - Hs - lo[None, :]) / (2 * Hs)).astype(np.int64)
        ok_d = bool(np.all(np.abs(2 * Hs * np.exp2(dd) - Ls[None, :]) <= 1e-9 * Ls[None, :]))
        ok_k = bool(np.all(np.abs((Cs - Hs - lo[None, :]) / (2 * Hs) - kk) <= 1e-6) and np.all(kk >= 0) and np.all(kk < 2 ** dd))
        Dsum = dd.sum(1)
        Dmax = int(Dsum.max())
        vol_ok = sum(1 << (Dmax - int(x)) for x in Dsum) == (1 << Dmax)
        KD = np.ascontiguousarray(np.concatenate([kk, dd], 1))
        uniq = np.unique(KD.view(np.dtype((np.void, KD.dtype.itemsize * KD.shape[1])))).shape[0] == KD.shape[0]
        ta = time.time()
        fails = rescued = fb = 0
        for q in range(Cs.shape[0]):
            c, h = [], []
            for i in range(dim):
                K = sk[i] * (1 << int(dd[q, i])) + int(kk[q, i])
                Dd = sd[i] + int(dd[q, i])
                c.append(elo[i] + eL[i] * arb(fmpq(2 * K + 1, 1 << (Dd + 1))))
                h.append(eL[i] * arb(fmpq(1, 1 << (Dd + 1))))
            if AG.bound_g(c, h, y) > 0:
                continue
            stack, used, good = [(c, h)], 0, True
            while stack:
                cc, hh = stack.pop()
                used += 1
                if used > 256:
                    good = False
                    break
                if AG.bound_g(cc, hh, y) > 0:
                    continue
                j = max(range(1, dim), key=lambda i: float(hh[i].mid()))
                if float(hh[0].mid()) * 30 > float(hh[j].mid()):
                    j = 0
                h2 = list(hh)
                h2[j] = hh[j] / 2
                c1, c2 = list(cc), list(cc)
                c1[j] = cc[j] - h2[j]
                c2[j] = cc[j] + h2[j]
                stack.append((c1, h2))
                stack.append((c2, list(h2)))
            fb += used
            rescued += good
            fails += not good
        rec.update({"leaves": int(Cs.shape[0]), "arb_failed": fails, "arb_rescued": rescued, "fallback_boxes": fb,
                    "dyadic_ok": ok_d and ok_k, "volume_ok": bool(vol_ok), "distinct": bool(uniq),
                    "arb_seconds": time.time() - ta})
        rec["status"] = "proved" if (fails == 0 and ok_d and ok_k and vol_ok and uniq) else "not proved"
        out.append(rec)
        print(json.dumps(rec), flush=True)
    json.dump(out, open(f"result_{tid}.json", "w"))


if __name__ == "__main__":
    main()
