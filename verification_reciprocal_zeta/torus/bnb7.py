"""Ball arithmetic replay of the branch and bound for the kinds m, M and l: every leaf is re-proved in Arb.

For each sub box of a task (see make_replay.py) the floating point branch and bound of the certified run is rerun on
the sub box (bnb.Enc for m, bnb_m.EncM with the certified m_low for M and l, bnb2.dfs), keeping every leaf. Each
leaf, a float box (c, h), is identified with an exact dyadic part of the piece: along coordinate i it is the k-th of
2^d equal parts, with k and d recovered from c and h (the float boxes differ from the exact ones only by rounding).
The exact box is built in Arb from the exact piece: x bounds are the decimals of the band (1, 6/5, 13/10, 3/2), torus
bounds are rational multiples of pi. The Arb enclosure of arb_disc.py (bound for m, bound_ratio_mlow for M and l,
the same formulas as the float code, with every rounding error enclosed) must then prove the target on the exact
box; a leaf that fails is bisected exactly in Arb (budget 256 boxes). Checks per sub box: the recovered (k, d) are
consistent, and the exact volumes of the leaves add up to the volume of the sub box (sum of 2^-D = 1 in integer
arithmetic), so that the leaves, being the leaves of a bisection tree, partition the sub box."""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, json, time, math
from fractions import Fraction as Fr
import numpy as np
from flint import arb, fmpq, ctx
from torus_common import Torus
from bnb import Enc
from bnb_m import EncM
from bnb2 import dfs
from arb_disc import ArbDisc, dec

ctx.prec = 128
PI = arb.pi()


class AllLeaves:
    def __init__(self):
        self.C, self.H = [], []

    def add(self, C, H, margin):
        self.C.append(C.copy())
        self.H.append(H.copy())


def exact_bound(v, dim):
    """Exact value intended by the float bound v of a piece: a decimal for x, a rational multiple of pi for a phase."""
    if dim == 0:
        f = Fr(repr(float(v)))
        return ("q", f)
    f = Fr(v / math.pi).limit_denominator(16)
    assert abs(float(f) * math.pi - v) <= 4e-15, (v, f)
    return ("pi", f)


def to_arb(e):
    kind, f = e
    a = arb(fmpq(f.numerator, f.denominator))
    return a * PI if kind == "pi" else a


def main():
    tid = int(os.environ.get("SLURM_ARRAY_TASK_ID", sys.argv[1]))
    items = json.load(open(sys.argv[2]))[tid]
    src_tasks = {"full": json.load(open("tasks_full.json")), "mlow": json.load(open("tasks_mlow.json"))}
    out = []
    cache = {}
    for it in items:
        t0 = time.time()
        t = src_tasks[it["src"]][it["index"]]
        n, kind = t["n"], t["kind"]
        if n not in cache:
            cache[n] = (Torus(n), ArbDisc(n))
        T, AD = cache[n]
        dim = len(t["lo"])
        # exact piece and the sub box (float and exact)
        elo = [to_arb(exact_bound(v, i)) for i, v in enumerate(t["lo"])]
        ehi = [to_arb(exact_bound(v, i)) for i, v in enumerate(t["hi"])]
        eL = [b - a for a, b in zip(elo, ehi)]
        sk, sd = [0] * dim, [0] * dim
        for j, k, d in it["cuts"]:
            sk[j], sd[j] = k, d
        flo = [t["lo"][i] + (t["hi"][i] - t["lo"][i]) * sk[i] / 2 ** sd[i] for i in range(dim)]
        fhi = [t["lo"][i] + (t["hi"][i] - t["lo"][i]) * (sk[i] + 1) / 2 ** sd[i] for i in range(dim)]
        lo, hi = np.array(flo), np.array(fhi)
        C0, H0 = (0.5 * (lo + hi))[None, :], (0.5 * (hi - lo))[None, :]
        keep = AllLeaves()
        if kind == "m":
            E = Enc(T)
            res = dfs(lambda C, H: E.bound("m", C, H), True, t["target"] ** 2, C0, H0, keeper=keep, report_every=1e9)
        else:
            E = EncM(T, t["m_low"])
            tgt = t["target"] ** 2 if kind == "l" else t["target"]
            res = dfs(lambda C, H: E.bound(kind, C, H), False, tgt, C0, H0, keeper=keep, report_every=1e9)
        rec = {"src": it["src"], "index": it["index"], "n": n, "kind": kind, "cuts": it["cuts"], "float_status": res["status"]}
        if res["status"] != "certified":
            rec["status"] = "float run not certified"
            out.append(rec)
            continue
        Cs, Hs = np.concatenate(keep.C), np.concatenate(keep.H)
        Ls = hi - lo
        # recover the dyadic position of every leaf inside the sub box
        dd = np.rint(np.log2(Ls[None, :] / (2 * Hs))).astype(np.int64)
        kk = np.rint((Cs - Hs - lo[None, :]) / (2 * Hs)).astype(np.int64)
        ok_d = np.all(np.abs(2 * Hs * np.exp2(dd) - Ls[None, :]) <= 1e-9 * Ls[None, :])
        ok_k = np.all(np.abs((Cs - Hs - lo[None, :]) / (2 * Hs) - kk) <= 1e-6) and np.all(kk >= 0) and np.all(kk < 2 ** dd)
        Dsum = dd.sum(1)
        Dmax = int(Dsum.max())
        vol = sum(1 << (Dmax - int(x)) for x in Dsum)
        vol_ok = vol == (1 << Dmax)
        KD = np.ascontiguousarray(np.concatenate([kk, dd], 1))
        uniq = np.unique(KD.view(np.dtype((np.void, KD.dtype.itemsize * KD.shape[1])))).shape[0] == KD.shape[0]
        t_target = dec(t["target"])
        m_low = t.get("m_low")

        def prove(c, h):
            if kind == "m":
                return AD.bound("m", c, h) > t_target * t_target
            b = AD.bound_ratio_mlow(kind, m_low, c, h)
            return b < (t_target if kind == "M" else t_target * t_target)

        fails, fallback_boxes, rescued = 0, 0, 0
        for q in range(Cs.shape[0]):
            c, h = [], []
            for i in range(dim):
                K = sk[i] * (1 << int(dd[q, i])) + int(kk[q, i])
                Dd = sd[i] + int(dd[q, i])
                c.append(elo[i] + eL[i] * arb(fmpq(2 * K + 1, 1 << (Dd + 1))))
                h.append(eL[i] * arb(fmpq(1, 1 << (Dd + 1))))
            if prove(c, h):
                continue
            # exact bisection in Arb
            stack, used, good = [(c, h)], 0, True
            while stack:
                cc, hh = stack.pop()
                used += 1
                if used > 256:
                    good = False
                    break
                if prove(cc, hh):
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
            fallback_boxes += used
            if good:
                rescued += 1
            else:
                fails += 1
        rec.update({"leaves": int(Cs.shape[0]), "arb_failed": fails, "arb_rescued": rescued, "fallback_boxes": fallback_boxes,
                    "dyadic_ok": bool(ok_d and ok_k), "volume_ok": bool(vol_ok), "distinct": bool(uniq),
                    "seconds": time.time() - t0})
        rec["status"] = "proved" if (fails == 0 and ok_d and ok_k and vol_ok and uniq) else "not proved"
        out.append(rec)
        print(json.dumps(rec), flush=True)
    json.dump(out, open(f"result_{tid}.json", "w"))


if __name__ == "__main__":
    main()
