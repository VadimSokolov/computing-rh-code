"""Parallel pieces of the torus branch and bound (numpy floating point), with an Arb spot check of leaves.

One task = one sub box of the domain. Kinds: "m" (prove |P|^2 >= target^2 on band x torus), "M" (prove
-Re(P'/P) <= target), "l" (prove |P'/P|^2 <= target^2), "q" (prove q(sigma,w) <= target for the second condition of
Theorem rz:thm:torus at height y, with |P(1+sigma,w)| - r(1+sigma)/y > 0). The torus is halved by the symmetry
w -> -w (P -> conj P), which leaves |P|, Re(P'/P), |P'/P| and q unchanged: w_2 runs over [0, pi] only.
Enclosures: bnb.Enc and bnb_q.EncQ (mean value form with ball enclosures); a box is discarded when its bound clears
the target by SLACK = 1e-9. Of the discarded boxes, the 40 with the smallest float margin and 120 random ones are
re-checked independently in Arb (arb_leafcheck.py) for the kinds m, M, l.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, time, json, math
import numpy as np
from torus_common import Torus
from bnb import Enc, SLACK
from bnb_q import EncQ


class LeafKeeper:
    def __init__(self, nsmall=40, nrand=120, seed=1):
        self.nsmall, self.nrand = nsmall, nrand
        self.small = None  # (margin, C, H)
        self.rand = []
        self.rng = np.random.default_rng(seed)
        self.seen = 0

    def add(self, C, H, margin):
        # smallest margins
        if self.small is None:
            M, CC, HH = margin, C, H
        else:
            M = np.concatenate([self.small[0], margin])
            CC = np.concatenate([self.small[1], C])
            HH = np.concatenate([self.small[2], H])
        idx = np.argsort(M)[:self.nsmall]
        self.small = (M[idx], CC[idx], HH[idx])
        # random: reservoir over batches (each leaf kept with probability ~ nrand/seen overall)
        for i in range(C.shape[0]):
            self.seen += 1
            if len(self.rand) < self.nrand:
                self.rand.append((C[i].copy(), H[i].copy()))
            else:
                j = self.rng.integers(0, self.seen)
                if j < self.nrand:
                    self.rand[j] = (C[i].copy(), H[i].copy())
            if i > 200:  # cap the per batch Python loop: sample the rest
                rest = C.shape[0] - i - 1
                self.seen += rest
                k = self.rng.binomial(rest, min(1.0, self.nrand / max(self.seen, 1)))
                for _ in range(min(k, self.nrand)):
                    jj = self.rng.integers(i + 1, C.shape[0])
                    self.rand[self.rng.integers(0, self.nrand)] = (C[jj].copy(), H[jj].copy())
                break


def dfs(bound, lower_kind, tgt, C0, H0, batch=100000, time_limit=None, keeper=None, report_every=60.0, label=""):
    stack = [(C0, H0)]
    processed = leaves = 0
    best = None
    best_pt = None
    t0 = time.time()
    last = t0
    while stack:
        C, H = stack.pop()
        if C.shape[0] > batch:
            stack.append((C[batch:], H[batch:]))
            C, H = C[:batch], H[:batch]
        out = bound(C, H)
        F, bd, contrib = out[0], out[1], out[2]
        processed += C.shape[0]
        if len(out) > 3:  # q: centre values of |P(x2)| - r/y must be positive
            if not np.all(out[4] > 0):
                j = int(np.argmin(out[4]))
                return {"status": "B infinite", "at": C[j].tolist(), "processed": processed}
        if lower_kind:
            j = int(np.argmin(F))
            if best is None or F[j] < best:
                best, best_pt = float(F[j]), C[j].copy()
            if F[j] < tgt:
                return {"status": "violated", "at": C[j].tolist(), "value": float(F[j]), "processed": processed}
            done = bd >= tgt + SLACK
            margin = bd - tgt
        else:
            j = int(np.argmax(F))
            if best is None or F[j] > best:
                best, best_pt = float(F[j]), C[j].copy()
            if F[j] > tgt:
                return {"status": "violated", "at": C[j].tolist(), "value": float(F[j]), "processed": processed}
            done = bd <= tgt - SLACK
            margin = tgt - bd
        nd = int(done.sum())
        leaves += nd
        if keeper is not None and nd:
            keeper.add(C[done], H[done], margin[done])
        rest = ~done
        if rest.any():
            Cr, Hr, cr = C[rest], H[rest], contrib[rest]
            inf_rows = ~np.isfinite(cr).all(1)
            cr = np.where(np.isfinite(cr), cr, 0.0)
            jj = np.argmax(cr, axis=1)
            if inf_rows.any():
                sc = Hr.copy()
                sc[:, 0] *= 30.0
                jj = np.where(inf_rows, np.argmax(sc, axis=1), jj)
            rows = np.arange(Cr.shape[0])
            Hn = Hr.copy()
            Hn[rows, jj] *= 0.5
            C1 = Cr.copy()
            C2 = Cr.copy()
            C1[rows, jj] -= Hn[rows, jj]
            C2[rows, jj] += Hn[rows, jj]
            stack.append((np.concatenate([C1, C2]), np.concatenate([Hn, Hn])))
        now = time.time()
        if now - last > report_every:
            last = now
            print(f"  {label}: processed={processed:.3e} leaves={leaves:.3e} pending={sum(c.shape[0] for c, _ in stack)} "
                  f"best={best:.6f} t={now-t0:.0f}s", flush=True)
        if time_limit is not None and now - t0 > time_limit:
            return {"status": "gave up", "processed": processed, "leaves": leaves, "best": best, "seconds": now - t0}
    return {"status": "certified", "processed": processed, "leaves": leaves, "best": best,
            "best_pt": None if best_pt is None else best_pt.tolist(), "seconds": time.time() - t0}


def run_task(t):
    n, kind = t["n"], t["kind"]
    T = Torus(n)
    lo = np.array(t["lo"], dtype=float)
    hi = np.array(t["hi"], dtype=float)
    C0 = (0.5 * (lo + hi))[None, :]
    H0 = (0.5 * (hi - lo))[None, :]
    label = f"n={n} {kind} piece={t.get('piece')}"
    keeper = LeafKeeper(seed=t.get("piece", 0) + 17 * n) if kind in ("m", "M", "l") else None
    if kind == "q":
        E = EncQ(T, t["y"])
        res = dfs(E.bound, False, t["target"], C0, H0, time_limit=t.get("time_limit"), label=label)
    else:
        E = Enc(T)
        tgt = t["target"] ** 2 if kind in ("m", "l") else t["target"]
        res = dfs(lambda C, H: E.bound(kind, C, H), kind == "m", tgt, C0, H0, time_limit=t.get("time_limit"),
                  keeper=keeper, label=label)
        if res["best"] is not None and kind in ("m", "l"):
            res["best"] = math.sqrt(res["best"])
        if res.get("value") is not None and kind in ("m", "l"):
            res["value"] = math.sqrt(res["value"])
    res["task"] = t
    # Arb spot check of kept leaves
    if keeper is not None and res["status"] == "certified" and keeper.small is not None:
        from arb_leafcheck import ArbTorus, check
        AT = ArbTorus(n)
        leafs = [(keeper.small[1][i], keeper.small[2][i]) for i in range(keeper.small[1].shape[0])] + keeper.rand
        ok = bad = used = 0
        ta = time.time()
        for c, h in leafs:
            good, nb = check(AT, kind, t["target"], list(map(float, c)), list(map(float, h)))
            used += nb
            if good:
                ok += 1
            else:
                bad += 1
        res["arb_check"] = {"checked": ok + bad, "ok": ok, "failed": bad, "arb_boxes": used,
                            "smallest_float_margin": float(keeper.small[0][0]), "seconds": time.time() - ta}
    return res


if __name__ == "__main__":
    tid = int(os.environ.get("SLURM_ARRAY_TASK_ID", sys.argv[1] if len(sys.argv) > 1 else 0))
    tasks = json.load(open(sys.argv[2]))
    t = tasks[tid]
    r = run_task(t)
    json.dump(r, open(f"result_{tid}.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in r.items() if k != "task"}))
