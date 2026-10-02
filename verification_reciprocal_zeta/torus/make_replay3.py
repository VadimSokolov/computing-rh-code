"""Sub boxes for the ball arithmetic replay (bnb7.py) of the two pieces whose grid replay was too large.

The grid of make_replay.py cut the pieces mlow 2298 (n=32, l) and mlow 2426 (n=36, l) along the widest torus
coordinates, which belong to the large primes; the search itself hardly cuts them, so the grid multiplied the number of
final boxes near the maximum of |P'/P| (several sub boxes ran out of 12 GB). Here the sub boxes are the top of the
search tree itself: starting from the piece, a box is split by the rule of bnb2.dfs (halve the coordinate with the
largest term sup|d_j F| h_j) until the search from the box, run in counting mode, has at most LCAP final boxes; a box
cleared on the way is a final box of the tree and a sub box of its own. The sub boxes are thus the leaves of a
bisection tree of the piece, and bnb7.py reproduces below each of them the same search as the certified run. Each
sub box is written as cuts [[dim, k, d], ...] (along dim the k-th of 2^d equal parts of the piece; bnb7.py builds the
float box from these, and so does the counting here). Task argument: an index into [2298, 2426]; output
items_<index>.json with the counted final boxes of every sub box."""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, json, time
import numpy as np
from torus_common import Torus
from bnb_m import EncM
from bnb2 import dfs

PIECES = [2298, 2426]
LCAP = 400000


class TooMany(Exception):
    pass


class Cap:
    def __init__(self, cap):
        self.cap, self.count = cap, 0

    def add(self, C, H, margin):
        self.count += C.shape[0]
        if self.count > self.cap:
            raise TooMany()


def main():
    tid = int(os.environ.get("SLURM_ARRAY_TASK_ID", sys.argv[1]))
    index = PIECES[tid]
    t = json.load(open("tasks_mlow.json"))[index]
    n, kind = t["n"], t["kind"]
    assert kind in ("M", "l")
    E = EncM(Torus(n), t["m_low"])
    tgt = t["target"] ** 2 if kind == "l" else t["target"]
    bound = lambda C, H: E.bound(kind, C, H)
    lo, hi = t["lo"], t["hi"]
    dim = len(lo)

    def box(kk, dd):
        flo = [lo[i] + (hi[i] - lo[i]) * kk[i] / 2 ** dd[i] for i in range(dim)]
        fhi = [lo[i] + (hi[i] - lo[i]) * (kk[i] + 1) / 2 ** dd[i] for i in range(dim)]
        a, b = np.array(flo), np.array(fhi)
        return (0.5 * (a + b))[None, :], (0.5 * (b - a))[None, :]

    t0 = time.time()
    queue = [([0] * dim, [0] * dim)]
    items, aborted = [], 0
    while queue:
        kk, dd = queue.pop()
        C0, H0 = box(kk, dd)
        cap = Cap(LCAP)
        try:
            res = dfs(bound, False, tgt, C0, H0, keeper=cap, report_every=1e9)
            assert res["status"] == "certified", res
            items.append({"src": "mlow", "index": index, "n": n, "kind": kind,
                          "cuts": [[i, kk[i], dd[i]] for i in range(dim) if dd[i] > 0], "est": res["leaves"]})
            continue
        except TooMany:
            aborted += 1
        # split by the rule of bnb2.dfs
        out = bound(C0, H0)
        F, bd, contrib = out[0], out[1], out[2]
        assert F[0] <= tgt, ("violated", C0.tolist())
        cr = contrib[0]
        if np.isfinite(cr).all():
            j = int(np.argmax(cr))
        else:
            sc = H0[0].copy()
            sc[0] *= 30.0
            j = int(np.argmax(sc))
        for b in (0, 1):
            k2, d2 = list(kk), list(dd)
            k2[j], d2[j] = 2 * kk[j] + b, dd[j] + 1
            queue.append((k2, d2))
        if aborted % 20 == 0:
            print(f"  piece {index}: sub boxes {len(items)} counted {sum(e['est'] for e in items)} aborted {aborted} "
                  f"queue {len(queue)} t={time.time() - t0:.0f}s", flush=True)
    # the sub boxes are leaves of a bisection tree: their volumes add up to the piece
    Ds = [sum(d for _, _, d in e["cuts"]) for e in items]
    Dm = max(Ds)
    assert sum(1 << (Dm - x) for x in Ds) == 1 << Dm
    json.dump(items, open(f"items_{tid}.json", "w"))
    print(json.dumps({"index": index, "n": n, "kind": kind, "sub_boxes": len(items), "counted_leaves": sum(e["est"] for e in items),
                      "max_leaves": max(e["est"] for e in items), "aborted": aborted, "seconds": round(time.time() - t0)}))


if __name__ == "__main__":
    main()
