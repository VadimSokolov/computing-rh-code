"""Parallel version of make_replay3.py (the sub boxes of the replay along the top of the search tree), for one piece.

Argument "frontier INDEX": starting from the piece INDEX of tasks_mlow.json, split boxes by the rule of bnb2.dfs
(halve the coordinate with the largest term sup|d_j F| h_j), breadth first, until at least F boxes remain to be split;
a box cleared on the way (bound within the target by SLACK) is a final box of the tree and a sub box of its own.
Writes frontier.json (the boxes left, as [k, d] lists) and items_cleared.json.
Array task tid (argument: the task file, frontier.json): the recursion of make_replay3.py below frontier box tid,
splitting by the same rule until the search from a box, counted, has at most LCAP final boxes; writes items_<tid>.json.
Argument "refine FILE SIZE TID...": replace the boxes TID of the frontier file FILE by frontiers of SIZE boxes (frontier2.json,
items_cleared2.json), for boxes too heavy for one task; the array then runs on frontier2.json.
Together the sub boxes are the leaves of one bisection tree of the piece, the same tree as the search from the piece."""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, json, time
import numpy as np
from torus_common import Torus
from bnb import SLACK
from bnb_m import EncM
from bnb2 import dfs

F = 256
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


class Piece:
    def __init__(self, index):
        self.index = index
        t = json.load(open("tasks_mlow.json"))[index]
        self.n, self.kind = t["n"], t["kind"]
        assert self.kind in ("M", "l")
        self.E = EncM(Torus(self.n), t["m_low"])
        self.tgt = t["target"] ** 2 if self.kind == "l" else t["target"]
        self.lo, self.hi = t["lo"], t["hi"]
        self.dim = len(self.lo)

    def bound(self, C, H):
        return self.E.bound(self.kind, C, H)

    def box(self, kk, dd):
        lo, hi = self.lo, self.hi
        flo = [lo[i] + (hi[i] - lo[i]) * kk[i] / 2 ** dd[i] for i in range(self.dim)]
        fhi = [lo[i] + (hi[i] - lo[i]) * (kk[i] + 1) / 2 ** dd[i] for i in range(self.dim)]
        a, b = np.array(flo), np.array(fhi)
        return (0.5 * (a + b))[None, :], (0.5 * (b - a))[None, :]

    def item(self, kk, dd, est):
        return {"src": "mlow", "index": self.index, "n": self.n, "kind": self.kind,
                "cuts": [[i, kk[i], dd[i]] for i in range(self.dim) if dd[i] > 0], "est": est}

    def split(self, kk, dd):
        """None if the box is cleared, else its two children by the rule of bnb2.dfs."""
        C0, H0 = self.box(kk, dd)
        out = self.bound(C0, H0)
        Fv, bd, contrib = out[0], out[1], out[2]
        assert Fv[0] <= self.tgt, ("violated", C0.tolist())
        if bd[0] <= self.tgt - SLACK:
            return None
        cr = contrib[0]
        if np.isfinite(cr).all():
            j = int(np.argmax(cr))
        else:
            sc = H0[0].copy()
            sc[0] *= 30.0
            j = int(np.argmax(sc))
        kids = []
        for b in (0, 1):
            k2, d2 = list(kk), list(dd)
            k2[j], d2[j] = 2 * kk[j] + b, dd[j] + 1
            kids.append((k2, d2))
        return kids


def expand(P, start, size=F):
    front, cleared = [start], []
    while front and len(front) < size:
        new = []
        for kk, dd in front:
            kids = P.split(kk, dd)
            if kids is None:
                cleared.append(P.item(kk, dd, 1))
            else:
                new.extend(kids)
        front = new
    return front, cleared


def frontier(index):
    P = Piece(index)
    front, cleared = expand(P, ([0] * P.dim, [0] * P.dim))
    json.dump({"index": index, "boxes": front}, open("frontier.json", "w"))
    json.dump(cleared, open("items_cleared.json", "w"))
    print(json.dumps({"index": index, "frontier": len(front), "cleared": len(cleared)}))


def refine(path, size, tids):
    """Replace the boxes tids of the frontier file path by their own frontiers (frontier2.json, items_cleared2.json)."""
    fr = json.load(open(path))
    P = Piece(fr["index"])
    front, cleared = [], []
    for tid in tids:
        f, c = expand(P, tuple(fr["boxes"][tid]), size)
        front += f
        cleared += c
    json.dump({"index": fr["index"], "boxes": front, "refines": {"file": path, "tids": tids}}, open("frontier2.json", "w"))
    json.dump(cleared, open("items_cleared2.json", "w"))
    print(json.dumps({"index": fr["index"], "refined": tids, "frontier": len(front), "cleared": len(cleared)}))


def count(tid, path, boxes=None):
    """The recursion below the frontier boxes 'boxes' (default: box tid) of the frontier file path."""
    fr = json.load(open(path))
    P = Piece(fr["index"])
    t0 = time.time()
    queue = [tuple(fr["boxes"][b]) for b in (boxes if boxes is not None else [tid])]
    items, aborted = [], 0
    while queue:
        kk, dd = queue.pop()
        C0, H0 = P.box(kk, dd)
        try:
            res = dfs(P.bound, False, P.tgt, C0, H0, keeper=Cap(LCAP), report_every=1e9)
            assert res["status"] == "certified", res
            items.append(P.item(kk, dd, res["leaves"]))
            continue
        except TooMany:
            aborted += 1
        kids = P.split(kk, dd)
        assert kids is not None
        queue.extend(kids)
    json.dump(items, open(f"items_{tid}.json", "w"))
    print(json.dumps({"tid": tid, "sub_boxes": len(items), "counted_leaves": sum(e["est"] for e in items),
                      "aborted": aborted, "seconds": round(time.time() - t0)}))


if __name__ == "__main__":
    if sys.argv[1] == "frontier":
        frontier(int(sys.argv[2]))
    elif sys.argv[1] == "refine":
        refine(sys.argv[2], int(sys.argv[3]), [int(a) for a in sys.argv[4:]])
    else:
        # array task: argv[2] is the task list (the boxes of the frontier file, one per task)
        tid = int(os.environ.get("SLURM_ARRAY_TASK_ID", sys.argv[1]))
        if sys.argv[2].startswith("frontier2"):
            # task list of frontier2.json: a list of lists of box indices
            count(tid, "frontier2.json", json.load(open(sys.argv[2]))[tid])
        else:
            count(tid, "frontier.json")
