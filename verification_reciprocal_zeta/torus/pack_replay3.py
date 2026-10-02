"""Task lists of the replay (bnb7.py) of the pieces mlow 2298 and mlow 2426 (arguments: the pieces, comma separated, and
the output file; tasks_replay3.json for 2298, run N3-replay3, and tasks_replay4.json for 2426, run N3-replay4), whose sub
boxes follow the top of their own search tree (make_replay3.py, make_replay3b.py).

Piece 2298: the serial run of make_replay3.py (N3-mr3, task 0). Piece 2426: the serial run (N3-mr3, task 1) if it
finished; otherwise the parallel runs of make_replay3b.py: the boxes cleared above the frontier (N3-mr3b-front), the
sub boxes counted below each frontier box (N3-mr3b), and, if a frontier box was refined (N3-mr3c-front, N3-mr3c)
instead of counted, the sub boxes of the refinement. For each piece the exact volumes of the sub boxes must add up
to the piece (integer arithmetic). The sub boxes are packed into tasks of at most LMAX counted final boxes."""
import json, os

LMAX = 500000
D = "/scratch/vsokolov/rh_book_repro/agents"


def load(p):
    return json.load(open(p))


def volume_ok(items):
    Ds = [sum(d for _, _, d in e["cuts"]) for e in items]
    Dm = max(Ds)
    return sum(1 << (Dm - x) for x in Ds) == 1 << Dm


import sys
want = [int(a) for a in sys.argv[1].split(",")]
out = sys.argv[2]
pieces, source = {}, {}
if 2298 in want:
    pieces[2298], source[2298] = load(f"{D}/N3-mr3/items_0.json"), "N3-mr3 task 0"
if 2426 not in want:
    pass
elif os.path.exists(f"{D}/N3-mr3/items_1.json"):
    pieces[2426], source[2426] = load(f"{D}/N3-mr3/items_1.json"), "N3-mr3 task 1"
else:
    fr = load(f"{D}/N3-mr3b-front/frontier.json")
    nb = len(fr["boxes"])
    its = load(f"{D}/N3-mr3b-front/items_cleared.json")
    missing = [t for t in range(nb) if not os.path.exists(f"{D}/N3-mr3b/items_{t}.json")]
    refined = []
    if missing:
        f2 = load(f"{D}/N3-mr3c-front/frontier2.json")
        refined = f2["refines"]["tids"]
        assert set(missing) <= set(refined), (missing, refined)
        tasks2 = load(f"{D}/N3-mr3c/frontier2_tasks.json")
        its += load(f"{D}/N3-mr3c-front/items_cleared2.json")
        for t in range(len(tasks2)):
            its += load(f"{D}/N3-mr3c/items_{t}.json")
    for t in range(nb):
        if t not in refined:
            its += load(f"{D}/N3-mr3b/items_{t}.json")
    pieces[2426] = its
    source[2426] = "N3-mr3b" + (f", boxes {refined} refined in N3-mr3c" if refined else "")
items = []
for idx, its in pieces.items():
    assert all(e["index"] == idx for e in its)
    ok = volume_ok(its)
    print(json.dumps({"index": idx, "source": source[idx], "sub_boxes": len(its), "counted_leaves": sum(e["est"] for e in its),
                      "max": max(e["est"] for e in its), "volume_ok": ok}))
    assert ok
    items += its
items.sort(key=lambda e: -e["est"])
tasks, cur, load_ = [], [], 0
for e in items:
    if cur and load_ + e["est"] > LMAX:
        tasks.append(cur)
        cur, load_ = [], 0
    cur.append(e)
    load_ += e["est"]
if cur:
    tasks.append(cur)
json.dump(tasks, open(out, "w"))
print(len(items), "sub boxes in", len(tasks), "tasks")
