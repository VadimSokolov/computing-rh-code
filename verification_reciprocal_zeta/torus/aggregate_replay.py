"""Aggregate the ball arithmetic replay (bnb7.py) per group (n, kind, band).

Runs, in /scratch/vsokolov/rh_book_repro/agents: N3-replay (tasks_replay.json, the grid sub boxes of make_replay.py),
N3-replay2 (the tasks of N3-replay that ran out of memory, rerun with more), N3-replay3 and N3-replay4 (tasks_replay3.json and
tasks_replay4.json: the pieces mlow 2298 and mlow 2426, whose grid sub boxes were too large, cut instead along the top
of their own search tree by make_replay3.py and make_replay3b.py, packed by pack_replay3.py; these sub boxes replace
every grid sub box of the two pieces).

A group is proved when every one of its pieces is covered by proved sub boxes: for each piece the sub boxes are
dyadic boxes of the piece, the leaves of a bisection tree (their exact volumes add up to 1, checked here in integer
arithmetic), and each sub box must have status "proved" (Arb proof of every final box, dyadic positions consistent
and distinct, volumes adding up)."""
import json, glob
from collections import defaultdict

D = "/scratch/vsokolov/rh_book_repro/agents"
src = {"full": json.load(open(f"{D}/N3-full/tasks_full.json")), "mlow": json.load(open(f"{D}/N3-mlow/tasks_mlow.json"))}
pieces = defaultdict(list)
for items in json.load(open(f"{D}/N3-replay/tasks_replay.json")):
    for it in items:
        pieces[(it["src"], it["index"])].append(("grid", it))
tree = defaultdict(list)
for run, tf in (("N3-replay3", "tasks_replay3.json"), ("N3-replay4", "tasks_replay4.json")):
    for items in json.load(open(f"{D}/{run}/{tf}")):
        for it in items:
            tree[(it["src"], it["index"])].append(("tree", it))
for key, its in tree.items():
    pieces[key] = its
res = {}
for run, dirs in (("grid", ("N3-replay", "N3-replay2")), ("tree", ("N3-replay3", "N3-replay4"))):
    for d in dirs:
        for f in glob.glob(f"{D}/{d}/result_*.json"):
            for r in json.load(open(f)):
                res[(run, r["src"], r["index"], json.dumps(r["cuts"]))] = r
groups = defaultdict(lambda: {"pieces": 0, "pieces_proved": 0, "sub_boxes": 0, "leaves": 0, "arb_failed": 0,
                              "arb_rescued": 0, "seconds": 0.0, "missing": 0, "not_proved": 0, "cover_ok": True})
for (s, i), its in pieces.items():
    t = src[s][i]
    g = groups[(t["n"], t["kind"], tuple(t["band"]))]
    g["pieces"] += 1
    # the sub boxes of a piece are products of dyadic intervals; their volumes must add up to 1
    Ds = [sum(d for _, _, d in it["cuts"]) for _, it in its]
    Dm = max(Ds)
    cover = sum(1 << (Dm - x) for x in Ds) == (1 << Dm)
    ok = cover
    for run, it in its:
        g["sub_boxes"] += 1
        r = res.get((run, s, i, json.dumps(it["cuts"])))
        if r is None:
            g["missing"] += 1
            ok = False
            continue
        g["leaves"] += r.get("leaves", 0)
        g["arb_failed"] += r.get("arb_failed", 0)
        g["arb_rescued"] += r.get("arb_rescued", 0)
        g["seconds"] += r.get("seconds", 0.0)
        if r["status"] != "proved":
            g["not_proved"] += 1
            ok = False
    g["cover_ok"] &= cover
    g["pieces_proved"] += ok
out = []
for k in sorted(groups):
    g = groups[k]
    g["seconds"] = round(g["seconds"])
    row = {"n": k[0], "kind": k[1], "band": list(k[2]), "proved": g["pieces_proved"] == g["pieces"]}
    row.update(g)
    out.append(row)
    print(json.dumps(row))
tot = {x: sum(r[x] for r in out) for x in ("pieces", "pieces_proved", "sub_boxes", "leaves", "arb_failed", "arb_rescued", "seconds", "missing", "not_proved")}
tot["groups"] = len(out)
tot["groups_proved"] = sum(r["proved"] for r in out)
tot["tree_pieces"] = sorted(f"{s} {i}" for s, i in tree)
print("TOTAL", json.dumps(tot))
json.dump({"groups": out, "total": tot}, open("replay_summary.json", "w"), indent=1)
