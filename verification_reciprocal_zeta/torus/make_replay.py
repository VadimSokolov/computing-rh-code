"""Task list of the ball arithmetic replay (bnb7.py) of the kinds m, M and l, run on Hopper.

Every piece of the certified runs (kind m: tasks_full.json, /scratch/.../N3-full; kinds M and l: tasks_mlow.json,
/scratch/.../N3-mlow) is cut dyadically into 2^s sub boxes so that each sub box has about LMAX leaves or fewer
(estimated from the leaf count of the piece), halving the coordinates of largest width in turn; the sub boxes are
then packed into tasks of at most LMAX estimated leaves. A sub box is written as (source, index, cuts) with cuts a
list of [dim, k, d]: along dim it is the k-th of 2^d equal parts of the piece."""
import json, math

LMAX = 500000
D = "/scratch/vsokolov/rh_book_repro/agents"
items = []
for src, kinds in (("full", ("m",)), ("mlow", ("M", "l"))):
    tasks = json.load(open(f"{D}/N3-{src}/tasks_{src}.json"))
    for i, t in enumerate(tasks):
        if t["kind"] not in kinds:
            continue
        r = json.load(open(f"{D}/N3-{src}/result_{i}.json"))
        assert r["status"] == "certified", (src, i)
        leaves = r["leaves"]
        s = max(0, math.ceil(math.log2(max(leaves, 1) / LMAX)))
        width = [hi - lo for lo, hi in zip(t["lo"], t["hi"])]
        cuts = {}
        for _ in range(s):
            j = max(range(1, len(width)), key=lambda q: width[q])  # the torus coordinates only
            width[j] /= 2
            cuts[j] = cuts.get(j, 0) + 1
        dims = sorted(cuts)
        est = leaves / 2 ** s
        # enumerate the 2^s sub boxes
        for code in range(2 ** s):
            cc, rest = [], code
            for j in dims:
                d = cuts[j]
                cc.append([j, rest % 2 ** d, d])
                rest //= 2 ** d
            items.append({"src": src, "index": i, "n": t["n"], "kind": t["kind"], "cuts": cc, "est": est})
# pack
items.sort(key=lambda e: -e["est"])
tasks_out, cur, load = [], [], 0.0
for e in items:
    if cur and load + e["est"] > LMAX:
        tasks_out.append(cur)
        cur, load = [], 0.0
    cur.append(e)
    load += e["est"]
if cur:
    tasks_out.append(cur)
json.dump(tasks_out, open("tasks_replay.json", "w"))
print(len(items), "sub boxes in", len(tasks_out), "tasks; estimated leaves", round(sum(e["est"] for e in items)))
