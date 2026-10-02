"""Aggregate the leaf verification run (bnb3.py results in /scratch/.../N3-verify) per group (n, kind, band or y)."""
import json, glob
from collections import defaultdict

D = "/scratch/vsokolov/rh_book_repro/agents/N3-verify"
G = defaultdict(lambda: {"tasks": 0, "focus_tasks": 0, "statuses": defaultdict(int), "leaves": 0, "small": 0,
                         "mirror_fail": 0, "naive_fail": 0, "point_viol": 0, "min_float_margin": None,
                         "max_bound_diff": 0.0, "max_rel_bound_diff": 0.0, "naive_boxes_max": 0, "no_leaves": 0})
fails = []
for f in glob.glob(f"{D}/result_*.json"):
    r = json.load(open(f))
    t = r["task"]
    key = (t["n"], t["kind"], tuple(t["band"]) if "band" in t else ("y", t["y"]))
    g = G[key]
    g["tasks"] += 1
    if t.get("label", "").startswith("focus"):
        g["focus_tasks"] += 1
    g["statuses"][r["status"]] += 1
    rows = r.get("leaf_rows", [])
    if not rows:
        g["no_leaves"] += 1
    for row in rows:
        g["leaves"] += 1
        if row["cls"] == "small":
            g["small"] += 1
            fm = row["float_margin"]
            if fm is not None and (g["min_float_margin"] is None or fm < g["min_float_margin"]):
                g["min_float_margin"] = fm
        if not row["mirror_ok"]:
            g["mirror_fail"] += 1
            fails.append({"file": f, "kind": "mirror", "row": row})
        if not row["naive_ok"]:
            g["naive_fail"] += 1
            fails.append({"file": f, "kind": "naive", "row": {k: row[k] for k in ("cls", "float_margin", "float_bound", "mirror_bound", "naive_boxes", "h")}})
        if row["point_violation"]:
            g["point_viol"] += 1
            fails.append({"file": f, "kind": "point", "row": row})
        if row.get("bound_diff") is not None:
            g["max_bound_diff"] = max(g["max_bound_diff"], row["bound_diff"])
            den = max(1.0, abs(row["mirror_bound"]))
            g["max_rel_bound_diff"] = max(g["max_rel_bound_diff"], row["bound_diff"] / den)
        g["naive_boxes_max"] = max(g["naive_boxes_max"], row["naive_boxes"])
out = []
tot = defaultdict(float)
for key in sorted(G, key=lambda k: (k[0], k[1], str(k[2]))):
    g = G[key]
    g["statuses"] = dict(g["statuses"])
    row = {"n": key[0], "kind": key[1], "domain": key[2]}
    row.update(g)
    out.append(row)
    for k in ("tasks", "leaves", "small", "mirror_fail", "naive_fail", "point_viol"):
        tot[k] += g[k]
    tot["max_bound_diff"] = max(tot["max_bound_diff"], g["max_bound_diff"])
    print(json.dumps(row))
print("TOTAL", json.dumps(tot))
json.dump({"groups": out, "total": tot, "failures": fails[:200]}, open("verify_aggregate.json", "w"), indent=1)
