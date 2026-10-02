"""Tasks for the leaf verification run (bnb3.py).

sweep  : every piece of tasks_full.json again, with a 240 s branch and bound, so that every piece contributes kept
         leaves (the 40 with the smallest float margin found in that time and 60 random ones);
focus  : for each group (n, kind, band) of kinds m, M, l a small box around the extreme point of explore.json
         (mapped into the half torus w_2 in [0, pi] by w -> -w), and for each q group a small box around the
         centre with the largest q among the finished pieces of the full run; 900 s each.
Reads tasks_full.json, explore.json and the full run results in /scratch/.../N3-full (for q only).
"""
import json, glob, os, math
from collections import defaultdict

FULL = "/scratch/vsokolov/rh_book_repro/agents/N3-full"
tasks = json.load(open("tasks_full.json"))
EX = json.load(open("explore.json"))
PI = math.pi
out = []
for i, t in enumerate(tasks):
    s = dict(t)
    s["time_limit"] = 240
    s["label"] = f"sweep piece {t['piece']}"
    s["source_task"] = i
    out.append(s)
groups = {}
for t in tasks:
    key = (t["n"], t["kind"], tuple(t["band"]) if "band" in t else ("y", t["y"]))
    groups.setdefault(key, t)
nf = 0
for (n, kind, dom), t in sorted(groups.items(), key=lambda kv: str(kv[0])):
    if kind == "q":
        best = None
        for f in glob.glob(f"{FULL}/result_*.json"):
            r = json.load(open(f))
            rt = r["task"]
            if rt["n"] == n and rt["kind"] == "q" and r.get("best_pt") is not None:
                if best is None or r["best"] > best[0]:
                    best = (r["best"], r["best_pt"])
        if best is None:
            continue
        c = best[1]
        hs, hw = 0.2 / 32, PI / 32
        lo = [max(0.3, c[0] - hs), max(0.0, c[1] - hw)] + [v - hw for v in c[2:]]
        hi = [min(0.5, c[0] + hs), min(PI, c[1] + hw)] + [v + hw for v in c[2:]]
        label = f"focus q best {best[0]:.9g}"
    else:
        e = [o for o in EX if o["n"] == n and o["kind"] == kind and tuple(o["band"]) == dom][0]
        w = [v % (2 * PI) for v in e["w"]]
        if w[0] > PI:
            w = [(-v) % (2 * PI) for v in w]
        a, b = dom
        hx, hw = (b - a) / 32, PI / 32
        x = min(max(e["x"], a), b)
        lo = [max(a, x - hx), max(0.0, w[0] - hw)] + [v - hw for v in w[1:]]
        hi = [min(b, x + hx), min(PI, w[0] + hw)] + [v + hw for v in w[1:]]
        label = f"focus explore value {e['value']:.12g}"
    s = dict(t)
    s.pop("piece", None)
    s["lo"], s["hi"] = lo, hi
    s["time_limit"] = 900
    s["label"] = label
    out.append(s)
    nf += 1
json.dump(out, open("tasks_verify.json", "w"))
print("verify tasks", len(out), "focus", nf)
