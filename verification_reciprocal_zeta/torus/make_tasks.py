"""Generate the branch and bound pieces (tasks_full.json) from the floating point extremes of explore.json."""
import json, math
import numpy as np
from torus_common import primes_upto

EX = json.load(open("explore.json"))
TH = {o["n"]: o for o in json.load(open("thresholds.json"))}
tasks = []


def pieces(n, a, b, splits):
    K = len(primes_upto(n))
    lo = [a, 0.0] + [0.0] * (K - 1)
    hi = [b, math.pi] + [2 * math.pi] * (K - 1)
    out = [(lo, hi)]
    for dim, s in splits:
        if dim > K:
            continue
        new = []
        for L, H in out:
            w = (H[dim] - L[dim]) / s
            for i in range(s):
                L2, H2 = list(L), list(H)
                L2[dim] = L[dim] + i * w
                H2[dim] = L[dim] + (i + 1) * w
                new.append((L2, H2))
        out = new
    return out


for r in EX:
    n, (a, b), kind, v = r["n"], r["band"], r["kind"], r["value"]
    if n == 40:
        continue
    if a == 1.0 and n >= 16:
        continue  # the single band is used only for n = 7, 9..15 (m = 0 for n >= 18)
    if kind == "m":
        target = math.floor((v - 3e-4) * 1000) / 1000
        splits = []
    elif kind == "M":
        target = math.ceil((v + 3e-4) * 1000) / 1000
        splits = [(1, 4), (2, 4), (3, 4)]
    else:
        target = math.ceil((v + 3e-4) * 1000) / 1000
        splits = [(1, 4), (2, 4), (3, 4)]
    for i, (lo, hi) in enumerate(pieces(n, a, b, splits)):
        tasks.append({"n": n, "kind": kind, "band": [a, b], "target": target, "float_value": v,
                      "lo": lo, "hi": hi, "piece": i, "time_limit": 100 * 60})
# second condition at y = floor(float y1), band [1.2,1.3] cases
for n, o in TH.items():
    if "y1_unrounded" not in o or n == 40:
        continue
    y = float(math.floor(o["y1_unrounded"]))
    target = math.log(y / math.pi) - 1 / y - 1e-12
    K = len(primes_upto(n))
    lo = [0.3, 0.0] + [0.0] * (K - 1)
    hi = [0.5, math.pi] + [2 * math.pi] * (K - 1)
    base = [(lo, hi)]
    for dim, s in [(1, 4), (2, 4)]:
        new = []
        for L, H in base:
            w = (H[dim] - L[dim]) / s
            for i in range(s):
                L2, H2 = list(L), list(H)
                L2[dim] = L[dim] + i * w
                H2[dim] = L[dim] + (i + 1) * w
                new.append((L2, H2))
        base = new
    for i, (L, H) in enumerate(base):
        tasks.append({"n": n, "kind": "q", "y": y, "target": target, "lo": L, "hi": H, "piece": i,
                      "time_limit": 100 * 60})
json.dump(tasks, open("tasks_full.json", "w"))
from collections import Counter
print(len(tasks), Counter((t["kind"]) for t in tasks))
print(sorted(set((t["n"], t["kind"], t["band"][0] if "band" in t else t["y"], t["target"]) for t in tasks)))
