"""Status of every group (n, kind, band or y) and certified_bounds.json (input of arb_thresholds.py).

Kind m comes from the full run (bnb2.py, /scratch/.../N3-full). Kinds M and l come from the run with the certified
lower bound m_low of |P_n| (bnb4.py, /scratch/.../N3-mlow), which supersedes the M and l pieces of the full run.
Every final box of m, M and l is proved again in Arb by the replay (bnb7.py, aggregate_replay.py). The second
condition comes from the run with G (bnb8.py, /scratch/.../N3-g), in which every final box is proved in Arb; it
supersedes the searches for q in logarithmic form (N3-full, N3-q2, N3-q3), which were slow and are not used.
A group counts as certified when every one of its pieces returned "certified" (m, M, l) or "proved" (G)."""
import json, glob, os
from collections import defaultdict

D = "/scratch/vsokolov/rh_book_repro/agents"


def load(run, tasks_name):
    tasks = json.load(open(f"{D}/N3-{run}/{tasks_name}"))
    res = {}
    for f in glob.glob(f"{D}/N3-{run}/result_*.json"):
        res[int(os.path.basename(f)[7:-5])] = json.load(open(f))
    return tasks, res


groups = defaultdict(list)
full_tasks, full_res = load("full", "tasks_full.json")
for i, t in enumerate(full_tasks):
    if t["kind"] == "m":
        groups[(t["n"], "m", tuple(t["band"]), t["target"])].append((t, full_res.get(i)))
mlow_tasks, mlow_res = load("mlow", "tasks_mlow.json")
for i, t in enumerate(mlow_tasks):
    groups[(t["n"], t["kind"], tuple(t["band"]), t["target"])].append((t, mlow_res.get(i)))
g_tasks, g_res = load("g", "tasks_g.json")
for i, e in enumerate(g_tasks):
    t = full_tasks[e["q_index"]]
    r = g_res.get(i)
    r = r[0] if r else None
    groups[(t["n"], "G", ("y", t["y"]), t["target"])].append((t, r))
summary = []
cb = defaultdict(dict)
for g in sorted(groups, key=lambda k: (k[0], k[1], str(k[2]))):
    n, kind, dom, target = g
    items = groups[g]
    rs = [r for _, r in items if r]
    if kind == "G":
        st = [r["status"] if r else "missing" for _, r in items]
        ok = all(s == "proved" for s in st)
        row = {"n": n, "kind": "G (second condition)", "y": dom[1], "q_target": target, "pieces": len(st),
               "proved": sum(s == "proved" for s in st), "all_proved": ok, "other": sorted(set(s for s in st if s != "proved")),
               "boxes": sum(r.get("processed") or 0 for r in rs), "leaves": sum(r.get("leaves", 0) for r in rs),
               "arb_failed": sum(r.get("arb_failed", 0) for r in rs), "arb_rescued": sum(r.get("arb_rescued", 0) for r in rs),
               "float_seconds": round(sum(r.get("float_seconds", 0) for r in rs), 1),
               "arb_seconds": round(sum(r.get("arb_seconds", 0) for r in rs), 1),
               "min_centre_G": min((r["min_centre_G"] for r in rs if r.get("min_centre_G") is not None), default=None)}
        summary.append(row)
        if ok:
            cb[(n, (1.2, 1.3))]["q_y"] = dom[1]
            cb[(n, (1.2, 1.3))]["q_target"] = target
        continue
    st = [r["status"] if r else "missing" for _, r in items]
    ok = all(s == "certified" for s in st)
    bests = [r["best"] for r in rs if r.get("best") is not None]
    ext = (min(bests) if kind == "m" else max(bests)) if bests else None
    row = {"n": n, "kind": kind, "band": list(dom), "target": target, "pieces": len(st),
           "certified": sum(s == "certified" for s in st), "all_certified": ok,
           "other": sorted(set(s for s in st if s != "certified")),
           "boxes": sum(r.get("processed", 0) for r in rs), "leaves": sum(r.get("leaves", 0) for r in rs),
           "cpu_seconds": round(sum(r.get("seconds", 0) for r in rs), 1), "extreme_centre_value": ext}
    ac = [r["arb_check"] for r in rs if "arb_check" in r and "leaves" in r["arb_check"]]
    if ac:
        row["sample_arb_leaves"] = sum(a["leaves"] for a in ac)
        row["sample_arb_mirror_ok"] = sum(a["mirror_ok"] for a in ac)
        row["sample_arb_naive_ok"] = sum(a["naive_ok"] for a in ac)
        row["sample_arb_max_bound_diff"] = max(a["max_bound_diff"] for a in ac)
    if kind in ("M", "l"):
        row["m_low"] = items[0][0]["m_low"]
    summary.append(row)
    if ok:
        cb[(n, dom)][kind] = repr(float(target))
bounds = []
for (n, dom), v in sorted(cb.items()):
    e = {"n": n, "band": list(dom)}
    e.update(v)
    bounds.append(e)
json.dump(bounds, open("certified_bounds.json", "w"), indent=1)
json.dump(summary, open("full_status.json", "w"), indent=1)
for s in summary:
    print(json.dumps(s))
