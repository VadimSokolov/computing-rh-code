"""Thresholds of Theorem rz:thm:torus (first condition) in Arb from directed bounds on the torus extremes.

Input certified_bounds.json: list of {n, band: [a,b], m, M, l} (m a lower bound, M and l upper bounds, as decimal
strings). r(a) and r'(a) are exact finite sums, enclosed in Arb. For each entry we give the real root y* of
  g(y) = log(y/pi) - 2/y - M - (l r + r' + r/y) / (y (m - r/y)),
which is increasing on y > r/m, enclosed by Arb bisection, and the smallest integer Y with g(Y) > 0 certified.
Also evaluated: the same with the two decimal directed values (m down, M and l up), and with the printed m, M of
Table rz:tab:torus (and our l). For the entries with band [1.2,1.3] the second condition is checked separately.
"""
import json, math
from fractions import Fraction as Fr
from flint import arb, fmpq, ctx
from torus_common import exact_coeffs

ctx.prec = 128
PI = arb.pi()
PRINTED = {16: ("1.06", "0.86", 30), 17: ("1.00", "0.89", 34), 18: ("0.95", "0.93", 37), 19: ("0.90", "0.96", 42),
           20: ("0.85", "0.98", 46), 21: ("0.80", "1.01", 50), 22: ("0.76", "1.04", 56), 23: ("0.72", "1.06", 62),
           24: ("0.68", "1.09", 69), 32: ("0.36", "1.26", 214), 36: ("0.20", "1.33", 566), 40: ("0.056", "1.40", 4898)}
SINGLE = {7: 14, 9: 23, 10: 31, 11: 42, 12: 57, 13: 85, 14: 131, 15: 229}
HEIGHT = {7: 15, 9: 24, 10: 32, 11: 43, 12: 58, 13: 86, 14: 132, 15: 230, 16: 500, 17: 35, 18: 250, 19: 43,
          20: 47, 21: 51, 22: 57, 23: 62, 24: 110, 32: 215, 36: 570}


def rr(n, a):
    ab = exact_coeffs(n)
    fa = Fr(repr(float(a)))
    xa = arb(fmpq(fa.numerator, fa.denominator))  # exact 6/5 or 1
    r = arb(0)
    rp = arb(0)
    for k, (al, be) in enumerate(ab, 1):
        w = (-2 * xa * arb(k).log()).exp()
        A = abs(arb(fmpq(al.numerator, al.denominator)))
        r += A * w
        rp += 2 * arb(k).log() * A * w
    return r, rp


def up3(b):
    """Decimal string of an arb rounded up at the third decimal (for printing r and r' in a table)."""
    v = b.upper()
    k = int(float((v * 1000).floor().mid()))
    while arb(k) / 1000 < v:
        k += 1
    return f"{k // 1000}.{k % 1000:03d}"


def g(y, m, M, l, r, rp):
    y = arb(y)
    d = m - r / y
    if not (d > 0):
        return None
    return (y / PI).log() - 2 / y - M - (l * r + rp + r / y) / (y * d)


def solve(m, M, l, r, rp):
    m, M, l = arb(m), arb(M), arb(l)
    lo = float((r / m).upper()) * (1 + 1e-9)
    hi = max(2 * lo, 10.0)
    while True:
        v = g(hi, m, M, l, r, rp)
        if v is not None and v > 0:
            break
        hi *= 2
    # bisection on floats with Arb sign decisions
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        v = g(mid, m, M, l, r, rp)
        if v is not None and v > 0:
            hi = mid
        else:
            lo = mid
        if hi - lo < 1e-7:
            break
    Y = math.ceil(hi)
    # certify g(Y) > 0 and report g(Y-1)
    gY = g(Y, m, M, l, r, rp)
    gY1 = g(Y - 1, m, M, l, r, rp)
    return {"root_in": [lo, hi], "Y": Y, "g(Y)": gY.str(8), "g(Y)>0": bool(gY > 0),
            "g(Y-1)": None if gY1 is None else gY1.str(8)}


def down2(s):
    return str(math.floor(float(s) * 100 + 1e-12) / 100)


def up2(s):
    return str(math.ceil(float(s) * 100 - 1e-12) / 100)


ATT = {}
for o in json.load(open("arb_points.json")):
    key = (o["n"], tuple(o["band"]))
    if o["kind"] == "m":
        ATT.setdefault(key, {})["m"] = arb(o["optimiser_point_|P|"]).upper()
    elif o["kind"] == "M":
        ATT.setdefault(key, {})["M"] = arb(o["optimiser_point_-Re(P'/P)"]).lower()
    else:
        ATT.setdefault(key, {})["l"] = arb(o["optimiser_point_|P'/P|"]).lower()


def lower_threshold(att, r, rp):
    """Values attained at points (m above the true minimum, M and l below the true maxima) make the first
    condition weaker than with the true extremes, so its root is a lower bound for the true threshold root."""
    s = solve(att["m"], att["M"], att["l"], r, rp)
    return {"root_in": s["root_in"], "smallest_integer_possible": s["Y"], "g(Y-1) with attained values": s["g(Y-1)"]}


out = []
for e in json.load(open("certified_bounds.json")):
    n, (a, b) = e["n"], e["band"]
    if not all(k in e for k in ("m", "M", "l")):
        print("incomplete", e)
        continue
    r, rp = rr(n, a)
    rec = {"n": n, "band": [a, b], "r(a)": r.str(15), "r'(a)": rp.str(15), "bounds": {k: e[k] for k in ("m", "M", "l")}}
    rec["certified_3dp"] = solve(e["m"], e["M"], e["l"], r, rp)
    if (n, (a, b)) in ATT and len(ATT[(n, (a, b))]) == 3:
        rec["attained"] = {k: v.str(15) for k, v in ATT[(n, (a, b))].items()}
        rec["threshold_lower_bound"] = lower_threshold(ATT[(n, (a, b))], r, rp)
    r3, rp3 = up3(r), up3(rp)
    rec["table_r_up3"] = [r3, rp3]
    rec["table_check"] = solve(e["m"], e["M"], e["l"], arb(r3), arb(rp3))
    rec["directed_2dp"] = {"m": down2(e["m"]), "M": up2(e["M"]), "l": up2(e["l"])}
    rec["certified_2dp"] = solve(rec["directed_2dp"]["m"], rec["directed_2dp"]["M"], rec["directed_2dp"]["l"], r, rp)
    if a == 1.2 and n in PRINTED:
        pm, pM, pth = PRINTED[n]
        rec["printed"] = {"m": pm, "M": pM, "threshold": pth}
        rec["printed_mM_our_l"] = solve(pm, pM, e["l"], r, rp)
    if a == 1.0 and n in SINGLE:
        rec["printed"] = {"threshold": SINGLE[n]}
    if n in HEIGHT:
        rec["count_height"] = HEIGHT[n]
    if "q_y" in e:
        Y = rec["certified_3dp"]["Y"]
        rec["second_condition"] = {"certified_at_y": e["q_y"], "q_target": e["q_target"], "covers_threshold": e["q_y"] <= Y}
    out.append(rec)
    print(json.dumps(rec), flush=True)
json.dump(out, open("arb_thresholds.json", "w"), indent=1)
