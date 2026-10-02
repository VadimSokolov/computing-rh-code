"""Arb values of |P_n|, -Re(P_n'/P_n) and |P_n'/P_n| at explicit points of the band and torus.

A value at a point is a certified bound on the extreme on the other side: |P_n(pt)| >= m_n, -Re(P'/P)(pt) <= M_n,
|P'/P|(pt) <= l_n. Points: the optimiser points of explore.json (their float coordinates, read exactly), and the
clean points x = 1.2 (or 1 for the single band) with every phase 0 (for M) or every phase pi (the Liouville point,
for m and l). Also the printed values of Table rz:tab:torus are compared.
"""
import json
from fractions import Fraction as Fr
from flint import arb, acb, fmpq, ctx
from torus_common import exact_coeffs, primes_upto, vp

ctx.prec = 200
PRINTED = {16: (1.06, 0.86), 17: (1.00, 0.89), 18: (0.95, 0.93), 19: (0.90, 0.96), 20: (0.85, 0.98),
           21: (0.80, 1.01), 22: (0.76, 1.04), 23: (0.72, 1.06), 24: (0.68, 1.09), 32: (0.36, 1.26),
           36: (0.20, 1.33), 40: (0.056, 1.40)}
PI = arb.pi()


def A(q):
    return arb(fmpq(q.numerator, q.denominator))


def evalP(n, x, w):
    """x, w: arb values (w list over primes up to n). Returns P, -P' (acb)."""
    ab = exact_coeffs(n)
    ps = primes_upto(n)
    P = acb(0)
    N = acb(0)
    for k, (a, b) in enumerate(ab, 1):
        phi = arb(0)
        for p, wp in zip(ps, w):
            e = vp(k, p)
            if e:
                phi += e * wp
        lk = arb(k).log()
        t = A(b) * (-2 * x * lk).exp() * acb(0, -phi).exp()
        P += t
        N += 2 * lk * t
    return P, N


EX = json.load(open("explore.json"))
out = []
for r in EX:
    n, (a, b), kind = r["n"], r["band"], r["kind"]
    x = arb(r["x"])
    w = [arb(t) for t in r["w"]]
    P, N = evalP(n, x, w)
    R = N / P
    rec = {"n": n, "band": [a, b], "kind": kind, "float_value": r["value"]}
    if kind == "m":
        rec["optimiser_point_|P|"] = abs(P).str(15)
        # Liouville point: x at the lower edge, every phase pi
        xl = arb(6) / 5 if a == 1.2 else arb(1)
        P2, N2 = evalP(n, xl, [PI] * len(w))
        rec["liouville_point_|P|"] = abs(P2).str(15)
        rec["liouville_point_|P'/P|"] = abs(N2 / P2).str(15)
        up = min(float(abs(P).upper()), float(abs(P2).upper()))
        rec["certified_upper_bound_m"] = up
        if a == 1.2 and n in PRINTED:
            rec["printed_m"] = PRINTED[n][0]
            rec["printed_m_exceeds_true_min"] = bool(min(abs(P), abs(P2)) < arb(str(PRINTED[n][0])))
    elif kind == "M":
        rec["optimiser_point_-Re(P'/P)"] = R.real.str(15)
        xl = arb(6) / 5 if a == 1.2 else arb(1)
        P0, N0 = evalP(n, xl, [arb(0)] * len(w))
        rec["zero_phase_point_-Re(P'/P)"] = (N0 / P0).real.str(20)
        if a == 1.2 and n in PRINTED:
            rec["printed_M"] = PRINTED[n][1]
            rec["true_max_exceeds_printed_M"] = bool((N0 / P0).real > arb(str(PRINTED[n][1])))
    else:
        rec["optimiser_point_|P'/P|"] = abs(R).str(15)
    out.append(rec)
    print(json.dumps(rec), flush=True)
json.dump(out, open("arb_points.json", "w"), indent=1)
