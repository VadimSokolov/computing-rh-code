# hcm_arb.py: certified signs for the finite HCM test of Section ch:weak ("The HCM form of the Thorin condition").
# For alpha in {1/2, 1}, u in {1, 25, 400} and w in {2.05, 3, 5, 10, 20, 50} it encloses the 216 normalised derivatives
#   d_k = (-1)^k H_u^{(k)}(w) / H_u(w),  k = 1..6,  H_u(w) = phi(u v) phi(u / v),  phi(x) = xi(alpha) / xi(alpha + sqrt x),
# with v = (w + sqrt(w^2 - 4)) / 2, in Arb ball arithmetic (python-flint).  The Taylor coefficients of H_u in w are
# power series with ball coefficients (arb_series, 320 bits), computed from xi(s) = s (s - 1) / 2 pi^(-s/2) Gamma(s/2) zeta(s)
# with Arb's rigorous series for Gamma and zeta, so the enclosures do not use the Polya route of hcm.py.
# It reads hcm.json, the values of hcm.py by numerical differentiation, compares them with the enclosures, prints a
# summary and writes hcm_arb.json: the summary and, for each (alpha, u, w), the enclosures of H_u(w) and of d_1..d_6
# and whether each d_k is certified positive.
# Run on Hopper (one core, under a second), from the folder data/, which receives hcm_arb.json:
#   cd data && bash ../misc/tools/hopper_run.sh -c 1 -m 2G -t 30 -g hcm_arb.json ../code/hcm_arb.py hcm.json
import json, math, time
from flint import arb, arb_series, ctx

ctx.prec = 320
N = 7                      # series length: coefficients of eps^0 .. eps^6
LOGPI = arb.pi().log()
t0 = time.time()

def xi_series(S):
    return S * (S - 1) / 2 * (-(S / 2) * LOGPI).exp() * (S / 2).gamma() * arb_series.zeta(S)

def xi_point(a):
    if a == 1:
        return arb(1) / 2
    return a * (a - 1) / 2 * (-(a / 2) * LOGPI).exp() * (a / 2).gamma() * a.zeta()

book = json.load(open("hcm.json"))
rows, allvals = [], []
worst_rel_rad, worst_rel_diff = 0.0, 0.0
for astr in ["0.5", "1"]:
    a = arb(astr)
    xa = xi_point(a)
    for u in [1, 25, 400]:
        for wstr in ["2.05", "3", "5", "10", "20", "50"]:
            W = arb_series([arb(wstr), 1], prec=N)          # w + eps
            V = (W + (W * W - 4).sqrt()) / 2
            S1 = a + (u * V).sqrt()
            S2 = a + (u * V.inv()).sqrt()
            H = xa * xa * (xi_series(S1) * xi_series(S2)).inv()
            c = H.coeffs()
            ds = [(-1) ** k * math.factorial(k) * c[k] / c[0] for k in range(1, 7)]
            brow = [r for r in book[str(float(astr))] if r["u"] == u and abs(r["w"] - float(wstr)) < 1e-12][0]
            for k, (dk, bd) in enumerate(zip(ds, brow["d"]), 1):
                worst_rel_rad = max(worst_rel_rad, float(dk.rad()) / abs(float(dk.mid())))
                worst_rel_diff = max(worst_rel_diff, abs(float(dk.mid()) - bd) / abs(float(dk.mid())))
                allvals.append((float(dk.mid()), astr, u, wstr, k, dk, ds[0]))
            rows.append(dict(alpha=astr, u=u, w=wstr, H=c[0].str(25), d=[dk.str(25) for dk in ds],
                             certified_positive=[bool(dk > 0) for dk in ds]))
npos = sum(sum(r["certified_positive"]) for r in rows)
allvals.sort(key=lambda t: t[0])
smallest = []
for mid, astr, u, wstr, k, dk, d1 in allvals[:3]:
    smallest.append(dict(alpha=astr, u=u, w=wstr, k=k, d_k=dk.str(15), d_1=d1.str(15), d_k_over_d_1_to_k=(dk / d1 ** k).str(10)))
summary = dict(values=6 * len(rows), certified_positive=npos, precision_bits=ctx.prec,
               worst_relative_radius="%.2e" % worst_rel_rad,
               largest_relative_difference_from_hcm_json="%.2e" % worst_rel_diff, smallest=smallest)
print("values:", summary["values"], " certified positive:", npos)
print("worst relative radius of an enclosure:", summary["worst_relative_radius"])
print("largest relative difference between hcm.json and the enclosures:", summary["largest_relative_difference_from_hcm_json"])
for s in smallest:
    print("small value: alpha=%s u=%d w=%s k=%d d_k=%s d_1=%s d_k/d_1^k=%s" % (s["alpha"], s["u"], s["w"], s["k"], s["d_k"], s["d_1"], s["d_k_over_d_1_to_k"]))
json.dump(dict(summary=summary, rows=rows), open("hcm_arb.json", "w"), indent=1)
print("time: %.1f s" % (time.time() - t0))
