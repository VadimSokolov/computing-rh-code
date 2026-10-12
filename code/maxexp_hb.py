# maxexp_hb.py: the Hermite Biehler test with Hecke's cut and with a moving cut, for Section sec:ch16:max (ch/ch16.tex,
# Section 8.5.4) and Figure fig:maxexp:dominance. Added for the book in October 2026, from the first author's note on maxima of
# exponentials (misc/reviews/nick-maxexp/, Section 6.4 and Figure 3 there), recomputed independently.
# For c > 0, xi(s, Delta) = A_c(s) + A_{1/c}(1 - s) with A_c(s) = int_c^oo Delta(iy) y^{s+9/2} dy, and with s = 1/2 + w, y = e^v,
#   A_c(1/2 + w) = int_{log c}^oo K(v) e^{w v} dv,   A_{1/c}(1/2 - w) = int_{-log c}^oo K(v) e^{-w v} dv,
# K(v) = e^{6v} Delta(i e^v) = K(-v) the kernel of Section sec:ch16:dbn. Gauss Legendre with 40 nodes on panels of width at most
# 1/4 up to v = 4.5, in double precision (the ratios below need a few digits). The script computes
#  1. with Hecke's cut c = 1: the ratio |A(1-s)|/|A(s)| at theta = 0 for sigma = 0.52 and 4.5, the zeros of xi(1/2 + i theta)
#     = 2 Re A(1/2 + i theta) below 25, and the ordinates theta <= 25 where |A(s)| = |A(1-s)| on the lines Re s = 0.51, ..., 2.5
#     and 4.5, by a scan with step 0.01 and Brent's method;
#  2. with a moving cut: for theta in [0, 25] (step 0.05), rho(theta) = max over c in [0.35, 1.2] (36 values) of the least ratio
#     |A_c(s)|/|A_{1/c}(1-s)| over 120 values of sigma in (1/2, 14.5]; rho > 1 means that one cut c makes A_c dominate on the
#     whole horizontal ray, which excludes zeros of xi(s, Delta) on it.
# Run from any folder: python3 maxexp_hb.py; writes maxexp_hb.json there.
import json, time
import numpy as np
from scipy.optimize import brentq

t0 = time.time()
out = {}
def put(k, v):
    out[k] = v
    print(k, v, round(time.time() - t0, 1), flush=True)
def K(v):
    y = np.exp(np.abs(v))
    p = np.ones_like(y)
    for n in range(1, 40):
        p *= (1 - np.exp(-2 * np.pi * n * y)) ** 24
    return np.exp(6 * np.abs(v) - 2 * np.pi * y) * p
X, WX = np.polynomial.legendre.leggauss(40)
NODES = {}
def nodes(a, b=4.5, h=0.25):
    key = round(float(a), 12)
    if key not in NODES:
        m = int(np.ceil((b - a) / h)); edges = np.linspace(a, b, m + 1)
        v = np.concatenate([(edges[i] + edges[i + 1]) / 2 + (edges[i + 1] - edges[i]) / 2 * X for i in range(m)])
        w = np.concatenate([(edges[i + 1] - edges[i]) / 2 * WX for i in range(m)])
        NODES[key] = (v, w * K(v))
    return NODES[key]
def halves(c, sig, th):
    # A_c(s) and A_{1/c}(1 - s) on the grid s = sig + i th (arrays; sig and th broadcast)
    v1, k1 = nodes(np.log(c)); v2, k2 = nodes(-np.log(c))
    w = (np.asarray(sig) - 0.5)[..., None] + 1j * np.asarray(th)[..., None]
    return (np.exp(w * v1) * k1).sum(-1), (np.exp(-w * v2) * k2).sum(-1)
A1, A2 = halves(1.0, np.array([0.5]), np.array([0.0]))
put('xi(1/2, Delta) = 2 A(1/2) (maxlaw.py: 0.00154487936039502720604)', float(2 * A1[0].real))
a1, a2 = halves(1.0, np.array([0.52, 4.5]), np.array([0.0, 0.0]))
put('|A(1-s)|/|A(s)| at theta = 0, sigma = 0.52 and 4.5', [float(abs(a2[0]) / abs(a1[0])), float(abs(a2[1]) / abs(a1[1]))])
# zeros on the critical line
f = lambda th: float(halves(1.0, np.array([0.5]), np.array([th]))[0][0].real)
ths = np.arange(0.01, 25.0, 0.01); vals = np.array([f(t) for t in ths])
zs = [brentq(f, ths[i], ths[i + 1], xtol=1e-13) for i in range(len(ths) - 1) if vals[i] * vals[i + 1] < 0]
put('zeros of xi(1/2 + i theta, Delta) below 25', [round(z, 8) for z in zs])
# the ordinates where the two halves have equal modulus
def lr(sig, th):
    a, b = halves(1.0, np.array([sig]), np.array([th]))
    return float(np.log(abs(b[0])) - np.log(abs(a[0])))
cross = {}
for sig in [0.51, 0.55, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5, 4.5]:
    a, b = halves(1.0, np.full(ths.shape, sig), ths)
    g = np.log(np.abs(b)) - np.log(np.abs(a))
    cross[str(sig)] = [round(brentq(lambda t: lr(sig, t), ths[i], ths[i + 1], xtol=1e-12), 4) for i in range(len(ths) - 1) if g[i] * g[i + 1] < 0]
put('ordinates <= 25 with |A(s)| = |A(1-s)|, by Re s', cross)
# the moving cut
TH = np.round(np.arange(0, 25.0001, 0.05), 2)
SIG = 0.5 + np.geomspace(0.002, 14, 120)
CS = np.round(np.concatenate([np.linspace(0.35, 1.2, 35), [1.0]]), 4)
best = np.full(TH.shape, -np.inf); bestc = np.zeros(TH.shape); fixed = None
for c in CS:
    v1, k1 = nodes(np.log(c)); v2, k2 = nodes(-np.log(c))
    E1 = np.exp(1j * TH[:, None] * v1[None, :]); E2 = np.exp(-1j * TH[:, None] * v2[None, :])
    m = np.full(TH.shape, np.inf)
    for sg in SIG:   # least ratio over sigma, for each theta
        a = E1 @ (k1 * np.exp((sg - 0.5) * v1)); b = E2 @ (k2 * np.exp(-(sg - 0.5) * v2))
        m = np.minimum(m, np.abs(a) / np.abs(b))
    if c == 1.0: fixed = m
    upd = m > best; best[upd] = m[upd]; bestc[upd] = c
def at(t):
    i = int(np.argmin(np.abs(TH - t)))
    return {'theta': float(TH[i]), 'rho': round(float(best[i]), 5), 'best c': float(bestc[i]), 'least ratio with c = 1': round(float(fixed[i]), 5)}
put('moving cut at theta = 5, 10.5, 12', [at(5), at(10.5), at(12)])
runs = []; i = 0
while i < len(TH):
    if best[i] < 1:
        j = i
        while j + 1 < len(TH) and best[j + 1] < 1: j += 1
        runs.append([float(TH[i]), float(TH[j]), round(float(best[i:j + 1].min()), 5)])
        i = j + 1
    else:
        i += 1
put('theta intervals (grid) where no cut in [0.35, 1.2] dominates, with the least rho', runs)
zz = [0.0] + zs + [25.0]
put('largest rho between consecutive zeros, and where', [{'between': [round(zz[k], 3), round(zz[k + 1], 3)], 'max rho': round(float(best[(TH > zz[k]) & (TH < zz[k + 1])].max()), 5),
     'at theta': float(TH[(TH > zz[k]) & (TH < zz[k + 1])][np.argmax(best[(TH > zz[k]) & (TH < zz[k + 1])])])} for k in range(len(zz) - 1)])
out['rho'] = {'theta': TH.tolist(), 'rho': [round(float(x), 5) for x in best], 'best c': bestc.tolist(), 'least ratio with c = 1': [round(float(x), 5) for x in fixed]}
out['seconds'] = time.time() - t0
json.dump(out, open('maxexp_hb.json', 'w'), indent=1)
print('done', round(time.time() - t0, 1))
