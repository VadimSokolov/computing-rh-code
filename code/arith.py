# code/arith.py: the polar, archimedean and prime terms of the explicit formula eq:ch2:Warith for the heat trace W(t),
# and W(t) as the sum over the first 1700 zeros, at nine values of t; writes arith.json, the data of Table tab:ch14:W
# (ch/ch14.tex) and of fig/heattrace.pdf (code/plots3.py, section 2). This version (September 2026, from
# misc/repro/reconstruct/p13_arith_exact) replaces the first archived script and writes the same format. Every change
# removes a loss of precision:
#  1. the weights Lambda(p^k) = log p are mp.log(p) at working precision; the first version took np.log(p) in double
#     precision and then mp.mpf(float(L[i])), an error of up to 9e-16 in each weight (half an ulp of log p < 12.3),
#     which put an absolute error of 2e-19 to 2e-18 into P, hence into Pi + A - P, for t >= 0.05: relative 4e-14 at
#     t = 0.05, 4e-2 at t = 0.2, and the wrong sign at t = 0.3, where W = 9.3e-27;
#  2. the zero ordinates are mp.mpf of the 70 digit strings in zeros_1700.txt (written by zeros_hp.py with Arb), where
#     the first version took mp.mpf(float(x)) of the float64 lists g_*.npy (errors up to half an ulp: 2.3e-13 near
#     height 2100, 8.4e-16 for gamma_1), which moved the zero sum by up to 7e-15 relative, in its fifteenth digit
#     for t >= 0.1;
#  3. t is mp.mpf of a decimal string, not of the binary double nearest to it (a relative shift of W up to 2e-15);
#  4. the working precision is 60 digits by default (first version: 40). At t = 0.3 the explicit formula cancels
#     26.4 digits: with changes 1 to 3 at 40 digits, Pi + A - P = 9.32464654810623e-27 against 9.32464654810622e-27
#     from the zeros (14.9 digits of agreement); at 60 and at 80 digits the two agree to 26.8 digits, the size of
#     the truncation of the prime sum at 2*10^5 (1.4e-53);
#  5. the quadrature for A is that of the first version (same breakpoints, default degree), but the script prints the
#     error estimate of mpmath and would raise the degree if the estimate exceeded 10^(5 - dps); at 40, 60 and 80
#     digits the default degree sufficed.
# The prime powers n <= 2*10^5, the cutoff (log n)^2/4t < 300 and the output format are those of the first version.
# Usage: python3 arith.py [dps] [output]   (defaults 60 and arith.json); diagnostics go to <output>.diag.json.
import sys, json, time, numpy as np, mpmath as mp

dps = int(sys.argv[1]) if len(sys.argv) > 1 else 60
out = sys.argv[2] if len(sys.argv) > 2 else 'arith.json'
mp.mp.dps = dps
t0 = time.time()
g = [mp.mpf(s) for s in open('zeros_1700.txt').read().split()]
assert len(g) == 1700
N = 200000
isp = np.ones(N + 1, bool); isp[:2] = False
for p in range(2, int(N**0.5) + 1):
    if isp[p]: isp[p*p::p] = False
pp = []                                    # (n, log p, log n) for the prime powers n = p^k <= N
for p in np.nonzero(isp)[0]:
    p = int(p); lp = mp.log(p); q = p; k = 1
    while q <= N: pp.append((q, lp, k * lp)); q *= p; k += 1
pp.sort()
print('dps %d, %d prime powers up to %d, %d zeros up to %s' % (dps, len(pp), N, len(g), mp.nstr(g[-1], 12)), flush=True)

def m(r): return mp.re(mp.digamma(mp.mpf(1)/4 + 1j*r/2)) - mp.log(mp.pi)

def quad(f, pts):
    tol = mp.mpf(10)**(5 - dps)
    deg = None
    while True:
        v, e = mp.quad(f, pts, error=True, maxdegree=deg)
        if e < tol or (deg is not None and deg >= 12): return v, e, deg
        deg = 8 if deg is None else deg + 2

rows = []; diag = []
for ts in ['0.002', '0.005', '0.01', '0.02', '0.05', '0.1', '0.15', '0.2', '0.3']:
    t = mp.mpf(ts)
    Pi = mp.exp(t/4)
    I, qerr, qdeg = quad(lambda r: mp.exp(-r*r*t)*m(r), [0, 5, 20, 80, mp.inf])
    A = I/(2*mp.pi)
    P = mp.fsum(lp/mp.sqrt(n)*mp.exp(-ln**2/(4*t)) for n, lp, ln in pp if ln**2/(4*t) < 300)/(2*mp.sqrt(mp.pi*t))
    Wa = Pi + A - P
    Wz = mp.fsum(mp.exp(-x**2*t) for x in g)
    dig = mp.log10((Pi + abs(A) + P)/abs(Wz))
    rows.append([float(t), mp.nstr(Pi, 12), mp.nstr(A, 12), mp.nstr(P, 12), mp.nstr(Wa, 15), mp.nstr(Wz, 15), float(dig)])
    agree = float(-mp.log10(abs(Wa - Wz)/abs(Wz))) if Wa != Wz else float(dps)
    diag.append(dict(t=ts, Pi=mp.nstr(Pi, 30), A=mp.nstr(A, 30), P=mp.nstr(P, 30), Wa=mp.nstr(Wa, 30), Wz=mp.nstr(Wz, 30),
                     Wa_minus_Wz=mp.nstr(Wa - Wz, 5), digits_of_agreement=agree, quad_error_estimate=mp.nstr(qerr/(2*mp.pi), 5),
                     quad_maxdegree=qdeg, digits_lost=mp.nstr(dig, 20)))
    print(rows[-1], flush=True)
    print('   Wa = %s\n   Wz = %s\n   Wa - Wz = %s (agreement %.1f digits), quadrature error estimate %s (maxdegree %s), %.1f s'
          % (mp.nstr(Wa, 30), mp.nstr(Wz, 30), mp.nstr(Wa - Wz, 5), agree, mp.nstr(qerr/(2*mp.pi), 5), qdeg, time.time() - t0), flush=True)
json.dump(rows, open(out, 'w'))
json.dump(diag, open(out.replace('.json', '') + '.diag.json', 'w'), indent=1)
print('done in %.1f s' % (time.time() - t0))
