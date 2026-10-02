# New diagnostic script (item arith_exact): splits the change from data/arith.json to the patched arith.json into
# the three input roundings of code/arith.py, at 60 digits, for each t of Table tab:ch14:W:
#   weights: P with mp.mpf(float(np.log(p))) minus P with mp.log(p); Pi + A - P moves by minus this amount;
#   zeros:   the zero sum with mp.mpf(float(x)) of the float64 lists g_*.npy minus the sum with the 70 digit strings;
#   t:       the zero sum at mp.mpf(float(t)), the binary double nearest t, minus the sum at t itself (both sides of
#            the explicit formula move together, so this shifts the two columns alike).
# Usage: python3 decompose.py   (reads zeros_1700.txt and g_*.npy; writes decompose.json)
import json, numpy as np, mpmath as mp
mp.mp.dps = 60
g_hp = [mp.mpf(s) for s in open('zeros_1700.txt').read().split()]
g_f = [mp.mpf(float(x)) for x in np.concatenate([np.load(f) for f in ['g_1_400.npy', 'g_401_1000.npy', 'g_1001_1700.npy']])]
N = 200000
isp = np.ones(N + 1, bool); isp[:2] = False
for p in range(2, int(N**0.5) + 1):
    if isp[p]: isp[p*p::p] = False
pp = []
for p in np.nonzero(isp)[0]:
    p = int(p); lp = mp.log(p); lpf = mp.mpf(float(np.log(p))); q = p; k = 1
    while q <= N: pp.append((q, lp, lpf, k * lp)); q *= p; k += 1
out = []
for ts, tf in [('0.002', 0.002), ('0.005', 0.005), ('0.01', 0.01), ('0.02', 0.02), ('0.05', 0.05), ('0.1', 0.1), ('0.15', 0.15), ('0.2', 0.2), ('0.3', 0.3)]:
    t = mp.mpf(ts); td = mp.mpf(tf)
    c = 1/(2*mp.sqrt(mp.pi*t))
    P = c*mp.fsum(lp/mp.sqrt(n)*mp.exp(-ln**2/(4*t)) for n, lp, lpf, ln in pp if ln**2/(4*t) < 300)
    Pf = c*mp.fsum(lpf/mp.sqrt(n)*mp.exp(-ln**2/(4*t)) for n, lp, lpf, ln in pp if ln**2/(4*t) < 300)
    W = mp.fsum(mp.exp(-x**2*t) for x in g_hp)
    Wf = mp.fsum(mp.exp(-x**2*t) for x in g_f)
    Wd = mp.fsum(mp.exp(-x**2*td) for x in g_hp)
    r = dict(t=ts, W=mp.nstr(W, 20), weights_dP=mp.nstr(Pf - P, 5), weights_rel_to_W=mp.nstr((P - Pf)/W, 5),
             zeros_dW=mp.nstr(Wf - W, 5), zeros_rel=mp.nstr((Wf - W)/W, 5), t_double_minus_t=mp.nstr(td - t, 5),
             t_dW=mp.nstr(Wd - W, 5), t_rel=mp.nstr((Wd - W)/W, 5))
    out.append(r)
    print(r, flush=True)
json.dump(out, open('decompose.json', 'w'), indent=1)
