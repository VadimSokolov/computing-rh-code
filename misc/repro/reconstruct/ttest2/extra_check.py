# Reconstructed check, not the authors' script. Two side questions about the archived code/ttest.py and
# code/ttest2.py (Table tab:tilt and Section sec:ch4:balls of ch/ch04.tex), using the archived code/tilted.py.
# (1) What ttest2.py would print if it converted the tilted value with 60 digits instead of 45 (tomp uses
#     str(45)), against the same 50 digit reference: the errors the table would show without the 45 digit floor.
# (2) Why the archived ttest.py prints an error of 2.0e-15 at theta = 14: it builds s from the double 0.05
#     (acb(0.05, th)) and the reference from the double 0.55 (mp.mpf(0.55)), two points 4.2e-17 apart.
#     We compare the change of xi'/xi between these two points with the error ttest.py prints, and recompute
#     ttest.py's Horner value at theta = 14 with the exact inputs of ttest2.py.
import json
import mpmath as mp
from flint import acb, arb
import tilted

def dl(s): return 1/s+1/(s-1)-mp.log(mp.pi)/2+mp.digamma(s/2)/2+mp.zeta(s,derivative=1)/mp.zeta(s)
def tomp(z, n): return mp.mpc(mp.mpf(z.real.mid().str(n,radius=False)), mp.mpf(z.imag.mid().str(n,radius=False)))
out = dict(conv60=[], ttest14={})
for delta, prec, h, xmax, ths in [(0.10, 200, 0.004, 4.0, [200, 500, 1000]), (0.02, 220, 0.001, 4.2, [1000, 2000, 5000])]:
    S = tilted.setup(delta, prec, h, xmax)
    for th in ths:
        A, B = tilted.xi_pair_direct(S, acb(arb('0.05'), th)); r = B/A
        mp.mp.dps = 50; ref50 = dl(mp.mpf('0.55') + 1j*mp.mpf(th))
        mp.mp.dps = 70; e60 = abs(tomp(r, 60) - ref50); e60re = abs(tomp(r, 60).real - ref50.real)
        row = dict(delta=delta, theta=th, err_conv60_ref50=mp.nstr(e60, 4), err_re_conv60_ref50=mp.nstr(e60re, 4))
        out['conv60'].append(row); print(json.dumps(row), flush=True)
mp.mp.dps = 50
sig_t = mp.mpf('0.5') + mp.mpf(0.05)          # where ttest.py evaluates: 1/2 + double(0.05)
sig_r = mp.mpf(0.55)                          # where ttest.py's reference is taken: double(0.55)
d = abs(dl(sig_r + 14j) - dl(sig_t + 14j))
S = tilted.setup(0.1, 200, 0.004, 4.0)
A, B = tilted.xi_pair(S, acb(0.05, 14)); r_float = B/A          # exactly as ttest.py
A2, B2 = tilted.xi_pair(S, acb(arb('0.05'), 14)); r_exact = B2/A2  # Horner, exact inputs as in ttest2.py
mp.mp.dps = 40                                                   # ttest.py works at 40 digits
e_float = abs(tomp(r_float, 40) - dl(mp.mpf(0.55) + 1j*mp.mpf(14)))
e_exact = abs(tomp(r_exact, 40) - dl(mp.mpf('0.55') + 1j*mp.mpf(14)))
out['ttest14'] = dict(sigma_gap=mp.nstr(sig_r - sig_t, 4), change_of_ratio_over_gap=mp.nstr(d, 4),
                      ttest_err_float_inputs=mp.nstr(e_float, 4), ttest_err_exact_inputs=mp.nstr(e_exact, 4))
print(json.dumps(out['ttest14']))
json.dump(out, open('extra_check.json', 'w'), indent=1)
