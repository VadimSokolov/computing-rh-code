# Certificate for the cutoff in Theorem fm:thm:zc(i): 3 s_16 < Var X < 3 s_17, where s_K = sum_{k<=K} gamma_k^{-2} over the
# first K zeros on the critical line and Var X = (log xi)''(1/2) is the variance of the Polya variable, E e^{sX} = xi(1/2+s)/xi(1/2).
# Ball arithmetic (python-flint).  The ordinates are Arb's certified zeros (acb.zeta_zero, which isolates them and counts them
# by Turing's method), so no decimal table enters.  Var X is twice the coefficient of x^2 in the Taylor series of
#   log xi(1/2 + x) = log s + log(1-s) - (s/2) log pi + log Gamma(s/2) + log(-zeta(s)) + const,  s = 1/2 + x,
# where log(1-s) and log(-zeta(s)) differ from log(s-1) and log zeta(s) by constants, so the series stays real.
# Writes fm_cutoff.json.  Run on Hopper from ferro/: bash ../misc/tools/hopper_run.sh -c 1 -t 20 -g fm_cutoff.json fm_cutoff.py
import json, time
from flint import arb, acb, acb_series, ctx

ctx.prec = 200
t0 = time.time()
gam = [acb.zeta_zero(k).imag for k in range(1, 18)]
s = [arb(0)]
for g in gam:
    s.append(s[-1] + 1 / g ** 2)

x = acb_series([acb(arb(1) / 2), acb(1)], prec=3)
LOGPI = arb.pi().log()
lx = x.log() + (1 - x).log() - x * (LOGPI / 2) + (x / 2).lgamma() + (-x.zeta()).log()
var = 2 * lx.coeffs()[2].real

m16 = var - 3 * s[16]                      # must be positive: the Maxwell construction works for K <= 16
m17 = 3 * s[17] - var                      # must be positive: it fails at K = 17
lower = lambda b: float(b.mid() - b.rad())
out = {'prec_bits': ctx.prec,
       'gamma_1': gam[0].str(30), 'gamma_16': gam[15].str(30), 'gamma_17': gam[16].str(30),
       'var_X': var.str(25), '3s_16': (3 * s[16]).str(25), '3s_17': (3 * s[17]).str(25),
       'margin_16_lower': lower(m16), 'margin_17_lower': lower(m17),
       'two_beta_16': (var - 2 * s[16]).str(15), 's_16': s[16].str(15),
       'ok': lower(m16) > 0 and lower(m17) > 0,
       'seconds': round(time.time() - t0, 2)}
print(json.dumps(out, indent=1))
json.dump(out, open('fm_cutoff.json', 'w'), indent=1)
