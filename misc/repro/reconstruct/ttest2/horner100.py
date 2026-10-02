# Reconstructed diagnostic, not the authors' script. Section sec:ch4:balls of ch/ch04.tex quotes a Horner
# experiment (2001 nodes, delta 0.1, 200 bits: radius growth of about 1.3 per step, final ball at theta 100
# of radius 10^109 with a correct midpoint, direct evaluation giving radius 10^-52 at the same point) that
# ch/appA.tex attributes to the archived code/ttest.py. But code/ttest.py evaluates only at the heights
# 14, 200, 500, 1000 and 2000 and prints only the radius of the real part of the ratio B/A. This script
# evaluates the archived code/tilted.py (copied beside it) with the parameters of ttest.py/ttest2.py for
# 2001 nodes (delta 0.1, 200 bits, h 0.004, xmax 4) at theta = 100 (and 14, 200 for comparison) through
# both evaluators, and replays Horner's rule step by step with the operations of Arb's
# _acb_poly_evaluate_horner (u = f[i] + u*x) to measure the growth of the radius per step. Writes horner100.json.
import math, json
from flint import acb, arb, ctx
import tilted

def lg(x):
    x = float(x)
    return round(math.log10(x), 2) if x > 0 else None
def rad(z):  # the larger of the real and imaginary radii of an acb
    return max(float(z.real.rad()), float(z.imag.rad()))

S = tilted.setup(0.1, 200, 0.004, 4.0)
J = S['J']; h = S['h']; c = S['P'].coeffs()
res = dict(nodes=2*J+1, prec=ctx.prec)
k2 = J + 500                                     # the coefficient at x = 2
res['log10_abs_coeff_at_x2'] = lg(abs(c[k2]).mid())
res['log10_max_abs_coeff'] = lg(max(float(abs(ck).mid()) for ck in c))
for th in [100, 14, 200]:
    s = acb(arb('0.05'), th)
    E = (s*h).exp()
    PE = S['P'](E)                               # Arb's Horner evaluation, as in tilted.xi_pair
    Ah, Bh = tilted.xi_pair(S, s)
    Ad, Bd = tilted.xi_pair_direct(S, s)
    u = c[-1]; rads = [rad(u)]
    for i in range(len(c) - 2, -1, -1):
        u = c[i] + u*E
        rads.append(rad(u))
    same = (u.real.mid() == PE.real.mid()) and (u.imag.mid() == PE.imag.mid()) and rad(u) == rad(PE)
    wrap = math.exp(0.05*0.004)*(abs(math.cos(th*0.004)) + abs(math.sin(th*0.004)))
    last = rads[-1501:]
    g1500 = (last[-1]/last[0])**(1/1500) if last[0] > 0 else None
    ratios = sorted(rads[i+1]/rads[i] for i in range(len(rads)-1501, len(rads)-1) if rads[i] > 0)
    rd = Bd/Ad
    try:
        rh = Bh/Ah; rh_rad = lg(rh.real.rad()) if rh.real.rad().is_finite() else 'inf'
    except Exception as ex:
        rh_rad = repr(ex)
    row = dict(theta=th, wrap_factor=round(wrap, 4), wrap_factor_1500=lg(wrap**1500),
               replay_equals_arb_horner=bool(same),
               growth_per_step_geomean_last1500=round(g1500, 4) if g1500 else None,
               growth_per_step_median_last1500=round(ratios[len(ratios)//2], 4) if ratios else None,
               log10_rad_P_horner=lg(rad(PE)), log10_abs_P_direct=lg(abs(Ad/((s*acb(0, S['eta'])).exp() * E**(-J))).mid()),
               log10_rad_A_horner=lg(rad(Ah)), log10_rad_A_direct=lg(rad(Ad)), log10_abs_A=lg(abs(Ad).mid()),
               horner_mid_rel_err=lg(abs(Ah.mid() - Ad.mid()).mid() / abs(Ad).mid()),
               log10_rad_re_ratio_direct=lg(rd.real.rad()), log10_rad_re_ratio_horner=rh_rad,
               log10_rad_B_direct=lg(rad(Bd)), log10_abs_B=lg(abs(Bd).mid()))
    res[f'theta_{th}'] = row
    print(json.dumps(row), flush=True)
json.dump(res, open('horner100.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in res.items() if not k.startswith('theta')}))
