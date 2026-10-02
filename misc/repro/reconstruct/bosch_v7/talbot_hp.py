# Supplementary check, not the authors' code (misc/repro/reconstruct/bosch_v7). Block 6 of
# code_bosch/verify_bosch_xi_v7.py inverts phi(s) = xi(1/2)/xi(1/2+sqrt s) by Talbot's method at 40 digits; at
# t=0.003 that gives a negative value (output line 60), because the density there is about 1e-72. This script runs
# the same call, mp.invertlaplace(phi, t, method='talbot'), at 150 digits, as a check of the residue series of
# residue_hp_collect.py that does not use the zeros. Usage: python talbot_hp.py T
import sys, time, json
import mpmath as mp
mp.mp.dps=150
def xi(z): return z*(z-1)/2*mp.pi**(-z/2)*mp.gamma(z/2)*mp.zeta(z)
half=mp.mpf(1)/2; x0=xi(half)
phi=lambda s: x0/xi(half+mp.sqrt(s))
ts=sys.argv[1]; t=mp.mpf(ts)
t0=time.time(); inv=mp.invertlaplace(phi,t,method='talbot')
s="t=%s: Talbot inverse Laplace at 150 digits %s (%.0fs)"%(ts,mp.nstr(inv,15),time.time()-t0)
print(s,flush=True)
json.dump({"t":ts,"talbot":mp.nstr(inv,30),"seconds":round(time.time()-t0)},open("talbot_%s.json"%ts,"w"))
