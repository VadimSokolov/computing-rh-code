# PATCHED EXCERPT of code/plots.py, not the authors' file: its imports, the line that loads plotcurves.npz, and
# the section "mass derivative check and prime table" (code/plots.py lines 1-5, 7 and 101-109) that writes
# data/extra.json, with the block that this folder proposes for the kappa_cell item inserted after
# out=dict(XipXi200=...). The figures of code/plots.py are left out. The output extra.json must hold the archived
# keys unchanged (XipXi200, pred_coeff, prime) and two new ones: kappa_cell, the first order cell constant
# (1/pi) sum_k |Xi'/Xi(b_k) - Xi'/Xi(b_{k-1})| of ch07 and ch15, and cell_leaks, the 79 first order cell errors
# (mu_eps(C_k)-1)/eps = -(Xi'/Xi(b_k) - Xi'/Xi(b_{k-1}))/pi, with b_0 = 0, b_k = (gamma_k+gamma_{k+1})/2, b_79 = 200.
# Inputs: polya_flint.py (code/), plotcurves.npz and results_prime.json of the earlier rerun (read only, copied
# into the run folder by the job). Run: python3 plots_extra_patched.py
import numpy as np, json
from flint import acb
from polya_flint import setup, xi_pair
C = np.load('plotcurves.npz'); Z = C['zeros']
# mass derivative check and prime table
S=setup(420,200,3)
a,b=xi_pair(S,acb(0,200)); XiprimeOverXi=float((-(b.imag))/a.real)  # d/dtheta Xi = -Im xi'
out=dict(XipXi200=XiprimeOverXi, pred_coeff=-XiprimeOverXi/np.pi)
# ---- new: first order cell constant of ch07 and ch15, Xi'/Xi at the cell boundaries by the Polya route ----
bc=[0.0]+[float((Z[k]+Z[k+1])/2) for k in range(78)]+[200.0]
Lp=[]
for t in bc:
    xa,xb=xi_pair(S,acb(0,t)); Lp.append(float((-(xb.imag))/xa.real))
leaks=[-(Lp[k+1]-Lp[k])/np.pi for k in range(79)]
out['kappa_cell']=float(np.sum(np.abs(leaks))); out['cell_leaks']=leaks
# ---- end of the new block ----
pr=json.load(open('results_prime.json')); new=[]
for q in pr:
    a,b=xi_pair(S,acb(q['alpha']-0.5,q['b'])); q['polya']=float((b/a).real); new.append(q)
out['prime']=new
json.dump(out,open('extra.json','w'),indent=1); print(out['XipXi200'], out['pred_coeff'], out['kappa_cell'])
