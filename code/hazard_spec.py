import mpmath as mp, json, numpy as np
from clock import *
mp.mp.dps=50
P1=json.load(open('part1.json'))
g=[mp.mpf(r) for r in P1['roots']]
R=[]
for gk in g:
    _,B = xi_pair(S, to_acb(mp.mpc(0,gk)))
    R.append(mp.re(2j*gk*to_mp(acb(XI0))/to_mp(B)))
def f_s(t): return mp.fsum(R[k]*mp.exp(-g[k]**2*t) for k in range(len(g)))
def S_s(t): return mp.fsum(R[k]/g[k]**2*mp.exp(-g[k]**2*t) for k in range(len(g)))
ts=np.round(np.arange(0.010,0.4001,0.0005),5)
rows=[]
for t in ts:
    t=mp.mpf(str(t)); fv=f_s(t); Sv=S_s(t); rows.append([float(t), mp.nstr(fv,30), mp.nstr(Sv,30), mp.nstr(fv/Sv,30)])
# moments check from residues: E T = sum R_k/g^4 ? (int t f = sum R/g^4)
ET=mp.fsum(R[k]/g[k]**4 for k in range(len(g))); m1=mp.fsum(1/gk**2 for gk in g)
json.dump(dict(rows=rows,R=[mp.nstr(r,25) for r in R],g=[mp.nstr(x,40) for x in g],ET_res=mp.nstr(ET,15),sum_inv_g2_79=mp.nstr(m1,15)),open('hazard_spec.json','w'))
print('ET from residues', ET, ' first residues', [mp.nstr(r,8) for r in R[:4]])
for r in rows[::60]: print(r[0], r[3][:22])
