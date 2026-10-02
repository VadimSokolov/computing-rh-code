import mpmath as mp, json, numpy as np, time
from clock import *
P1=json.load(open('part1.json'))
g=[mp.mpf(r) for r in P1['roots']]
R=[]
for gk in g:
    _,B = xi_pair(S, to_acb(mp.mpc(0,gk)))
    R.append(mp.re(2j*gk*to_mp(acb(XI0))/to_mp(B)))
def spec(t): return mp.fsum(R[k]*mp.exp(-g[k]**2*t) for k in range(len(g)))
def spec_surv(t): return mp.fsum(R[k]/g[k]**2*mp.exp(-g[k]**2*t) for k in range(len(g)))
out=dict(R=[mp.nstr(r,20) for r in R[:10]], rows=[])
ts=list(np.round(np.concatenate([np.arange(0.004,0.02,0.0005),np.arange(0.02,0.1,0.002),np.arange(0.1,0.301,0.01)]),5))
t0=time.time()
for t in ts:
    t=mp.mpf(str(t)); f=dens(t); Sv=surv(t)
    row=dict(t=float(t), f=mp.nstr(f,25), S=mp.nstr(Sv,25))
    if t>=0.012: row['f_spec']=mp.nstr(spec(t),25); row['S_spec']=mp.nstr(spec_surv(t),25)
    out['rows'].append(row)
json.dump(out,open('clockgrid.json','w'),indent=1)
print('done',len(ts),time.time()-t0)
for r in out['rows'][::10]: print(r['t'], r['f'][:14], r.get('f_spec','')[:14], r['S'][:14], r.get('S_spec','')[:14])
