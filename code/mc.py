import json, numpy as np, mpmath as mp, time
R=json.load(open('rsboundary.json')); T=np.array(R['t']); B=np.array(R['b']); t0=R['t0']
D=json.load(open('hazard_spec.json')); tt=np.array([r[0] for r in D['rows']]); SS=np.array([float(mp.mpf(r[2])) for r in D['rows']])
def bnd(t): return np.where(t<t0, B[0], np.interp(t,T,B))
rng=np.random.default_rng(1); n=200000; dt=2e-5; tmax=0.15; steps=int(tmax/dt)
x=np.zeros(n); hit=np.full(n,np.inf); alive=np.ones(n,bool); t=0.0; b_prev=bnd(0.0)
t1=time.time()
for k in range(steps):
    idx=np.nonzero(alive)[0]
    if idx.size==0: break
    xn=x[idx]+np.sqrt(dt)*rng.standard_normal(idx.size); b_new=bnd(t+dt)
    crossed=xn>=b_new
    # Brownian bridge crossing between grid points
    d0=b_prev-x[idx]; d1=b_new-xn
    pb=np.exp(-2*np.clip(d0,0,None)*np.clip(d1,0,None)/dt)
    crossed|= rng.random(idx.size)<pb
    hit[idx[crossed]]=t+dt*0.5; alive[idx[crossed]]=False; x[idx]=xn; t+=dt; b_prev=b_new
print('sim sec',time.time()-t1, 'unhit',alive.sum())
hs=np.sort(hit[np.isfinite(hit)])
grid=np.linspace(0.011,0.1,400); Fe=np.searchsorted(hs,grid)/n; Ft=1-np.interp(grid,tt,SS)
ks=np.max(np.abs(Fe-Ft)); print('KS', ks, 'MC noise ~', 1/np.sqrt(n))
print('mean sim', hs.mean(), 'true E T 0.0231050')
json.dump(dict(grid=list(grid),Fe=list(Fe),Ft=list(Ft),ks=ks,n=n,mean=float(hs.mean())),open('mc.json','w'))
