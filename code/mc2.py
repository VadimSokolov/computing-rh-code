import json, numpy as np, mpmath as mp, time
V=json.load(open('volterra.json')); T=np.array(V['t']); B=np.array(V['b']); t0=V['t0']
R=json.load(open('rsboundary.json'))
D=json.load(open('hazard_spec.json')); tt=np.array([r[0] for r in D['rows']]); SS=np.array([float(mp.mpf(r[2])) for r in D['rows']])
g1=14.134725141734693; beta=-np.sqrt(2)*g1; c=3/(2*abs(beta))
m=(T>=0.15)
a=np.mean(B[m]-beta*T[m]-c*np.log(T[m])); res=B[m]-(beta*T[m]+a+c*np.log(T[m]))
cfree=np.polyfit(np.column_stack([T[m],np.log(T[m])]).T[1],B[m]-beta*T[m],1)
print('tail fit a',a,'max residual',np.abs(res).max(),'free log coefficient',cfree[0],'predicted',c)
def sim(TB,BB,tb0,n=200000,dt=2e-5,tmax=0.15,seed=1):
    rng=np.random.default_rng(seed)
    bnd=lambda t: np.where(t<tb0,BB[0],np.interp(t,TB,BB))
    x=np.zeros(n); hit=np.full(n,np.inf); alive=np.ones(n,bool); t=0.0; bp=bnd(0.0)
    for k in range(int(tmax/dt)):
        idx=np.nonzero(alive)[0]
        if idx.size==0: break
        xn=x[idx]+np.sqrt(dt)*rng.standard_normal(idx.size); bn=bnd(t+dt)
        cr=xn>=bn; d0=bp-x[idx]; d1=bn-xn
        cr|=rng.random(idx.size)<np.exp(-2*np.clip(d0,0,None)*np.clip(d1,0,None)/dt)
        hit[idx[cr]]=t+dt/2; alive[idx[cr]]=False; x[idx]=xn; t+=dt; bp=bn
    return np.sort(hit[np.isfinite(hit)])
grid=np.linspace(0.011,0.1,400); Ft=1-np.interp(grid,tt,SS)
out={}
for name,TB,BB,tb in [('volterra',T,B,t0),('tangent',np.array(R['t']),np.array(R['b']),R['t0'])]:
    t1=time.time(); hs=sim(TB,BB,tb); n=200000
    Fe=np.searchsorted(hs,grid)/n; ks=np.abs(Fe-Ft).max()
    out[name]=dict(ks=float(ks),mean=float(hs.mean()),Fe=list(Fe))
    print(name,'KS',ks,'mean',hs.mean(),'sec',time.time()-t1,flush=True)
out['grid']=list(grid); out['Ft']=list(Ft); out['tail']=dict(a=a,c=c,beta=beta,maxres=float(np.abs(res).max()),cfree=float(cfree[0]))
json.dump(out,open('mc2.json','w'))
