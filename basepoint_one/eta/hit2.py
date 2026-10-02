import numpy as np, json, time
from scipy.optimize import brentq
H=json.load(open('hit.json')); t=np.array(H['t']); b=np.array(H['b']); F=np.array(H['F']); S=np.array(H['S']); i0=H['i0']
f=np.gradient(F,t); bt=b.copy()
for k in range(i0+1,len(t)):
    tk=t[k]; dt=tk-t[k-1]; fk=max(f[k],1e-300)
    Hf=lambda x: tk*(x-bt[k-1])/dt-(x-np.exp(np.log(np.sqrt(2*np.pi)*tk**1.5*fk)+x*x/(2*tk)))
    xs=bt[k-1]+np.linspace(-0.5,0.5,201); hv=np.array([Hf(x) for x in xs])
    ok=[q for q in range(200) if np.isfinite(hv[q]) and np.isfinite(hv[q+1]) and hv[q]*hv[q+1]<=0]
    if not ok: print('tangent stops at',tk); bt[k:]=np.nan; break
    q=min(ok,key=lambda q: abs(xs[q]-bt[k-1])); bt[k]=brentq(Hf,xs[q],xs[q+1])
for q in [0.08,0.12,0.2,0.3,0.4]:
    i=np.argmin(abs(t-q)); print('t %.3f exact %.5f tangent %.5f'%(t[i],b[i],bt[i]))
rng=np.random.default_rng(9); n=150000
grid=np.concatenate([np.arange(0,0.5,2e-4),np.arange(0.5,5+1e-9,2e-3)]); bnd=np.interp(grid,t,b,left=b[0])
x=np.zeros(n); alive=np.ones(n,bool); Sv=[]; gaps={}; t0=time.time()
for k in range(1,len(grid)):
    dt=grid[k]-grid[k-1]; idx=np.nonzero(alive)[0]
    xn=x[idx]+np.sqrt(dt)*rng.standard_normal(idx.size); d0=bnd[k-1]-x[idx]; d1=bnd[k]-xn
    cr=(xn>=bnd[k])|(rng.random(idx.size)<np.exp(-2*np.clip(d0,0,None)*np.clip(d1,0,None)/dt))
    alive[idx[cr]]=False; x[idx]=xn; Sv.append(alive.mean())
    for q in [1.0,2.0,5.0]:
        if abs(grid[k]-q)<1e-9: gaps[q]=float(np.sum(bnd[k]-x[alive])/n)
G=grid[1:]; Sv=np.array(Sv); ex=np.interp(G,t,S); ks=np.max(np.abs(Sv-ex))
print('MC sec',time.time()-t0,'KS',ks,'gaps',gaps,'target',0.16508331230553805/np.sqrt(2))
for q in [0.1,0.2,0.4,1,5]:
    i=np.argmin(abs(G-q)); print('t',q,'MC',Sv[i],'exact',ex[i])
H['bt']=bt.tolist(); H['mc']=dict(G=G[::20].tolist(),S=Sv[::20].tolist(),ex=ex[::20].tolist(),ks=float(ks),gaps=gaps)
json.dump(H,open('hit.json','w'))
