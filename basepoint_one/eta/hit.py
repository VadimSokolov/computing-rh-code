import numpy as np, mpmath as mp, json, time
from scipy.stats import norm
from scipy.optimize import brentq
mp.mp.dps=25
chi=[0,1,0,0,0,-1,0,-1,0,0,0,1]
Lam=lambda s: (12/mp.pi)**(s/2)*mp.gamma(s/2)*mp.dirichlet(s,chi)
g=np.load('g12.npy'); B=json.load(open('base.json')); b1=B['b']; p1=B['p1']; p2=B['p2']
lam=g**2+0.25
tm=p1-np.sum(1/g**2); tv=p2-np.sum(1/g**4)
rng=np.random.default_rng(21); n=300000
Tt=np.zeros(n)
for c0 in range(0,n,50000):
    Tt[c0:c0+50000]=(rng.exponential(size=(50000,len(g)))/lam).sum(1)+rng.gamma(tm*tm/tv,tv/tm,50000)
print('E Ttilde MC',Tt.mean(),' b1',b1)
Z=rng.standard_normal(n); X=Tt+Tt**2/(2*Z**2)
L1=Lam(mp.mpf(1)); lt=[]
for s in [0.1,0.5,1,5,20]:
    emp=np.mean(np.exp(-s*X)); se=np.std(np.exp(-s*X))/np.sqrt(n); ex=float(L1/Lam(1+mp.sqrt(s)))
    lt.append(dict(s=s,emp=float(emp),se=float(se),exact=ex)); print('s',s,'MC %.6f +- %.6f exact %.6f'%(emp,se,ex))
tails=[]
for x in [1e2,1e4,1e6]:
    pr=np.mean(X>x); tails.append(dict(x=x,scaled=float(pr*np.sqrt(x)))); print('x',x,'sqrt(x)P',pr*np.sqrt(x),'pred',b1/np.sqrt(np.pi))
# exact boundary
Ts=np.sort(Tt); c=Ts/np.sqrt(2)
tg=np.concatenate([np.linspace(0.03,0.4,800)[:-1],np.geomspace(0.4,1000,900)])
cdfL=lambda s: L1/Lam(1+mp.sqrt(s))/s
F=np.zeros(len(tg)); S=np.zeros(len(tg)); t0=time.time()
for i,x in enumerate(tg):
    if x<0.12: F[i]=float(mp.invertlaplace(cdfL,mp.mpf(x),method='talbot')); S[i]=1-F[i]
    else:
        m=Ts<x; d=x-Ts[m]
        F[i]=np.sum(2*norm.sf(c[m]/np.sqrt(d)))/n; S[i]=(np.sum(2*norm.cdf(c[m]/np.sqrt(d))-1)+np.sum(~m))/n
print('F time',time.time()-t0,'check continuity at 0.12:',F[np.searchsorted(tg,0.12)-1],F[np.searchsorted(tg,0.12)])
dF=np.diff(F); dS=-np.diff(S)
i0=int(np.argmax(F>1e-7)); b=np.zeros(len(tg)); b[:i0+1]=np.sqrt(tg[i0])*norm.isf(F[i0]/2); F0=F[i0]
for k in range(i0+1,len(tg)):
    tk=tg[k]; comp=S[k]<0.5
    sm=(tg[i0:k]+tg[i0+1:k+1])/2; w=(dS if comp else dF)[i0:k]
    def G(bk):
        bb=b.copy(); bb[k]=bk; bm=(bb[i0:k]+bb[i0+1:k+1])/2
        if not comp: return norm.sf(bk/np.sqrt(tk))-(np.sum(w*norm.sf((bk-bm)/np.sqrt(tk-sm)))+F0*norm.sf((bk-b[i0])/np.sqrt(tk-tg[i0])))
        return (np.sum(w*norm.cdf((bk-bm)/np.sqrt(tk-sm)))+F0*norm.cdf((bk-b[i0])/np.sqrt(tk-tg[i0])))-(norm.cdf(bk/np.sqrt(tk))-S[k])
    xs=b[k-1]+np.linspace(-0.3,0.3,121); gv=np.array([G(x) for x in xs])
    idx=[q for q in range(120) if gv[q]*gv[q+1]<=0 and gv[q]>0] or [q for q in range(120) if gv[q]*gv[q+1]<=0]
    if not idx: print('fail',k,tk); b=b[:k]; tg=tg[:k]; F=F[:k]; S=S[:k]; break
    q=min(idx,key=lambda q: abs(xs[q]-b[k-1])); b[k]=brentq(G,xs[q],xs[q+1],xtol=1e-12)
print('volterra done',time.time()-t0)
for q in [tg[i0],0.08,0.12,0.2,0.4,1,10,100,1000]:
    i=np.argmin(abs(tg-q)); print('t %.3f b %.5f S %.5f'%(tg[i],b[i],S[i]))
dP=np.diff(1-S); bm=(b[1:]+b[:-1])/2
Eb=np.sum(bm[i0:]*dP[i0:])+b[i0]*(1-S[i0])+b[-1]*S[-1]
print('E b(X1)',Eb,' target b1/sqrt2',b1/np.sqrt(2))
json.dump(dict(lt=lt,tails=tails,ET=float(Tt.mean()),t=tg.tolist(),b=b.tolist(),F=F.tolist(),S=S.tolist(),i0=i0,Eb=float(Eb)),open('hit.json','w'))
