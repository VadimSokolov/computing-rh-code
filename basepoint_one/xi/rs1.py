import numpy as np, json, time
from scipy.stats import norm
from scipy.optimize import brentq
g=np.concatenate([np.load(f) for f in ['g_1_400.npy','g_401_1000.npy','g_1001_1700.npy']])
p1=0.0231049931154189707888; p2=0.0000371725992852696862
lam=g**2+0.25; K=200; rest=lam[K:]
tm=np.sum(1/rest)+(p1-np.sum(1/g**2)); tv=np.sum(1/rest**2)+(p2-np.sum(1/g**4))
rng=np.random.default_rng(3); n=300000
Tt=(rng.exponential(size=(n,K))/lam[:K]).sum(1)+rng.gamma(tm*tm/tv,tv/tm,n)
Tt=np.sort(Tt); c=Tt/np.sqrt(2)
# grid: uniform-ish in log time
tg=np.concatenate([np.linspace(0.006,0.05,700)[:-1],np.geomspace(0.05,100,900)])
def FS(x):
    m=Tt<x; d=x-Tt[m]
    F=np.sum(2*norm.sf(c[m]/np.sqrt(d)))/n
    S=(np.sum(2*norm.cdf(c[m]/np.sqrt(d))-1)+np.sum(~m))/n
    return F,S
t0=time.time()
import mpmath as mp
mp.mp.dps=30
xi=lambda z: z*(z-1)/2*mp.pi**(-z/2)*mp.gamma(z/2)*mp.zeta(z)
cdfL=lambda s: mp.mpf(0.5)/xi(1+mp.sqrt(s))/s
F=np.zeros(len(tg)); S=np.zeros(len(tg))
for i,x in enumerate(tg):
    if x<0.03:
        F[i]=float(mp.invertlaplace(cdfL,mp.mpf(x),method='talbot')); S[i]=1-F[i]
    else: F[i],S[i]=FS(x)
print('F,S computed',time.time()-t0, 'S at 1,10,100:',[S[np.argmin(abs(tg-q))]*np.sqrt(q) for q in [1,10,100]])
dF=np.diff(F); dS=-np.diff(S)   # increments of crossing probability (use dS in tail for relative accuracy)
b=np.zeros(len(tg)); i0=np.argmax(F>1e-7); 
b[:i0+1]=np.sqrt(tg[i0])*norm.isf(F[i0]/2)
F0=F[i0]; start=i0
for k in range(start+1,len(tg)):
    tk=tg[k]; comp=S[k]<0.5
    sm=(tg[start:k]+tg[start+1:k+1])/2; w=(dF if not comp else dS)[start:k]
    def G(bk):
        bb=b.copy(); bb[k]=bk; bm=(bb[start:k]+bb[start+1:k+1])/2
        if not comp:
            return norm.sf(bk/np.sqrt(tk))-(np.sum(w*norm.sf((bk-bm)/np.sqrt(tk-sm)))+F0*norm.sf((bk-b[start])/np.sqrt(tk-tg[start])))
        return (np.sum(w*norm.cdf((bk-bm)/np.sqrt(tk-sm)))+F0*norm.cdf((bk-b[start])/np.sqrt(tk-tg[start])))-(norm.cdf(bk/np.sqrt(tk))-S[k])
    xs=b[k-1]+np.linspace(-0.2,0.2,81)*max(1,np.sqrt(tk)); gv=np.array([G(x) for x in xs])
    idx=[q for q in range(80) if gv[q]*gv[q+1]<=0 and gv[q]>0] or [q for q in range(80) if gv[q]*gv[q+1]<=0]
    if not idx: print('fail',k,tk); b=b[:k]; tg=tg[:k]; break
    q=min(idx,key=lambda q: abs(xs[q]-b[k-1])); b[k]=brentq(G,xs[q],xs[q+1],xtol=1e-12)
print('volterra time',time.time()-t0)
binf=0.0230957089661/np.sqrt(2)
for q in [0.009,0.012,0.02,0.05,0.1,1,10,50,100]:
    i=np.argmin(abs(tg-q)); print('t %.3f b %.6f'%(tg[i],b[i]))
print('b_inf predicted',binf)
json.dump(dict(t=tg.tolist(),b=b.tolist(),F=F[:len(tg)].tolist(),S=S[:len(tg)].tolist(),binf=binf,i0=int(i0)),open('rs1.json','w'))
