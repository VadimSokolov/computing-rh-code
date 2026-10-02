import numpy as np, json, time, sys, mpmath as mp
from scipy.stats import norm
from scipy.optimize import brentq
al=float(sys.argv[1]); a=al-0.5
mp.mp.dps=25
xi=lambda z: z*(z-1)/2*mp.pi**(-z/2)*mp.gamma(z/2)*mp.zeta(z)
g=np.concatenate([np.load(f) for f in ['g_1_400.npy','g_401_1000.npy','g_1001_1700.npy']])
p1=0.0231049931154189707888; p2=0.0000371725992852696862
lam=g**2+a*a; K=200; rest=lam[K:]
tm=np.sum(1/rest)+(p1-np.sum(1/g**2)); tv=np.sum(1/rest**2)+(p2-np.sum(1/g**4))
rng=np.random.default_rng(3); n=300000
Tt=(rng.exponential(size=(n,K))/lam[:K]).sum(1)+rng.gamma(tm*tm/tv,tv/tm,n)
Ts=np.sort(Tt); c=np.sqrt(2)*a*Ts
xa=xi(mp.mpf(al)); cdfL=lambda s: xa/xi(al+mp.sqrt(s))/s
tg=np.concatenate([np.linspace(0.006,0.05,700)[:-1],np.geomspace(0.05,100,900)])
F=np.zeros(len(tg)); S=np.zeros(len(tg)); t0=time.time()
for i,x in enumerate(tg):
    if x<0.03: F[i]=float(mp.invertlaplace(cdfL,mp.mpf(x),method='talbot')); S[i]=1-F[i]
    else:
        m=Ts<x; d=x-Ts[m]
        F[i]=np.sum(2*norm.sf(c[m]/np.sqrt(d)))/n; S[i]=(np.sum(2*norm.cdf(c[m]/np.sqrt(d))-1)+np.sum(~m))/n
dF=np.diff(F); dS=-np.diff(S)
i0=int(np.argmax(F>1e-7)); b=np.zeros(len(tg)); b[:i0+1]=np.sqrt(tg[i0])*norm.isf(F[i0]/2); F0=F[i0]
for k in range(i0+1,len(tg)):
    tk=tg[k]; comp=S[k]<0.5
    sm=(tg[i0:k]+tg[i0+1:k+1])/2; w=(dS if comp else dF)[i0:k]
    def G(bk):
        bb=b.copy(); bb[k]=bk; bm=(bb[i0:k]+bb[i0+1:k+1])/2
        if not comp: return norm.sf(bk/np.sqrt(tk))-(np.sum(w*norm.sf((bk-bm)/np.sqrt(tk-sm)))+F0*norm.sf((bk-b[i0])/np.sqrt(tk-tg[i0])))
        return (np.sum(w*norm.cdf((bk-bm)/np.sqrt(tk-sm)))+F0*norm.cdf((bk-b[i0])/np.sqrt(tk-tg[i0])))-(norm.cdf(bk/np.sqrt(tk))-S[k])
    xs=b[k-1]+np.linspace(-0.2,0.2,81)*max(1,np.sqrt(tk)); gv=np.array([G(x) for x in xs])
    idx=[q for q in range(80) if gv[q]*gv[q+1]<=0 and gv[q]>0] or [q for q in range(80) if gv[q]*gv[q+1]<=0]
    if not idx: print('fail',k,tk); b=b[:k]; tg=tg[:k]; F=F[:k]; S=S[:k]; break
    q=min(idx,key=lambda q: abs(xs[q]-b[k-1])); b[k]=brentq(G,xs[q],xs[q+1],xtol=1e-12)
ba=float(mp.diff(lambda z: mp.log(xi(z)),mp.mpf(al)))
dP=np.diff(1-S); bm=(b[1:]+b[:-1])/2; Eb=np.sum(bm[i0:]*dP[i0:])+b[i0]*(1-S[i0])+b[-1]*S[-1]
print('alpha',al,'time',time.time()-t0,'b at 0.1,1,10,100:',[round(float(np.interp(q,tg,b)),5) for q in [0.1,1,10,100]],'E b(X)',Eb,'target',ba/np.sqrt(2),'sqrt(100)S(100)',S[-1]*10,'pred',ba/np.sqrt(np.pi))
json.dump(dict(alpha=al,t=tg.tolist(),b=b.tolist(),S=S.tolist(),Eb=float(Eb),target=ba/np.sqrt(2),ba=ba),open('bnd_%s.json'%al,'w'))
