import numpy as np, json, mpmath as mp
from scipy.stats import norm
from scipy.optimize import brentq
from scipy.interpolate import CubicSpline
H=json.load(open('hazard_spec.json'))['rows']; Cg=json.load(open('clockgrid.json'))['rows']
t1=np.array([r[0] for r in H]); f1=np.array([float(mp.mpf(r[1])) for r in H]); S1=np.array([float(mp.mpf(r[2])) for r in H])
t2=np.array([r['t'] for r in Cg]); f2=np.array([float(mp.mpf(r['f'])) for r in Cg]); S2=np.array([float(mp.mpf(r['S'])) for r in Cg])
m2=(t2>=0.006)&(t2<0.00999)
t=np.concatenate([t2[m2],t1]); f=np.concatenate([f2[m2],f1]); S=np.concatenate([S2[m2],S1])
o=np.argsort(t); t,f,S=t[o],np.maximum(f[o],1e-300),S[o]
lf=CubicSpline(t,np.log(f)); F=lambda x: 1-np.interp(x,t,S)
# Volterra equation of the first kind (Peskir): P(B_t>=b(t)) = int_0^t f(s) P(B_t>=b(t) | B_s=b(s)) ds
t0=0.009; T=0.3; n=2900; tg=np.linspace(t0,T,n+1); dt=tg[1]-tg[0]
fg=np.exp(lf(tg)); F0=F(t0)
b=np.zeros(n+1); b[0]=np.sqrt(t0)*norm.isf(max(F0,1e-15)/2)
mass0=F0   # mass crossed before t0, placed at t0 with the level b[0]
for k in range(1,n+1):
    tk=tg[k]
    Sk=np.interp(tk,t,S)
    comp = Sk<0.5
    def G(bk):
        bb=b.copy(); bb[k]=bk
        sm=(tg[:k]+tg[1:k+1])/2; bm=(bb[:k]+bb[1:k+1])/2; fm=np.exp(lf(sm))
        if not comp:
            lhs=norm.sf(bk/np.sqrt(tk))
            ker=norm.sf((bk-bm)/np.sqrt(tk-sm))
            rhs=np.sum(fm*ker)*dt+mass0*norm.sf((bk-b[0])/np.sqrt(tk-t0))
            return lhs-rhs
        # complementary form, free of cancellation in the right tail
        lhs=norm.cdf(bk/np.sqrt(tk))-Sk
        ker=norm.cdf((bk-bm)/np.sqrt(tk-sm))
        rhs=np.sum(fm*ker)*dt+mass0*norm.cdf((bk-b[0])/np.sqrt(tk-t0))
        return rhs-lhs
    found=False
    for W in [0.05,0.2,0.8]:
        xs=b[k-1]+np.linspace(-W,W,41); gv=np.array([G(x) for x in xs])
        idx=[q for q in range(40) if gv[q]>0 and gv[q+1]<=0]
        if idx: lo,hi=xs[idx[0]],xs[idx[0]+1]; found=True; break
    if not found:
        print('fail',k,tk); b=b[:k]; tg=tg[:k]; break
    b[k]=brentq(G,lo,hi,xtol=1e-12) if G(hi)!=0 else hi
R=json.load(open('rsboundary.json'))
for q in [0.011,0.02,0.05,0.1,0.2,0.3]:
    i=np.argmin(abs(tg-q)); j=np.argmin(abs(np.array(R['t'])-q))
    print(q, round(b[i],5), round(R['b'][j],5) if abs(R['t'][j]-q)<1e-3 else None)
sl=np.polyfit(tg[tg>0.2],b[tg>0.2],1)[0]; print('asymptotic slope',sl,'target',-np.sqrt(2)*14.134725)
json.dump(dict(t=list(tg),b=list(b),t0=t0,F0=F0,slope=sl),open('volterra.json','w'))
