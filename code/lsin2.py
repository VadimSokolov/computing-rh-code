exec(open('lfun.py').read().split('res={}')[0])
mp.mp.dps=30
import sympy
from scipy.optimize import minimize_scalar
primes=list(sympy.primerange(2,NMAX+1))
def lam_f(name):
    out={}
    for p in primes:
        if name=='zeta': ap=None
        lp=mp.log(p); k=1; pk=p
        if name=='Delta':
            a=mp.mpf(tau[p])/mp.mpf(p)**mp.mpf(5.5); s0,s1=mp.mpf(2),a
        while pk<=NMAX:
            if name=='zeta': c=lp
            elif name=='chi12': c=chi12(p)**k*lp
            else:
                c=s1*lp; s0,s1=s1,a*s1-s0
            out[pk]=c; k+=1; pk*=p
    return out
setups={'zeta':(Phi,0.5,3.0),'chi12':(K12,2,5.0),'Delta':(KD,1,4.0)}
sigma=mp.mpf(2); eps=sigma-mp.mpf(0.5)
th=np.linspace(0.0,40,161)
res={}
for name,(K,scale,U) in setups.items():
    h=mp.mpf(1)/80; us,Ks=nodes(K,h,U); uKs=[u/scale*k for u,k in zip(us,Ks)]
    def dlog(w): return Lam(us,uKs,h,w,scale)/Lam(us,Ks,h,w,scale)
    L=lam_f(name)
    D=float(mp.re(dlog(eps)))
    if name=='zeta':
        arch=lambda s: mp.re(1/s+1/(s-1)-mp.log(mp.pi)/2+mp.digamma(s/2)/2)
        pole=lambda s: mp.re(1/(s-1))
    elif name=='chi12':
        arch=lambda s: mp.re(mp.log(12/mp.pi)/2+mp.digamma(s/2)/2); pole=lambda s:0
    else:
        arch=lambda s: mp.re(-mp.log(2*mp.pi)+mp.digamma(s+mp.mpf(11)/2)); pole=lambda s:0
    rows=[]
    for t in th:
        s=sigma+1j*mp.mpf(t)
        ex=float(mp.re(dlog(eps+1j*mp.mpf(t))))
        A=float(arch(s)-arch(sigma)-(pole(s)-pole(sigma)))
        Pl=float(pole(s)-pole(sigma))
        Pr=float(mp.fsum(c*mp.mpf(n)**(-sigma)*2*mp.sin(mp.mpf(t)*mp.log(n)/2)**2 for n,c in L.items()))
        rows.append((float(t),ex,D,A,Pl,Pr,D+A+Pl+Pr))
    R=np.array(rows)
    err=np.abs(R[:,1]-R[:,6]).max()
    # Extrema on [0,40]. The grid th, of step 0.25, is the grid of the figure. The extrema of the prime piece are located on
    # a grid of step 0.001 and refined by Newton's method on its derivative, and the minima of the left side on a grid of
    # step 0.05, refined by Brent's method; the endpoints are included. The pole piece decreases, so its minimum is at 40.
    wl=[(c*mp.mpf(n)**(-sigma),mp.log(n)) for n,c in L.items()]
    def Pk(t,k):   # the prime piece (k=0) and its first two derivatives in theta
        t=mp.mpf(t)
        if k==0: return mp.fsum(w*2*mp.sin(t*l/2)**2 for w,l in wl)
        if k==1: return mp.fsum(w*l*mp.sin(t*l) for w,l in wl)
        return mp.fsum(w*l*l*mp.cos(t*l) for w,l in wl)
    Wf=np.array([float(w) for w,l in wl]); Lf=np.array([float(l) for w,l in wl])
    tf=np.linspace(0.0,40,40001)
    Pf=np.concatenate([(2*np.sin(np.outer(tf[i:i+1000],Lf)/2)**2)@Wf for i in range(0,len(tf),1000)])
    cand=[(Pk(0,0),mp.mpf(0)),(Pk(40,0),mp.mpf(40))]
    for i in range(1,len(tf)-1):
        if (Pf[i]-Pf[i-1])*(Pf[i+1]-Pf[i])<=0:
            t=mp.mpf(tf[i])
            for it in range(50):
                d=Pk(t,1)/Pk(t,2); t-=d
                if abs(d)<mp.mpf(10)**-25: break
            assert abs(t-tf[i])<0.002 and 0<t<40, (name,tf[i],t)
            cand.append((Pk(t,0),t))
    Pmin=min(cand,key=lambda x:x[0]); Pmax=max(cand,key=lambda x:x[0])
    lhs=lambda t: float(mp.re(dlog(eps+1j*mp.mpf(t))))
    tg=np.linspace(0.0,40,801); eg=[lhs(t) for t in tg]
    loc=[]
    for i in range(1,len(tg)-1):
        if eg[i]<=eg[i-1] and eg[i]<=eg[i+1]:
            r=minimize_scalar(lhs,bounds=(tg[i-1],tg[i+1]),method='bounded',options=dict(xatol=1e-9)); loc.append((float(r.x),float(r.fun)))
    Emin=min([(0.0,eg[0]),(40.0,eg[-1])]+loc,key=lambda x:x[1])
    res[name]=dict(rows=R.tolist(),D=D,err=float(err),minP=float(Pmin[0]),maxP=float(Pmax[0]),thetaminP=float(Pmin[1]),thetamaxP=float(Pmax[1]),minPole=float(R[:,4].min()),minex=Emin[1],thetaminex=Emin[0],localminex=loc,neg_prime_weights=int(sum(1 for c in L.values() if c<0)),n_weights=len(L))
    print(name,'drift %.4f'%D,'max |exact-decomp| %.2e'%err,'prime piece range [%.6f, %.6f] at theta %.4f and %.4f'%(Pmin[0],Pmax[0],Pmin[1],Pmax[1]),'(on the grid of the figure [%.6f, %.6f])'%(R[:,5].min(),R[:,5].max()),'pole min %.6f'%R[:,4].min(),'min exact %.6f at theta %g, smallest interior local minimum %s'%(Emin[1],Emin[0],min(loc,key=lambda x:x[1]) if loc else None),'negative weights %d of %d'%(res[name]['neg_prime_weights'],len(L)),flush=True)
json.dump(res,open('lsin2.json','w'))
