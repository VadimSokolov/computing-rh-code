import numpy as np, mpmath as mp, json
mp.mp.dps=25
xi=lambda s: s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
dl=lambda s: mp.diff(lambda x: mp.log(xi(x)), s)
N=10**7
sieve=np.ones(N+1,bool); sieve[:2]=False
for p in range(2,int(N**0.5)+1):
    if sieve[p]: sieve[p*p::p]=False
lam=np.zeros(N+1)
for p in np.nonzero(sieve)[0]:
    lp=np.log(p); pk=p
    while pk<=N: lam[pk]=lp; pk*=p
n=np.nonzero(lam)[0].astype(float); L_=lam[n.astype(int)]; ln=np.log(n)
out=[]
for al in [1.0,0.9,0.75,0.6,0.5]:
    a=1-al
    b=float(mp.re(dl(mp.mpf(al)))) if al!=0.5 else 0.0
    for th in [5.0,14.134725]:
        s=mp.mpc(al,th)
        exact=float(mp.re(dl(s)))
        A=float(mp.re(1/s)-1/mp.mpf(al)+(mp.re(mp.digamma(s/2))-mp.digamma(mp.mpf(al)/2))/2)
        for X in [1e3,1e4,1e5,1e6,1e7]:
            m=n<=X; Lg=np.log(X)
            pr=np.sum(L_[m]*n[m]**(-al)*2*np.sin(th*ln[m]/2)**2)
            if a==0: pole=Lg-np.sin(th*Lg)/th
            else:
                pole=(np.exp(a*Lg)-1)/a-np.real((np.exp((a+1j*th)*Lg)-1)/(a+1j*th))
            tot=b+A+pr-pole
            out.append(dict(alpha=al,theta=th,X=X,prime=float(pr),pole=float(pole),total=float(tot),exact=exact,err=float(tot-exact)))
            print('alpha %.2f theta %6.3f X %.0e primes %10.3f pole %10.3f total %9.5f exact %9.5f err %9.2e'%(al,th,X,pr,pole,tot,exact,tot-exact),flush=True)
json.dump(out,open('below.json','w'))
