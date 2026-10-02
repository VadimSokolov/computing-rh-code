import numpy as np, mpmath as mp, json
mp.mp.dps=25
xi=lambda s: s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
b1=-mp.log(mp.pi)/2+mp.digamma(mp.mpf(1.5))/2+mp.euler
print('b1 =',mp.nstr(b1,15),' xi\'/xi(1)=sum 1/rho =',mp.nstr(1+mp.euler/2-mp.log(4*mp.pi)/2,15))
# archimedean Levy density e^{-3x}/(x(1-e^{-2x}))
arch=lambda s: mp.quad(lambda x:(mp.e**(-s*x)-1+s*x)*mp.e**(-3*x)/(x*(1-mp.e**(-2*x))),[0,1,mp.inf])
N=10**7
sieve=np.ones(N+1,bool); sieve[:2]=False
for p in range(2,int(N**0.5)+1):
    if sieve[p]: sieve[p*p::p]=False
lam=np.zeros(N+1)
for p in np.nonzero(sieve)[0]:
    lp=np.log(p); pk=p
    while pk<=N: lam[pk]=lp; pk*=p
n=np.nonzero(lam)[0].astype(float); c=lam[n.astype(int)]/(n*np.log(n)); ln=np.log(n)
out={'b1':float(b1),'levy':[],'thorin':[]}
for s in [0.5,1.0,2.0]:
    ex=float(mp.log(xi(1+mp.mpf(s))/mp.mpf(0.5)))
    A=float(arch(mp.mpf(s)))
    for X in [1e4,1e5,1e6,1e7]:
        m=n<=X; L=np.log(X)
        pr=np.sum(c[m]*(np.exp(-s*ln[m])-1+s*ln[m]))
        comp=float(mp.quad(lambda x:(mp.e**(-s*x)-1+s*x)/x,[0,1,L]))
        tot=float(b1)*s+A+pr-comp
        out['levy'].append(dict(s=s,X=X,arch=A,prime=float(pr),comp=comp,total=tot,exact=ex,err=tot-ex))
        print('s %.1f X %.0e drift %.5f arch %.6f prime %.5f minus %.5f total %.8f exact %.8f err %.1e'%(s,X,float(b1)*s,A,pr,comp,tot,ex,tot-ex))
# Thorin conversion: pi rho1 = b1 + A(theta) + P(theta), A as integral vs digamma form
for th in [5.0,14.134725,30.0]:
    Aint=float(mp.quad(lambda x:2*mp.sin(th*x/2)**2*mp.e**(-3*x)/(1-mp.e**(-2*x)),[0,1,mp.inf],maxdegree=10))
    Adig=float((mp.re(mp.digamma(mp.mpc(1.5,th/2)))-mp.digamma(mp.mpf(1.5)))/2)
    out['thorin'].append(dict(theta=th,Aint=Aint,Adig=Adig))
    print('theta',th,'A integral %.10f  A digamma %.10f'%(Aint,Adig))
# negative part threshold of x*Pi_1 density: e^{-3x}/(1-e^{-2x}) = 1
x0=float(mp.findroot(lambda x: mp.e**(-3*x)-(1-mp.e**(-2*x)),0.3)); print('ac part negative for x >',x0)
out['x0']=x0
json.dump(out,open('levy.json','w'))
