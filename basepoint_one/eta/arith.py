import numpy as np, mpmath as mp, json
mp.mp.dps=20
chi=[0,1,0,0,0,-1,0,-1,0,0,0,1]
Lam=lambda s: (12/mp.pi)**(s/2)*mp.gamma(s/2)*mp.dirichlet(s,chi)
dlog=lambda s: mp.log(12/mp.pi)/2+mp.digamma(s/2)/2+mp.dirichlet(s,chi,derivative=1)/mp.dirichlet(s,chi)
g=np.load('g12.npy'); B=json.load(open('base.json')); b1=B['b']
N=10**7
sieve=np.ones(N+1,bool); sieve[:2]=False
for p in range(2,int(N**0.5)+1):
    if sieve[p]: sieve[p*p::p]=False
lam=np.zeros(N+1)
for p in np.nonzero(sieve)[0]:
    lp=np.log(p); pk=p; k=1
    cp=chi[p%12]
    while pk<=N: lam[pk]=lp*cp**k; pk*=p; k+=1
n=np.nonzero(lam)[0].astype(float); L_=lam[n.astype(int)]; ln=np.log(n)
out={'levy':[],'sin2':[],'err':[]}
# 1. Levy form at alpha=1
arch=lambda s: mp.quad(lambda x:(mp.e**(-s*x)-1+s*x)*mp.e**(-x)/(x*(1-mp.e**(-2*x))),[0,1,mp.inf])
for s in [0.5,1.0,2.0]:
    ex=float(mp.log(Lam(1+mp.mpf(s))/Lam(mp.mpf(1)))); A=float(arch(mp.mpf(s)))
    for X in [1e5,1e7]:
        m=n<=X; pr=float(np.sum(L_[m]/(n[m]*ln[m])*(np.exp(-s*ln[m])-1+s*ln[m])))
        tot=b1*s+A+pr
        out['levy'].append(dict(s=s,X=X,A=A,prime=pr,total=tot,exact=ex))
        print('Levy s %.1f X %.0e drift %.5f gamma %.6f primes %.6f total %.7f exact %.7f'%(s,X,b1*s,A,pr,tot,ex))
# 2. sine squared at alpha and error formula
rho=np.concatenate([0.5+1j*g,0.5-1j*g])
def E(X,al,th):
    a=rho-al; t=-X**a/a+0.5*X**(a+1j*th)/(a+1j*th)+0.5*X**(a-1j*th)/(a-1j*th); return float(np.real(np.sum(t)))
for al in [1.0,0.9,0.75,0.6,0.5]:
    ba=float(mp.re(dlog(mp.mpf(al)))) if al!=0.5 else 0.0
    for th in [3.8046276,10.0]:
        ex=float(mp.re(dlog(mp.mpc(al,th)))) if al!=0.5 else float('nan')
        A=float((mp.re(mp.digamma(mp.mpc(al,th)/2))-mp.digamma(mp.mpf(al)/2))/2)
        for X in [1e3,1e5,1e7]:
            m=n<=X; pr=float(np.sum(L_[m]*n[m]**(-al)*2*np.sin(th*ln[m]/2)**2))
            tot=ba+A+pr; pred=E(X,al,th) if al>0.5 else float('nan')
            out['sin2'].append(dict(alpha=al,theta=th,X=X,b=ba,A=A,prime=pr,total=tot,exact=ex,err=tot-ex,pred=pred,corr=tot-pred-ex))
            print('alpha %.2f th %.4f X %.0e drift %.5f A %.5f primes %9.4f total %9.5f exact %9.5f err %9.2e pred %9.2e corrected %9.2e'%(al,th,X,ba,A,pr,tot,ex,tot-ex,pred,tot-pred-ex),flush=True)
json.dump(out,open('arith.json','w'))
