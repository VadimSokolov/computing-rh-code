import numpy as np, mpmath as mp, json
mp.mp.dps=25
xi=lambda s: mp.mpf(0.5) if s==1 else s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
dl=lambda s: mp.diff(lambda z: mp.log(xi(z)), s)
b1=float(1+mp.euler/2-mp.log(4*mp.pi)/2)
out={}
# Thorin density rho_1 on a grid (exact), Weyl tail beyond 200
th=np.concatenate([np.linspace(0,60,2401),np.linspace(60.05,200,1400)])
r=np.array([float(mp.re(dl(mp.mpc(1,t)))/mp.pi) if t>0 else b1/np.pi for t in th])
weyl=lambda t: mp.log(t/(2*mp.pi))/(2*mp.pi)
rows=[]
for k in range(1,7):
    w=(2*k-1)**2
    I=np.trapezoid(np.log1p(w/th[1:]**2)*r[1:],th[1:])+r[0]*float(mp.quad(lambda x: mp.log(1+w/x**2),[0,th[1]]))+float(mp.quad(lambda x: mp.log(1+w/x**2)*weyl(x),[200,mp.inf]))
    Ck=2*k*(2*k-1)*float(mp.factorial(k-1))/float(mp.pi)**k
    thor=Ck*np.exp(-I); exact=float(1/mp.zeta(2*k))
    stable=np.exp(-b1*(2*k-1)); rest=np.exp(-I)/stable
    rows.append(dict(k=k,thorin=thor,exact=exact,stable=float(stable),rest=float(rest),C=Ck))
    print('k',k,'1/zeta(2k)',exact,' Thorin form',thor,' stable factor',stable,' remaining GGC factor',rest)
out['coef']=rows
# prime side Levy form check: log xi(2k)/xi(1) = b1 s + int (e^{-sx}-1+sx) Pi_1(dx), s=2k-1
N=10**7
mu=None
sieve=np.ones(N+1,bool); sieve[:2]=False
for p in range(2,int(N**0.5)+1):
    if sieve[p]: sieve[p*p::p]=False
lam=np.zeros(N+1)
for p in np.nonzero(sieve)[0]:
    lp=np.log(p); pk=p
    while pk<=N: lam[pk]=lp; pk*=p
n=np.nonzero(lam)[0].astype(float); c=lam[n.astype(int)]/(n*np.log(n)); ln=np.log(n)
lev=[]
for k in [1,2,3]:
    s=2*k-1
    A=float(mp.quad(lambda x:(mp.e**(-s*x)-1+s*x)*mp.e**(-3*x)/(x*(1-mp.e**(-2*x))),[0,1,mp.inf]))
    pr=np.sum(c*(np.exp(-s*ln)-1+s*ln)); L=np.log(N)
    comp=float(mp.quad(lambda x:(mp.e**(-s*x)-1+s*x)/x,[0,1,L]))
    tot=b1*s+A+pr-comp; ex=float(mp.log(xi(mp.mpf(2*k))/mp.mpf(0.5)))
    lev.append(dict(k=k,levy=float(tot),exact=ex)); print('k',k,'Levy form (primes to 1e7)',tot,' exact log xi(2k)/xi(1)',ex)
out['levy']=lev
json.dump(out,open('ggc_coef.json','w'))
