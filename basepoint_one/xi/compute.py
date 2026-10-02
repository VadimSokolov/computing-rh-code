import numpy as np, mpmath as mp, json, time
mp.mp.dps=20
g=np.concatenate([np.load(f) for f in ['g_1_400.npy','g_401_1000.npy','g_1001_1700.npy']])
E=float(mp.euler)
def rho1(th):
    s=mp.mpc(1,th)
    z=mp.euler if th==0 else mp.re(mp.zeta(s,derivative=1)/mp.zeta(s))
    return float((mp.re(1/s)-mp.log(mp.pi)/2+mp.re(mp.digamma(s/2))/2+z)/mp.pi)
out={}
# 1. density on a grid and cumulative mass vs N(T)
th=np.linspace(0,100,2001); r=np.array([rho1(t) for t in th])
M=np.concatenate([[0],np.cumsum((r[1:]+r[:-1])/2*np.diff(th))])
out['grid']=dict(th=th.tolist(),rho=r.tolist(),M=M.tolist())
print('min rho on [0,100]',r.min(),'at',th[r.argmin()],' M(100)=',M[-1],' N(100)=',int(np.sum(g<100)))
# 2. sine squared lemma at alpha=1 : primes minus Lebesgue, Mertens constant
N=10**7
lam=np.zeros(N+1)
sieve=np.ones(N+1,bool); sieve[:2]=False
for p in range(2,int(N**0.5)+1):
    if sieve[p]: sieve[p*p::p]=False
P=np.nonzero(sieve)[0]
for p in P:
    lp=np.log(p); pk=p
    while pk<=N: lam[pk]=lp; pk*=p
n=np.nonzero(lam)[0]; w=lam[n]/n; ln=np.log(n)
rows=[]
for t in [5.0,14.134725,30.0]:
    exact=np.pi*rho1(t)
    arch=1/(1+t*t)-float(mp.log(mp.pi))/2+float(mp.re(mp.digamma(mp.mpc(0.5,t/2))))/2+E
    for X in [1e3,1e4,1e5,1e6,1e7]:
        m=n<=X; L=np.log(X)
        pr=np.sum(w[m]*2*np.sin(t*ln[m]/2)**2)
        pole=L-np.sin(t*L)/t
        approx=arch+pr-pole
        rows.append(dict(theta=t,X=X,prime=float(pr),pole=float(pole),approx=float(approx),exact=exact,err=float(approx-exact)))
        print('theta %.4f X %.0e prime %.4f pole %.4f sum %.6f exact %.6f err %.2e'%(t,X,pr,pole,approx,exact,approx-exact))
out['sin2']=rows
# 3. no-go: growth of |xi(1)/xi(1+it)|
xi=lambda s: s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
x1=mp.mpf(0.5); ng=[]
for t in [10,20,40,80]:
    v=abs(x1/xi(mp.mpc(1,t))); ng.append((t,float(v))); print('t',t,'|phi(-t^2)|',mp.nstr(v,5))
out['nogo']=ng
json.dump(out,open('alpha1.json','w'))
