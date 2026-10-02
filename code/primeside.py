import numpy as np, mpmath as mp, json
from scipy import integrate
from polya import xi_pair
mp.mp.dps=30
PI=np.pi
NP = 2_000_000
# von Mangoldt via sieve
lam = np.zeros(NP+1)
isp = np.ones(NP+1, bool); isp[:2]=False
for p in range(2, int(NP**0.5)+1):
    if isp[p]: isp[p*p::p]=False
for p in np.nonzero(isp)[0]:
    lp = np.log(p); q = p
    while q <= NP:
        lam[q] = lp; q *= p
n = np.nonzero(lam)[0]; L = lam[n]; logn = np.log(n)
def dlogxi(s):
    return 1/s + 1/(s-1) - mp.log(mp.pi)/2 + mp.digamma(s/2)/2 + mp.zeta(s,derivative=1)/mp.zeta(s)
out=[]
for alpha in [3.0, 2.0, 1.5]:
    c = float(dlogxi(mp.mpf(alpha)))
    for b in [2.0, 5.0, 14.134725, 20.0, 35.0, 50.0]:
        gam,_ = integrate.quad(lambda x: (1-np.cos(b*x))*np.exp(-alpha*x)/np.expm1(2*x), 0, 60, limit=2000, epsabs=1e-15)
        pole = -(1/(alpha-1) - (alpha-1)/((alpha-1)**2+b**2))
        primes = np.sum(L*np.exp(-alpha*logn)*(1-np.cos(b*logn)))
        bracket = c + gam + pole + primes
        exact = float(mp.re(dlogxi(mp.mpf(alpha)+1j*mp.mpf(b))))
        F,F1 = xi_pair(alpha-0.5, np.array([b]))
        pol = (F1[0]/F[0]).real
        out.append(dict(alpha=alpha,b=b,c=c,gamma=gam,pole=pole,primes=float(primes),bracket=float(bracket),polya=float(pol),exact=exact))
        print(alpha,b,bracket,pol,exact, bracket-exact, pol-exact, flush=True)
json.dump(out, open('results_prime.json','w'), indent=1)
