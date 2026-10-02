"""Riesz function R(x) = x * sum_n mu(n) n^-2 exp(-x/n^2) for large x.
Segmented Mobius sieve, extended precision (numpy longdouble) and exact handling of the tail:
  R(x)/x = sum_{n0<=n<=N} mu(n) n^-2 e^{-x/n^2} + [6/pi^2 - sum_{n<=N} mu(n)/n^2] + E,
  E = -x * sum_{n>N} mu(n)/n^4 + ... is the next term of the tail expansion. It is dropped: computing it as a
  difference of O(1) sums would amplify rounding by the factor x. With square root cancellation in the Mobius
  tail, |E| is about x N^(-3.5); K = 200 keeps it below 1e-4 of R(x)/x for x up to about 1e14.
with n0 = sqrt(x)/12 (earlier terms are below e^-144) and N = K*sqrt(x), default K = 200.
Twisted version (Mellin transform Gamma(1-s+iT)/zeta(2s)), sensitive to zeros near height 2T:
  R_T(x) = x * sum_n mu(n) n^-2 (x/n^2)^{iT} e^{-x/n^2}, tail via x^{iT}/zeta(2+2iT).
Usage: python3 riesz_large.py 1e12 [K] [T]
"""
import sys, time, numpy as np, mpmath as mp
LD=np.longdouble
def mobius_segments(N, block=10**7):
    # simple segmented Mobius sieve using primes up to sqrt(N)
    r=int(N**0.5)+1
    isp=np.ones(r+1,bool); isp[:2]=False
    for p in range(2,int(r**0.5)+1):
        if isp[p]: isp[p*p::p]=False
    primes=np.nonzero(isp)[0]
    for lo in range(1,N+1,block):
        hi=min(N,lo+block-1); L=hi-lo+1
        mu=np.ones(L,dtype=np.int8); rem=np.arange(lo,hi+1,dtype=np.int64)
        for p in primes:
            if p*p>hi: 
                pass
            s=((lo+p-1)//p)*p
            if s>hi: continue
            idx=np.arange(s-lo,L,p); mu[idx]*=-1; rem[idx]//=p
            p2=p*p; s2=((lo+p2-1)//p2)*p2
            if s2<=hi: mu[np.arange(s2-lo,L,p2)]=0
        mu[rem>1]*=-1          # one prime factor above sqrt(N) remains
        yield lo, mu
def riesz(x, K=200, T=0.0):
    x=float(x); N=int(K*np.sqrt(x)); n0=max(1,int(np.sqrt(x)/12))
    mp.mp.dps=40
    v=mp.mpf(6)/mp.pi**2; hi=float(v); lo=float(v-hi); six=LD(hi)+LD(lo)   # 6/pi^2 to extended precision
    if T!=0.0: return riesz_twisted(x,K,T)
    s_main=LD(0); s2=LD(0); xl=LD(x)
    for lo,mu in mobius_segments(N):
        n=np.arange(lo,lo+len(mu),dtype=LD); m=mu.astype(LD)
        s2+=np.sum(m/n**2)
        sel=n>=n0
        if np.any(sel):
            nn=n[sel]; s_main+=np.sum(m[sel]/nn**2*np.exp(-xl/nn**2))
    tail=six-s2
    return float(xl*(s_main+tail))
def riesz_twisted(x, K=200, T=50.0):
    x=float(x); N=int(K*np.sqrt(x)); n0=max(1,int(np.sqrt(x)/12))
    mp.mp.dps=40
    iz=complex(1/mp.zeta(2+2j*mp.mpf(T)))
    CL=np.clongdouble; xl=LD(x); lx=np.log(xl)
    s_main=CL(0); s2=CL(0)
    for lo,mu in mobius_segments(N):
        n=np.arange(lo,lo+len(mu),dtype=LD); m=mu.astype(LD); ln=np.log(n)
        s2+=np.sum(m/n**2*np.exp(-2j*LD(T)*ln))
        sel=n>=n0
        if np.any(sel):
            nn=n[sel]; s_main+=np.sum(m[sel]/nn**2*np.exp(1j*LD(T)*(lx-2*ln[sel]))*np.exp(-xl/nn**2))
    tail=np.exp(1j*LD(T)*lx)*(CL(iz)-s2)
    return complex(xl*(s_main+tail))
if __name__=='__main__':
    x=float(sys.argv[1]) if len(sys.argv)>1 else 1e10
    K=float(sys.argv[2]) if len(sys.argv)>2 else 200
    T=float(sys.argv[3]) if len(sys.argv)>3 else 0.0
    t0=time.time(); R=riesz(x,K,T)
    if T==0.0: print('x=%.3e  R(x)=%.10e  R/x^(1/4)=%.6e  (%.1fs)'%(x,R,R/x**0.25,time.time()-t0))
    else: print('x=%.3e  T=%g  R_T(x)=%s  |R_T|/x^(1/4)=%.6e  (%.1fs)'%(x,T,R,abs(R)/x**0.25,time.time()-t0))
