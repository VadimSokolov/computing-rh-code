# Diagnostic for the reconstruction of basepoint_one/wiener/riesz_explicit.json (not an authors' script).
# Question: what does maxdev = 3.11e-8, attained at the end point x = 1e10 of [1e3, 1e10], measure? ch/strip.tex (Section wr:sec:riesz)
# says that beyond x = 1e10 double precision cancellation limits the accuracy of the Mobius series. This script recomputes R(x) on the
# grid of salem_riesz.json with the formula and the sieve of basepoint_one/wiener/salem_riesz.py, once with the archived N = 1e7 (this
# must reproduce salem_riesz.json) and once with N = 1e8, and compares both with the explicit formula of riesz_explicit.py (60 zeros,
# three trivial terms). It also prints the neglected fourth trivial term, the Mertens values M(1e7), M(1e8) and sum_{1e7<n<=1e8} mu(n)/n^4,
# which gives the leading part -x^2 sum_{n>N} mu(n) n^-4 of the terms dropped by the N = 1e7 computation.
# Inputs beside it: salem_riesz.json, g_1_400.npy, check_large.json (riesz_program/data, riesz_large.py with K = 200, extended precision).
import json, time, numpy as np, mpmath as mp
mp.mp.dps=20
S=json.load(open('salem_riesz.json')); CL=json.load(open('check_large.json'))
g=np.load('g_1_400.npy')[:60]
coef=[complex(mp.gamma(1-mp.mpc(0.5,t)/2)/(2*mp.zeta(mp.mpc(0.5,t),derivative=1))) for t in g]
triv=[(k,float(mp.gamma(1+k)/(2*mp.zeta(-2*k,derivative=1)))) for k in range(1,4)]
pred=lambda x: sum(2*(c*x**(0.5j*t)).real for c,t in zip(coef,g))+sum(r*x**(-k) for k,r in triv)/x**0.25
triv4=float(mp.gamma(5)/(2*mp.zeta(-8,derivative=1)))
print('fourth trivial coefficient 4!/(2 zeta\'(-8)) =',triv4)
xs=np.array(S['riesz']['x']); Rarch=np.array(S['riesz']['R'])
N1=10**7; N2=10**8; t0=time.time()
# Mobius sieve exactly as in salem_riesz.py, to N2
mu=np.ones(N2+1,dtype=np.int8); isp=np.ones(N2+1,bool); isp[:2]=False
for p in range(2,N2+1):
    if isp[p]:
        if p*p<=N2: isp[p*p::p]=False
        mu[p::p]*=-1
        if p*p<=N2: mu[p*p::p*p]=0
print('mobius sieve to 1e8 %.1f s'%(time.time()-t0),flush=True)
del isp
M1=int(mu[1:N1+1].astype(np.int64).sum()); M2=int(mu[1:].astype(np.int64).sum())
n=np.arange(1,N2+1,dtype=np.float64); m=mu[1:].astype(np.float64); del mu
w=m/n**2
S4=float(np.sum(m[N1:]/n[N1:]**4))
print('Mertens M(1e7) = %d, M(1e8) = %d; sum_{1e7<n<=1e8} mu(n)/n^4 = %.6e; -M(1e7)/1e28 = %.3e'%(M1,M2,S4,-M1/1e28),flush=True)
out=[]
print('%-10s %-17s %-17s %-17s %-10s %-10s %-10s %-10s %-9s'%('x','R/x^1/4 N=1e7','R/x^1/4 N=1e8','explicit','dev N=1e7','dev N=1e8','x^1.75 S4','4th triv','same'))
for x,ra in zip(xs,Rarch):
    R1=x*(np.sum(w[:N1]*np.expm1(-x/n[:N1]**2))+6/np.pi**2)
    R2=x*(np.sum(w*np.expm1(-x/n**2))+6/np.pi**2)
    q1=R1/x**0.25; q2=R2/x**0.25; p=pred(x); t4=triv4*x**(-4)/x**0.25; tr=x**1.75*S4
    out.append(dict(x=float(x),R1=float(R1),R2=float(R2),pred=float(p),dev1=float(q1-p),dev2=float(q2-p),x175S4=tr,triv4=t4,same_as_archive=bool(R1==ra)))
    print('%-10.3e % .10e % .10e % .10e % .2e % .2e % .2e % .2e %s'%(x,q1,q2,p,q1-p,q2-p,tr,t4,R1==ra),flush=True)
for r in CL: print('check_large.json (riesz_large.py, K = 200): x %.0e  R/x^1/4 %.10e  explicit %.10e  diff %.2e'%(r['x'],r['ratio'],r['pred'],r['ratio']-r['pred']))
mk=lambda lo,hi,key: max(abs(r[key]) for r in out if lo<=r['x']<=hi)
for lo,hi in [(1e3,1e10),(1e3,1e9),(1e4,1e10),(1e3,1e11),(1e3,1e12)]:
    print('max |dev| on [%.0e,%.0e]: N=1e7 %.3e   N=1e8 %.3e'%(lo,hi,mk(lo,hi,'dev1'),mk(lo,hi,'dev2')))
json.dump(dict(M1e7=M1,M1e8=M2,S4_1e7_1e8=S4,triv4=triv4,rows=out),open('truncation_check.json','w'))
