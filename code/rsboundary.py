import json, numpy as np, mpmath as mp
from scipy.stats import norm
from scipy.optimize import brentq
from scipy.interpolate import CubicSpline
D=json.load(open('hazard_spec.json')); rows=D['rows']
t=np.array([r[0] for r in rows]); S=np.array([float(mp.mpf(r[2])) for r in rows]); h=np.array([float(mp.mpf(r[3])) for r in rows])
lh=CubicSpline(t,np.log(h)); lS=CubicSpline(t,np.log(S))
def lHtan(alpha,beta,tau):
    if alpha<=0: return np.nan
    x=(alpha+beta*tau)/np.sqrt(tau); y=(beta*tau-alpha)/np.sqrt(tau)
    lg=np.log(alpha)-0.5*np.log(2*np.pi*tau**3)-x*x/2
    ls=np.log(norm.cdf(x)-np.exp(-2*alpha*beta+norm.logcdf(y))) if x<8 else np.log1p(-norm.sf(x)-np.exp(-2*alpha*beta+norm.logcdf(y)))
    return lg-ls
def solve_beta(b,tau,guess):
    F=lambda be: lHtan(b-tau*be,be,tau)-lh(tau)
    d=max(1.0,abs(guess)*0.05)
    for k in range(40):
        lo,hi=guess-d,min(guess+d,b/tau*(1-1e-12))
        flo,fhi=F(lo),F(hi)
        if np.isfinite(flo) and np.isfinite(fhi) and flo*fhi<0: return brentq(F,lo,hi,xtol=1e-12)
        d*=1.6
    return np.nan
t0=0.011; F0=1-np.exp(lS(t0)); b=np.sqrt(t0)*norm.isf(F0/2); beta=-87.547
dt=1e-4; tau=t0; T=[tau]; B=[b]; BE=[solve_beta(b,tau,beta)]; beta=BE[0]
while tau<0.35-1e-12:
    k1=solve_beta(b,tau,beta); k2=solve_beta(b+dt/2*k1,tau+dt/2,k1); k3=solve_beta(b+dt/2*k2,tau+dt/2,k2); k4=solve_beta(b+dt*k3,tau+dt,k3)
    if not np.all(np.isfinite([k1,k2,k3,k4])): print('stop at',tau); break
    b+=dt*(k1+2*k2+2*k3+k4)/6; tau+=dt; beta=k4; T.append(tau); B.append(b); BE.append(solve_beta(b,tau,beta))
T=np.array(T);B=np.array(B);BE=np.array(BE)
for q in [0.011,0.015,0.02,0.03,0.05,0.1,0.2,0.3,0.35]:
    i=np.argmin(abs(T-q)); print(round(T[i],4), round(B[i],5), round(BE[i],4), round(B[i]-T[i]*BE[i],5))
print('target slope', np.sqrt(2)*14.134725141734693)
json.dump(dict(t=list(T),b=list(B),beta=list(BE),t0=t0),open('rsboundary.json','w'))
