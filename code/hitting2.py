import json, numpy as np, mpmath as mp
from scipy.stats import norm
from scipy.optimize import brentq
from scipy.interpolate import CubicSpline
from scipy.integrate import solve_ivp
D=json.load(open('hazard_spec.json')); rows=D['rows']
t=np.array([r[0] for r in rows]); S=np.array([float(mp.mpf(r[2])) for r in rows]); h=np.array([float(mp.mpf(r[3])) for r in rows])
hm=[mp.mpf(r[3]) for r in rows]; g1sq=mp.mpf(D['g'][0])**2; g2sq=mp.mpf(D['g'][1])**2
dev=np.array([float(abs(x-g1sq)) for x in hm])
sel=(t>=0.05)&(t<=0.12); slope=np.polyfit(t[sel],np.log(dev[sel]),1)[0]
print('decay slope',slope,'predicted',float(-(g2sq-g1sq)))
lh=CubicSpline(t,np.log(h)); lS=CubicSpline(t,np.log(S))
def Htan(alpha,beta,tau):
    x=(alpha+beta*tau)/np.sqrt(tau); y=(beta*tau-alpha)/np.sqrt(tau)
    lg=np.log(alpha)-0.5*np.log(2*np.pi*tau**3)-x*x/2
    surv=norm.cdf(x)-np.exp(-2*alpha*beta+norm.logcdf(y))
    return np.exp(lg)/surv
def beta_of(b,tau):
    target=np.log(np.exp(lh(tau)))
    F=lambda be: np.log(Htan(b-tau*be,be,tau))-target
    hi=b/tau*(1-1e-12)
    grid=np.concatenate([-np.logspace(5,-3,300),np.linspace(0,hi,600)])
    v=np.array([F(x) for x in grid])
    for i in range(len(grid)-1):
        if np.isfinite(v[i]) and np.isfinite(v[i+1]) and v[i]*v[i+1]<0: return brentq(F,grid[i],grid[i+1],xtol=1e-12)
    return np.nan
t0=0.011; F0=1-np.exp(lS(t0)); b0=np.sqrt(t0)*norm.isf(F0/2); print('t0',t0,'F0',F0,'b0',b0)
for b in [b0]:
    print('beta at t0', beta_of(b,t0))
sol=solve_ivp(lambda tau,b:[beta_of(b[0],tau)],[t0,0.3],[b0],max_step=0.0005,rtol=1e-9,atol=1e-11,dense_output=True)
print(sol.status, sol.message)
tb=np.linspace(t0,sol.t[-1],400); bb=sol.sol(tb)[0]; beta=np.array([beta_of(x,y) for x,y in zip(bb,tb)])
print('b at', [(round(a,3),round(c,4),round(d,3)) for a,c,d in zip(tb[::50],bb[::50],beta[::50])])
print('asymptotic slope target', float(mp.sqrt(2*g1sq)))
json.dump(dict(slope=slope,pred=float(-(g2sq-g1sq)),tb=list(tb),b=list(bb),beta=list(beta),t0=t0,b0=b0),open('hitting2.json','w'))
