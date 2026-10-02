import json, numpy as np, mpmath as mp
from scipy.stats import norm
from scipy.optimize import brentq
from scipy.interpolate import CubicSpline
from scipy.integrate import solve_ivp
H=json.load(open('hazard_spec.json'))
rows=H['rows']; t=np.array([r[0] for r in rows]); S=np.array([float(mp.mpf(r[2])) for r in rows]); h=np.array([float(mp.mpf(r[3])) for r in rows])
hm=[mp.mpf(r[3]) for r in rows]
lh=CubicSpline(t,np.log(h)); lS=CubicSpline(t,np.log(S))
g1sq=float(mp.mpf(H['g'][0])**2); g2sq=float(mp.mpf(H['g'][1])**2)
out={}
out['hazard']=dict(t=list(t), h=list(h))
# quasi stationarity rate
sel=[(r[0],abs(hm[i]-mp.mpf(H['g'][0])**2)) for i,r in enumerate(rows) if 0.05<=r[0]<=0.15]
tt=np.array([a for a,b in sel]); d=np.array([float(mp.log(b)) for a,b in sel]); slope=np.polyfit(tt,d,1)[0]
out['qs_slope']=slope; out['g1sq']=g1sq; out['g2sq']=g2sq; out['qs_diff']=[(a,float(b)) for a,b in sel]
print('slope of log|h-g1^2|', slope, 'predicted', -(g2sq-g1sq))
# Roberts Shortland tangent hazard of the line alpha + beta s at time s=tau
def Htan(alpha,beta,tau):
    x=(alpha+beta*tau)/np.sqrt(tau)
    lg=np.log(alpha)-0.5*np.log(2*np.pi*tau**3)-x*x/2
    # survival = Phi(x) - exp(-2 alpha beta) Phi((beta tau - alpha)/sqrt tau), computed stably
    y=(beta*tau-alpha)/np.sqrt(tau)
    surv=norm.cdf(x)-np.exp(-2*alpha*beta+norm.logcdf(y))
    return np.exp(lg)/surv
def beta_of(b,tau):
    target=np.exp(lh(tau))
    F=lambda be: np.log(Htan(b-tau*be,be,tau))-np.log(target)
    grid=np.concatenate([-np.logspace(4,-3,3000),np.linspace(1e-3,b/tau*(1-1e-9),3000)]); vals=[F(x) for x in grid]
    for i in range(len(grid)-1):
        if np.isfinite(vals[i]) and np.isfinite(vals[i+1]) and vals[i]*vals[i+1]<0:
            return brentq(F,grid[i],grid[i+1])
    raise RuntimeError('no root at %g'%tau)
t0=0.0105
F0=1-np.exp(lS(t0)); b0=np.sqrt(t0)*norm.isf(F0/2)
print('t0',t0,'F(t0)',F0,'b0',b0)
sol=solve_ivp(lambda tau,b:[beta_of(b[0],tau)],[t0,0.3],[b0],max_step=0.001,rtol=1e-8,atol=1e-10,dense_output=True)
tb=np.linspace(t0,0.3,600); bb=sol.sol(tb)[0]; beta=np.array([beta_of(x,y) for x,y in zip(bb,tb)])
out['boundary']=dict(t=list(tb),b=list(bb),beta=list(beta))
print('beta at 0.1,0.2,0.25:', np.interp([0.05,0.1,0.2,0.3],tb,beta), 'target sqrt(2 g1^2)=', np.sqrt(2*g1sq))
json.dump(out,open('hitting.json','w'))
