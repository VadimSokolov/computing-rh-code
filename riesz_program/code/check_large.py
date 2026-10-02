import numpy as np, mpmath as mp, json, time, sys
from riesz_large import riesz
mp.mp.dps=20
g=np.load('g_1_400.npy')[:60]
coef=[complex(mp.gamma(1-mp.mpc(0.5,t)/2)/(2*mp.zeta(mp.mpc(0.5,t),derivative=1))) for t in g]
triv=[(k,float(mp.gamma(1+k)/(2*mp.zeta(-2*k,derivative=1)))) for k in range(1,4)]
pred=lambda x: sum(2*(c*x**(0.5j*t)).real for c,t in zip(coef,g))+sum(r*x**(-k) for k,r in triv)/x**0.25
rows=[]
for x in [1e8,1e10,1e11,1e12]:
    t0=time.time(); R=riesz(x)
    rows.append(dict(x=x,R=R,ratio=R/x**0.25,pred=pred(x),sec=time.time()-t0))
    print('x %.0e  R/x^1/4 %.6e  explicit %.6e  diff %.1e  (%.1fs)'%(x,R/x**0.25,pred(x),R/x**0.25-pred(x),time.time()-t0),flush=True)
json.dump(rows,open('check_large.json','w'))
