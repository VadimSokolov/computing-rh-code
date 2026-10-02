import numpy as np, json, time
from flint import acb, arb, ctx
import tilted
S=tilted.setup(0.01,256,0.0006,3.9)
J=S['J']; h=S['h']; eta=S['eta']
nodes=[acb((k-J)*h, eta) for k in range(2*J+1)]
P0=list(S['P'].coeffs()) if hasattr(S['P'],'coeffs') else S['P']
def Z(t,theta):
    # H_t up to a positive factor: int Phi(v) exp(t v^2/4) exp(i theta v) dv on the tilted line
    s=acb(0,arb(str(theta))); tt=arb(str(t))/4; A=acb(0)
    for k,v in enumerate(nodes):
        A+=P0[k]*(tt*v*v+s*v).exp()
    return float(A.real*(arb.pi()*arb(str(theta))/4).exp())
V2=[v*v for v in nodes]; IV=[acb(0,1)*v for v in nodes]
def jet(t,theta):
    # the same sum without the positive factor, and its derivatives d/dtheta, d2/dtheta2 and d2/dtheta dt (real parts,
    # in ball arithmetic); d/dt = -(1/4) d2/dtheta2
    s=acb(0,arb(theta)); tt=arb(t)/4; g0=acb(0); g1=acb(0); g2=acb(0); g3=acb(0)
    for k in range(len(nodes)):
        e=P0[k]*(tt*V2[k]+s*nodes[k]).exp(); ie=IV[k]*e
        g0+=e; g1+=ie; g2+=IV[k]*ie; g3+=V2[k]*ie/4
    return g0.real, g1.real, g2.real, g3.real
grid=np.round(np.arange(7005.03,7005.1301,0.0025),4)
def nsign(t):
    v=[Z(t,x) for x in grid]; return sum(1 for i in range(len(v)-1) if v[i]*v[i+1]<0), v
t0=time.time(); n0,v0=nsign(0.0); print('t=0 sign changes',n0,time.time()-t0,flush=True)
a,b=0.0,-0.004
nb,_=nsign(b); print('t=-0.004',nb,flush=True)
va=v0
for it in range(12):
    m=(a+b)/2; n,v=nsign(m)
    if n>=2: a,va=m,v
    else: b=m
    print(it,m,n,flush=True)
print('collision by bisection',(a+b)/2,time.time()-t0,flush=True)
# Newton's method for the double zero, H_t(2 theta) = d/dtheta H_t(2 theta) = 0 in the unknowns theta and t, started from
# the bisection value and the grid point between the two sign changes at the upper end of the final bracket
sc=[i for i in range(len(grid)-1) if va[i]*va[i+1]<0]
thc=arb(str(round((grid[sc[0]]+grid[sc[-1]+1])/2,5))); tc=arb(str((a+b)/2))
for it in range(30):
    g0,g1,g2,g3=jet(tc,thc); gt=-g2/4; det=g1*g3-gt*g2
    dth=(gt*g1-g0*g3)/det; dt=(g0*g2-g1*g1)/det
    thc=arb((thc+dth).mid()); tc=arb((tc+dt).mid())
    print('newton',it,tc.str(20,radius=False),thc.str(20,radius=False),flush=True)
    if abs(dth)<arb(10)**-30 and abs(dt)<arb(10)**-30: break
json.dump(dict(collision=(a+b)/2,bracket=[a,b],collision_newton=dict(t=float(tc),theta=float(thc),iterations=it+1)),open('lehmer_dbn.json','w'))
print('collision',(a+b)/2,'by Newton',tc.str(20,radius=False),'at theta',thc.str(20,radius=False),time.time()-t0)
