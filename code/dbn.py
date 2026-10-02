import json, numpy as np, mpmath as mp
from flint import arb, acb, acb_poly, ctx
from polya_flint import setup
P1=json.load(open('part1.json')); g=np.array([float(x) for x in P1['roots']])
gaps=np.diff(g); k=int(np.argmin(gaps)); print('closest pair',k+1,k+2,g[k],g[k+1],gaps[k])
S0=setup(420,200,3); h=S0['h']
Phi=[c.real/h for c in S0["PA"].coeffs()]; Phi[0]=Phi[0]*2   # recover Phi_j
def Xi_t(t, theta):
    # Xi_t(theta) = int Phi(v) exp(t v^2/4) cos(theta v) dv, trapezoid
    tt=arb(t)/4; th=arb(theta); s=Phi[0]
    for j in range(1,len(Phi)):
        u=j*h; s+=2*Phi[j]*(tt*u*u).exp()*(th*u).cos()
    return float(s*h)
def Xi_t_jet(t, theta):
    # Xi_t and its derivatives d/dtheta, d2/dtheta2 and d2/dtheta dt, in ball arithmetic; d/dt = -(1/4) d2/dtheta2
    tt=arb(t)/4; th=arb(theta); g0=Phi[0]; g1=arb(0); g2=arb(0); g3=arb(0)
    for j in range(1,len(Phi)):
        u=j*h; w=2*Phi[j]*(tt*u*u).exp(); sn,cs=(th*u).sin_cos()
        g0+=w*cs; g1-=w*u*sn; g2-=w*u*u*cs; g3-=w*u*u*u*sn/4
    return g0*h, g1*h, g2*h, g3*h
def refine(t, a, b):
    # Newton's method for the zero of Xi_t in the grid cell (a,b), started at its midpoint
    th=arb((a+b)/2)
    for it in range(50):
        g0,g1,_,_=Xi_t_jet(t,th); d=g0/g1; th=arb((th-d).mid())
        if abs(d)<arb(10)**-40: break
    assert arb(a)<th and th<arb(b), (t,a,b)
    return th
out={'pair':[k+1,k+2,g[k],g[k+1]]}
lo,hi=g[k]-0.3,g[k+1]+0.3
def nsign(t):
    th=np.linspace(lo,hi,121); v=np.array([Xi_t(t,x) for x in th]); return int(np.sum(np.sign(v[:-1])!=np.sign(v[1:]))), th, v
traj=[]
for t in [0.2,0.1,0.0,-0.05,-0.1,-0.15,-0.2,-0.3]:
    n,th,v=nsign(t); cells=[(th[i],th[i+1]) for i in range(len(th)-1) if np.sign(v[i])!=np.sign(v[i+1])]
    mids=[float((a+b)/2) for a,b in cells]; zs=[float(refine(t,a,b)) for a,b in cells]
    traj.append(dict(t=t,n=n,z=zs,mid=mids)); print(t,n,'refined',[round(z,4) for z in zs],'grid midpoints',[round(z,3) for z in mids],flush=True)
# bisection on the collision time
a,b=0.0,-0.3
for it in range(18):
    m=(a+b)/2; n,_,_=nsign(m)
    if n>=2: a=m
    else: b=m
out['collision']=(a+b)/2; out['traj']=traj
print('collision time (our normalisation)',(a+b)/2, 'de Bruijn Newman normalisation t/4 ?')
# Newton's method for the double zero, Xi_t = d/dtheta Xi_t = 0 in the unknowns theta and t, started from the bisection
# value and the grid point between the two sign changes at the upper end of the final bracket
n,th,v=nsign(a); sc=[i for i in range(len(th)-1) if np.sign(v[i])!=np.sign(v[i+1])]
thc=arb((th[sc[0]]+th[sc[-1]+1])/2); tc=arb((a+b)/2)
for it in range(50):
    g0,g1,g2,g3=Xi_t_jet(tc,thc); gt=-g2/4; det=g1*g3-gt*g2
    dth=(gt*g1-g0*g3)/det; dt=(g0*g2-g1*g1)/det
    thc=arb((thc+dth).mid()); tc=arb((tc+dt).mid())
    if abs(dth)<arb(10)**-40 and abs(dt)<arb(10)**-40: break
out['collision_newton']=dict(t=float(tc),theta=float(thc),iterations=it+1)
print('collision time by Newton',tc.str(20,radius=False),'at theta',thc.str(20,radius=False),flush=True)
json.dump(out,open('dbn.json','w'))
