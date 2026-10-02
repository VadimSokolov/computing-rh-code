import numpy as np, mpmath as mp, json, time
mp.mp.dps=15
out={}
# arithmetic factor a_k = prod_p (1-1/p)^{k^2} sum_m (Gamma(m+k)/(m! Gamma(k)))^2 p^{-m}
def ak(k,P=200000):
    ps=[p for p in range(2,P) if mp.isprime(p)] if False else None
    import sympy
    s=mp.mpf(0)
    for p in sympy.primerange(2,P):
        p=mp.mpf(p); ser=mp.hyp2f1(k,k,1,1/p)   # sum (k)_m^2/(m!)^2 p^-m
        s+=k*k*mp.log(1-1/p)+mp.log(ser)
    return mp.e**s
A={1:1.0000000000000107,2:0.607927332969108,3:0.04932184233405926}
print('6/pi^2 =',float(6/mp.pi**2))
G=lambda k: float(mp.barnesg(1+k)**2/mp.barnesg(1+2*k))
lead={k:A[k]*G(k) for k in A}; print('leading KS coefficients',lead,' 1/(2 pi^2)=',1/(2*np.pi**2))
out['ak']=A; out['lead']=lead
# moments of zeta on the critical line
t0=time.time()
ts=np.arange(1.0,1000.0,0.05)
z=np.array([abs(complex(mp.zeta(mp.mpc(0.5,t)))) for t in ts])
print('zeta values computed in %.0fs'%(time.time()-t0),flush=True)
mom=[]
dt=0.05
for T in [125,250,500,1000]:
    m=ts<=T
    row=dict(T=T)
    for k in [1,2,3]:
        M=np.sum(z[m]**(2*k))*dt/T
        row['m%d'%k]=float(M); row['lead%d'%k]=float(lead[k]*np.log(T/(2*np.pi))**(k*k))
    row['exact1']=float(np.log(T/(2*np.pi))+2*np.euler_gamma-1)
    mom.append(row); print(T,{k:round(v,4) for k,v in row.items() if k!='T'})
out['moments']=mom
# Selberg central limit theorem: log|zeta(1/2+it)| on [1000,2000]
m=(ts>=500)
L=np.log(z[m][z[m]>0]); v=0.5*np.log(np.log(750))
out['selberg']=dict(mean=float(L.mean()),var=float(L.var()),pred_var=float(v))
print('log|zeta| on [500,1000]: mean',L.mean(),'var',L.var(),' Selberg variance 1/2 loglog t =',v)
# Thorin density as the sigma derivative of log|xi|: log|xi(1/2+it)| = log|xi(1+it)| - pi int_{1/2}^1 rho_a(t) da
xi=lambda s: s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
dl=lambda s: mp.diff(lambda w: mp.log(xi(w)), s)
chk=[]
for t in [10.0,14.134725,20.0,100.0]:
    lhs=float(mp.log(abs(xi(mp.mpc(0.5,t))))) if abs(t-14.134725)>1e-3 else float('-inf')
    integ=float(mp.quad(lambda a: mp.re(dl(mp.mpc(a,t))),[0.5,1]))
    rhs=float(mp.log(abs(xi(mp.mpc(1,t)))))-integ
    mono=all(abs(xi(mp.mpc(a,t)))<abs(xi(mp.mpc(a+0.05,t))) for a in np.arange(0.5,1.0,0.05))
    chk.append(dict(t=t,lhs=lhs,rhs=rhs,monotone=bool(mono))); print('t',t,'log|xi(1/2+it)|',lhs,' via Thorin integral',rhs,' |xi(sigma+it)| increasing on [1/2,1]:',mono)
out['thorin_logxi']=chk
json.dump(out,open('ks.json','w'))
