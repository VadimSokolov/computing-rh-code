import mpmath as mp, numpy as np, time, json
from flint import fmpz_poly
mp.mp.dps=60
NMAX=3000
# eta q-series prod(1-q^n) via pentagonal numbers
e=[0]*(NMAX+1); k=0
while True:
    for kk in ([k,-k] if k else [0]):
        m=kk*(3*kk-1)//2
        if m<=NMAX: e[m]+=(-1)**abs(kk)
    k+=1
    if k*(3*k-1)//2>NMAX: break
E=fmpz_poly(e)
def mullow(a,b,n): c=a*b; return fmpz_poly(list(c.coeffs())[:n+1])
P=fmpz_poly([1])
for i in range(24): P=mullow(P,E,NMAX)
tau=[0]+[int(x) for x in list(P.coeffs())[:NMAX]]   # tau(n)=coeff of q^(n-1) in prod^24
assert tau[1]==1 and tau[2]==-24 and tau[3]==252 and tau[11]==534612
print('tau ok')
def chi12(n): return {1:1,5:-1,7:-1,11:1}.get(n%12,0)
# kernels on the real line (even), Pólya type
def KD(u):  # e^{6u} Delta(i e^u), use symmetry for u<0
    u=abs(u); y=mp.exp(u); q=mp.exp(-2*mp.pi*y)
    s=mp.mpf(0); qn=q
    for n in range(1,80):
        s+=tau[n]*qn; qn*=q
        if abs(qn)<mp.mpf(10)**(-mp.mp.dps-5): break
    return mp.exp(6*u)*s
def K12(u):  # e^{u/4} eta(i e^u)
    u=abs(u); y=mp.exp(u); s=mp.mpf(0)
    for n in range(1,200):
        if chi12(n): 
            t=mp.exp(-mp.pi*n*n*y/12)
            s+=chi12(n)*t
            if t<mp.mpf(10)**(-mp.mp.dps-5): break
    return mp.exp(u/4)*s
def Phi(u):
    u=abs(u); s=mp.mpf(0)
    for n in range(1,60):
        a=mp.pi*n*n*mp.exp(4*u); t=(2*a*a-3*a)*mp.exp(u)*mp.exp(-a)
        s+=t
        if abs(t)<mp.mpf(10)**(-mp.mp.dps-5): break
    return 2*s
def nodes(K,h,U):
    us=[mp.mpf(j)*h for j in range(-int(U/h),int(U/h)+1)]
    return us,[K(u) for u in us]
def Lam(us,Ks,h,z,scale=1):   # integral K(u) e^{z u/scale} du
    return h*mp.fsum(Kv*mp.exp(z*u/scale) for u,Kv in zip(us,Ks))
res={}
setups={'zeta':(Phi,0.5,3.0),'chi12':(K12,2,5.0),'Delta':(KD,1,4.0)}
for name,(K,scale,U) in setups.items():
    h=mp.mpf(1)/80
    t0=time.time(); us,Ks=nodes(K,h,U); tn=time.time()-t0
    # Z on the critical line: real function Lambda(1/2+it)
    T=np.arange(0.5,50,0.05); Z=[]
    t0=time.time()
    for t in T: Z.append(float(mp.re(Lam(us,Ks,h,1j*mp.mpf(t),scale))))
    tp=(time.time()-t0)/len(T)
    Z=np.array(Z); sc=np.where(np.diff(np.sign(Z))!=0)[0]
    zs=[]
    for i in sc:
        f=lambda t: mp.re(Lam(us,Ks,h,1j*t,scale))
        zs.append(float(mp.findroot(f,(mp.mpf(T[i]),mp.mpf(T[i+1])),solver='anderson')))
    dec=[float(mp.log10(abs(Lam(us,Ks,h,1j*mp.mpf(t),scale))+mp.mpf(10)**-300)) for t in [10.3,20.3,30.3,40.3]]
    res[name]=dict(nodes=len(us),node_time=tn,per_point=tp,zeros=zs,log10abs=dec,K0=float(Ks[len(Ks)//2]))
    print(name,'nodes',len(us),'per point %.4fs'%tp,'zeros<50:',len(zs),'first',[round(x,6) for x in zs[:5]],'log10|Lambda| at 10,20,30,40:',[round(x,1) for x in dec],flush=True)
json.dump(res,open('lfun.json','w'))
