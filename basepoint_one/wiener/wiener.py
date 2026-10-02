import numpy as np, json, time
from scipy.optimize import linprog
import mpmath as mp
def Phi(u):
    u=np.abs(u); s=np.zeros_like(u)
    for n in range(1,8):
        a=np.pi*n*n*np.exp(2*u)
        s+=2*(2*np.pi**2*n**4*np.exp(4.5*u)-3*np.pi*n*n*np.exp(2.5*u))*np.exp(-a)
    return s
u=np.linspace(-3,3,1201); du=u[1]-u[0]
sig=0.15; g=np.exp(-u**2/(2*sig**2))/np.sqrt(2*np.pi*sig**2)
xi=lambda s: s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
print('check: int Phi = xi(1/2)?',np.sum(Phi(u))*du,float(xi(mp.mpf(0.5))))
g1=14.134725141734693; lb=np.exp(-sig**2*g1**2/2)
print('lower bound at the centre |g^(gamma_1)| =',lb)
out={'lb':lb,'rows':[]}
for al in [0.5,0.6,0.75,1.0]:
    for N in [10,20,30,40,60,80]:
        xs=np.linspace(-1.5,1.5,N)
        A=np.array([np.exp((al-0.5)*(u-x))*Phi(u-x) for x in xs]).T
        sc=np.abs(A).sum(0)*du; A=A/sc
        M=len(u)
        # variables: c (N free), e (M >=0); minimize du*sum e ; constraints: A c - e <= g ; -A c - e <= -g
        cvec=np.concatenate([np.zeros(N),du*np.ones(M)])
        Aub=np.block([[A,-np.eye(M)],[-A,-np.eye(M)]]); bub=np.concatenate([g,-g])
        t0=time.time()
        res=linprog(cvec,A_ub=Aub,b_ub=bub,bounds=[(None,None)]*N+[(0,None)]*M,method='highs')
        if res.x is None: print('LP failed',al,N); continue
        c=res.x[:N]; h=A@c; c=c/sc; err=float(np.sum(np.abs(g-h))*du)
        hhat=abs(np.sum(h*np.exp(1j*g1*u))*du); cn=float(np.sum(np.abs(c)))
        out['rows'].append(dict(alpha=al,N=N,err=err,hhat_g1=float(hhat),coefnorm=cn))
        print('alpha %.2f N %3d L1 error %.5f  |h^(gamma1)| %.4f  sum|c| %.3e  (%.1fs)'%(al,N,err,hhat,cn,time.time()-t0),flush=True)
# |xi(alpha+i gamma_1)| and |xi'(1/2+i gamma_1)|, quoted in Section wr:sec:wiener of the book
# (the archived wiener.json holds them, but earlier versions of this script did not compute them)
mp.mp.dps=20
gm=mp.zetazero(1).imag
out['amp']={'0.5':0.0}  # 1/2+i gamma_1 is a zero
for al in [0.6,0.75,1.0]: out['amp'][str(al)]=float(abs(xi(mp.mpc(al,gm))))
out['dxi']=float(abs(mp.diff(xi,mp.mpc(0.5,gm))))
print('|xi(alpha+i gamma_1)|',out['amp'],"|xi'(1/2+i gamma_1)|",out['dxi'])
json.dump(out,open('wiener.json','w'))
