import numpy as np, mpmath as mp, json
from scipy.special import digamma
from flint import acb
from polya_flint import setup, xi_pair
S=setup(300,200,3)
g=np.concatenate([np.load(f) for f in ['g_1_400.npy','g_401_1000.npy','g_1001_1700.npy']])
def vm(N):
    lam=np.zeros(N+1); isp=np.ones(N+1,bool); isp[:2]=False
    for p in range(2,int(N**0.5)+1):
        if isp[p]: isp[p*p::p]=False
    for p in np.nonzero(isp)[0]:
        lp=np.log(p); q=p
        while q<=N: lam[q]=lp; q*=p
    return lam
eps=0.25; alpha=0.5+eps
th=np.linspace(10,60,501); s=alpha+1j*th
exact=np.array([float((b/a).real) for a,b in (xi_pair(S,acb(eps,t)) for t in th)])
arch=np.real(1/s+1/(s-1)-0.5*np.log(np.pi))+np.array([float(mp.re(mp.digamma(z/2)))/2 for z in s])
out={}
for x in [30,300,3000]:
    N=int(x*x); lam=vm(N); n=np.nonzero(lam)[0]; L=lam[n]
    w=np.where(n<=x,L,L*np.log(x*x/n)/np.log(x))
    lx=np.log(x)
    zp=np.array([-np.sum(w*np.exp(-z*np.log(n))) for z in s])
    zp+= (x**(2*(1-s))-x**(1-s))/((1-s)**2*lx)
    zp+= sum((x**(-2*q-s)-x**(-2*(2*q+s)))/((2*q+s)**2) for q in range(1,60))/lx
    rho=0.5+1j*g
    zs=np.array([np.sum((x**(rho-z)-x**(2*(rho-z)))/(z-rho)**2 + (x**(np.conj(rho)-z)-x**(2*(np.conj(rho)-z)))/(z-np.conj(rho))**2)/lx for z in s])
    approx=arch+np.real(zp); full=approx+np.real(zs)
    e1=np.abs(approx-exact); e2=np.abs(full-exact)
    out['curve_%d'%x]=list(e1); out['curvez_%d'%x]=list(e2)
    out[x]=dict(sup_primes=float(e1.max()),rel=float((e1/exact).max()),sup_with_zeros=float(e2.max()),xeps=x**(-eps))
    print(x, 'primes only sup err',e1.max(),'max rel',(e1/exact).max(),' with zero sum',e2.max(),' x^-eps',x**(-eps), flush=True)
out['theta']=list(th); out['exact']=list(exact)
json.dump(out,open('selberg.json','w'))
