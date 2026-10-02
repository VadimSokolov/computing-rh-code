import sys, numpy as np, mpmath as mp
# python3 weil.py     weights log p in double precision, 20 digits (Table sb:tab:bdlp)
# python3 weil.py mp  weights log p in multiple precision, 30 digits
MP=len(sys.argv)>1 and sys.argv[1]=='mp'
mp.mp.dps=30 if MP else 20
g=np.concatenate([np.load(f) for f in ['g_1_400.npy','g_401_1000.npy','g_1001_1700.npy']])
import sympy
P=list(sympy.primerange(2,20000))
lam={}
for p in P:
    pk=p
    while pk<20000: lam[pk]=mp.log(p) if MP else np.log(p); pk*=p
def weil(t,deriv=False):
    t=mp.mpf(t)
    if not deriv:
        h=lambda r: mp.e**(-t*r*r); hpole=2*mp.e**(t/4)
        g0=1/mp.sqrt(4*mp.pi*t); gf=lambda x: mp.e**(-x*x/(4*t))/mp.sqrt(4*mp.pi*t)
    else: # -d/dt
        h=lambda r: r*r*mp.e**(-t*r*r); hpole=2*(-mp.mpf(1)/4)*mp.e**(t/4)
        gf=lambda x: mp.diff(lambda tt: -mp.e**(-x*x/(4*tt))/mp.sqrt(4*mp.pi*tt), t)
        g0=gf(0)
    arch=mp.quad(lambda r: h(r)*mp.re(mp.digamma(mp.mpf(1)/4+1j*r/2)),[-mp.inf,-20,0,20,mp.inf])/(2*mp.pi)
    pr=2*sum(l/mp.sqrt(n)*gf(mp.log(n)) for n,l in lam.items())
    return (hpole-g0*mp.log(mp.pi)+arch-pr)/2   # sum over all rho counts +-gamma; divide by 2
for t in [0.005,0.01,0.02,0.05,0.1]:
    Wz=np.sum(np.exp(-g*g*t)); dWz=np.sum(g*g*np.exp(-g*g*t))
    print(t,'W zeros %.12e weil %s | -W\' zeros %.12e weil %s'%(Wz,mp.nstr(weil(t),13),dWz,mp.nstr(weil(t,True),13)))
