# Proposed code/perturb.py (written in round 2 of the editing, misc/repro/reconstruct/p13_perturb_fig/code_proposal/;
# it replaces section 8 of code/plots3.py, which plots an empty panel). Figure fig:ch13:perturb, fig/perturb.pdf: the
# relative change |W_F-W|/W of the heat trace when the closest pair of zeros in (1000, 2000) is moved to distance 1/4
# from the critical line (Theorem thm:ch13:perturb). W_F-W is evaluated in mpmath at 40 digits, because in double
# precision its terms underflow for t > 2e-4; W is the sum over the first 1700 zeros, as in the archived figure, which
# this script redraws exactly (with Matplotlib 3.10.8 the page content of the PDF is byte for byte the archived one).
# perturb.json records the pair, the curve, the maximum, the zero of W_F-W, where the change falls below 1e-1000, and,
# at a few small t, W from the explicit formula eq:ch2:Warith, which the sum over 1700 zeros falls short of below
# t = 1e-6. Writes perturb.json and fig/perturb.pdf; under three minutes on one core, most of it for the explicit formula.
import json, os, numpy as np, mpmath as mp, matplotlib
from sympy import primerange
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':0.25,'axes.spines.top':False,'axes.spines.right':False})
mp.mp.dps=40
g=np.concatenate([np.load(f) for f in ['g_1_400.npy','g_401_1000.npy','g_1001_1700.npy']])
gp=g[(g>1000)&(g<2000)]; k=int(np.argmin(np.diff(gp))); A_,B_=gp[k],gp[k+1]; gg=(A_+B_)/2
A,B,G=mp.mpf(float(A_)),mp.mpf(float(B_)),mp.mpf(float(gg)); a=G**2-mp.mpf(1)/16-1j*G/2   # a=-w^2, w=1/4+ig
G2=[mp.mpf(float(x))**2 for x in g]
dW=lambda t: 2*mp.re(mp.exp(-a*t))-mp.exp(-A*A*t)-mp.exp(-B*B*t)
W=lambda t: mp.fsum(mp.exp(-x*t) for x in G2)
rel=lambda t: mp.log10(abs(dW(t))/W(t))
PP=[(mp.log(p),p**j) for p in primerange(2,1000) for j in range(1,10) if p**j<1000]
def Wexp(t):   # e^{t/4} + archimedean integral - prime sum; prime powers below 1000 suffice for t <= 0.01
    m=lambda r: mp.re(mp.digamma(mp.mpf(1)/4+0.5j*r))-mp.log(mp.pi)
    R=14/mp.sqrt(t); pts=[mp.mpf(0)]+[mp.mpf(2)**j for j in range(-1,80) if mp.mpf(2)**j<R]+[R]
    ar=mp.quad(lambda r: mp.exp(-r*r*t)*m(r),pts)/(2*mp.pi)
    pr=mp.fsum(L/mp.sqrt(n)*mp.exp(-mp.log(n)**2/(4*t)) for L,n in PP)/(2*mp.sqrt(mp.pi*t))
    return mp.exp(t/4)+ar-pr
tt=np.logspace(-8,-2,400); y=[float(rel(mp.mpf(t))) for t in tt]
def argmax(f,lo=mp.mpf(-6.6),hi=mp.mpf(-5.6)):   # maximum of f(log10 t), by golden section
    q=(mp.sqrt(5)-1)/2; u1,u2=hi-q*(hi-lo),lo+q*(hi-lo); f1,f2=f(u1),f(u2)
    for _ in range(60):
        if f1>f2: hi,u2,f2=u2,u1,f1; u1=hi-q*(hi-lo); f1=f(u1)
        else: lo,u1,f1=u1,u2,f2; u2=lo+q*(hi-lo); f2=f(u2)
    return 10**((lo+hi)/2)
tmax=argmax(lambda u: rel(10**u)); tm2=argmax(lambda u: mp.log10(abs(dW(10**u))/Wexp(10**u)))
t0=mp.findroot(dW,(mp.mpf('5e-8'),mp.mpf('1.5e-7')),solver='anderson')   # W_F=W here, near 1/(2g^2)
t1000=10**mp.findroot(lambda u: rel(10**u)+1000,-3.2)
chk={'%g'%t: dict(W1700=float(W(mp.mpf(t))),Wexp=float(Wexp(mp.mpf(t))),log10_rel_1700=float(rel(mp.mpf(t))),
     log10_rel_exp=float(mp.log10(abs(dW(mp.mpf(t)))/Wexp(mp.mpf(t))))) for t in [1e-8,1e-7,7e-7,1e-6,1e-5]}
out=dict(pair=[float(A_),float(B_)],index=[int(np.nonzero(g==A_)[0][0])+1,int(np.nonzero(g==B_)[0][0])+1],g=float(gg),
         t=tt.tolist(),log10_rel=y,max_rel=float(10**rel(tmax)),t_max=float(tmax),
         max_rel_Wexp=float(abs(dW(tm2))/Wexp(tm2)),t_max_Wexp=float(tm2),
         t_zero=float(t0),t_below_1e_minus_1000=float(t1000),W_check=chk)
json.dump(out,open('perturb.json','w'),indent=1)
fig,ax=plt.subplots(figsize=(6.5,3.5))
ax.semilogx(tt,y,lw=1.6)
ax.set_ylim(-3000,5); ax.set_xlabel('t'); ax.set_ylabel(r'$\log_{10}|W_F-W|/W$')
ax.set_title('Moving the pair near height %.0f off the line:\nrelative change in the heat trace'%gg,fontsize=9)
os.makedirs('fig',exist_ok=True); fig.tight_layout(); fig.savefig('fig/perturb.pdf'); plt.close()
print('pair',A_,B_,'max %.4g at t=%.4g (W from the explicit formula: %.4g at t=%.4g), W_F=W at t=%.4g, below 1e-1000 from t=%.4g'
      %(out['max_rel'],out['t_max'],out['max_rel_Wexp'],out['t_max_Wexp'],out['t_zero'],out['t_below_1e_minus_1000']))
