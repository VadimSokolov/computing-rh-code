import json, numpy as np, mpmath as mp, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from flint import acb, arb, ctx
import tilted
plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':0.25,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':200,'axes.prop_cycle':plt.cycler(color=['#1b4f72','#c0392b','#27ae60','#d68910','#7d3c98','#17a589','#566573'])})
F='book/fig/'
g=np.concatenate([np.load(f) for f in ['g_1_400.npy','g_401_1000.npy','g_1001_1700.npy']])
P=lambda a,x: a/(np.pi*(a*a+x*x))
# 1 Polya kernel on lines in the strip
ctx.prec=100
fig,ax=plt.subplots(1,2,figsize=(9.5,3.5))
u=np.linspace(-1.2,1.2,601)
Phi=[float(tilted.Phi_c(acb(float(x),0),120).real) for x in u]
ax[0].plot(u,Phi,lw=2); ax[0].fill_between(u,0,Phi,alpha=0.15)
ax[0].set_title('Pólya kernel $\Phi(u)$'); ax[0].set_xlabel('u')
xx=np.linspace(-3,3,601)
for i,d in enumerate([0.6854,0.3927,0.1,0.02]):
    eta=arb.pi()/4-arb(d)
    v=[abs(complex(tilted.Phi_c(acb(float(x),eta),120))) for x in xx]
    ax[1].semilogy(xx,np.maximum(v,1e-30),lw=1.4,label=r'$\eta=%.4f$'%(np.pi/4-d))
ax[1].set_ylim(1e-12,1e5); ax[1].set_title(r'$|\Phi(x+i\eta)|$ inside the strip'); ax[1].set_xlabel('x'); ax[1].legend(fontsize=8)
fig.tight_layout(); fig.savefig(F+'kernel.pdf'); plt.close()
# 2 heat trace and arithmetic components
A=json.load(open('arith.json')); t=np.array([r[0] for r in A])
Pi=np.array([float(r[1]) for r in A]); Ar=np.array([float(r[2]) for r in A]); Pr=np.array([float(r[3]) for r in A]); Wz=np.array([float(r[5]) for r in A]); dig=np.array([r[6] for r in A])
tt=np.linspace(0.002,0.3,300); Wt=[np.sum(np.exp(-g**2*x)) for x in tt]
fig,ax=plt.subplots(1,2,figsize=(9.5,3.5))
ax[0].semilogy(tt,Wt,lw=2,label='W(t) from 1700 zeros')
ax[0].semilogy(tt,np.exp(-g[0]**2*tt),'--',lw=1,label=r'$e^{-\gamma_1^2 t}$')
ax[0].semilogy(t,Pi,'o',ms=4,label='polar'); ax[0].semilogy(t,np.abs(Ar),'s',ms=4,label='|archimedean|'); ax[0].semilogy(t[t>0.004],Pr[t>0.004],'^',ms=4,label='prime')
ax[0].set_xlabel('t'); ax[0].legend(fontsize=7); ax[0].set_title('Heat trace and its arithmetic pieces')
ax[1].plot(t,dig,'o-',lw=1.8); ax[1].plot(tt,tt*g[0]**2/np.log(10),'--',color='gray',label=r'$\gamma_1^2t/\log 10$')
ax[1].set_xlabel('t'); ax[1].set_ylabel('decimal digits lost'); ax[1].legend(fontsize=8); ax[1].set_title('Cancellation in the explicit formula')
fig.tight_layout(); fig.savefig(F+'heattrace.pdf'); plt.close()
# 3 precision budget
fig,ax=plt.subplots(1,2,figsize=(9.5,3.5))
th=np.linspace(0,8000,200)
ax[0].plot(th,0.341*th,lw=2,label='untilted, $0.341\\theta$')
for d in [0.1,0.02,0.01]: ax[0].plot(th,d*th/np.log(10),lw=1.6,label=r'tilted, $\delta=%s$'%d)
ax[0].scatter([1000,5000,7005],[0.1*1000/np.log(10),0.02*5000/np.log(10),0.01*7005/np.log(10)],color='k',zorder=5,s=18,label='runs in the book')
ax[0].set_yscale('log'); ax[0].set_ylim(0.5,5000); ax[0].set_xlabel(r'height $\theta$'); ax[0].set_ylabel('digits lost'); ax[0].legend(fontsize=7); ax[0].set_title('Precision cost of height')
dd=np.array([0.6854,0.3927,0.2,0.1,0.05,0.02,0.01]); M=np.array([0.554,2.52,19.0,126,774,7836,43503])
ax[1].loglog(dd,M,'o-',lw=1.8,label=r'$M_\eta$'); ax[1].loglog(dd,126*(dd/0.1)**-2.5,'--',color='gray',label=r'$\propto\delta^{-5/2}$')
ax[1].set_xlabel(r'$\delta=\pi/4-\eta$'); ax[1].legend(fontsize=8); ax[1].set_title('Size of the kernel near the edge of the strip')
fig.tight_layout(); fig.savefig(F+'budget.pdf'); plt.close()
# 4 decomposition near gamma1
eps=0.1; x=np.linspace(10,24,1400)
tot=np.sum(P(eps,x[:,None]-g[None,:])+P(eps,x[:,None]+g[None,:]),axis=1)
own=P(eps,x-g[0]); nb=P(eps,x-g[1]); rest=tot-own-nb
fig,ax=plt.subplots(figsize=(7.5,3.5))
ax.semilogy(x,tot,lw=2,label=r'$\rho_{0.1}$ (all zeros)'); ax.semilogy(x,own,'--',lw=1.2,label=r'kernel of $\gamma_1$'); ax.semilogy(x,nb,':',lw=1.4,label=r'kernel of $\gamma_2$'); ax.semilogy(x,rest,lw=1.2,label='all other zeros')
ax.set_xlabel(r'$\theta$'); ax.legend(fontsize=8); ax.set_title('Poisson structure of the Thorin density near the first two zeros')
fig.tight_layout(); fig.savefig(F+'decomp.pdf'); plt.close()
# 5 sine squared decomposition
R=json.load(open('results_prime.json')); R2=[r for r in R if abs(r['b']-2)<1e-9]
al=[r['alpha'] for r in R2]; w=0.2; xi=np.arange(len(al))
fig,ax=plt.subplots(1,2,figsize=(9.5,3.5))
ax[0].bar(xi-1.5*w,[r['c'] for r in R2],w,label='drift'); ax[0].bar(xi-0.5*w,[r['gamma'] for r in R2],w,label='archimedean')
ax[0].bar(xi+0.5*w,[r['pole'] for r in R2],w,label='pole'); ax[0].bar(xi+1.5*w,[r['primes'] for r in R2],w,label='primes')
ax[0].plot(xi,[r['exact'] for r in R2],'kD',ms=6,label=r'$v_\alpha(2)$')
ax[0].axhline(0,color='k',lw=0.6); ax[0].set_xticks(xi); ax[0].set_xticklabels([r'$\alpha=%g$'%a for a in al]); ax[0].legend(fontsize=7); ax[0].set_title(r'Sine squared pieces at $\theta=2$')
aa=np.linspace(1.05,3,100)
pole=aa**0*4/((aa-1)*((aa-1)**2+4))
ax[1].semilogy(aa,pole,lw=2,label='|pole piece| at $\\theta=2$')
kap=[(abs(r['pole'])+abs(r['primes']))/r['exact'] for r in R2]
ax[1].semilogy(al,kap,'o',ms=7,color='#c0392b',label=r'condition number $\kappa_\alpha(2)$')
ax[1].set_xlabel(r'$\alpha$'); ax[1].legend(fontsize=8); ax[1].set_title('The cancellation grows like $1/(\\alpha-1)$')
fig.tight_layout(); fig.savefig(F+'sin2pieces.pdf'); plt.close()
# 6 Turing and Thorin: M_eps - N
C=np.load('plotcurves.npz'); Z=C['zeros']
fig,ax=plt.subplots(1,2,figsize=(9.5,3.5))
for i,e in enumerate([0.5,0.1,0.02]):
    th_=C['th_%s'%e]; Mv=C['M_%s'%e]; N=np.searchsorted(Z,th_)
    m=(th_>20)&(th_<60); ax[0].plot(th_[m],Mv[m]-N[m],lw=1.1,label=r'$\varepsilon=%s$'%e)
ax[0].axhline(0,color='k',lw=0.6); ax[0].set_xlabel(r'$\theta$'); ax[0].set_title(r'$M_\varepsilon(\theta)-N(\theta)$: Poisson averaged counting error'); ax[0].legend(fontsize=8)
TI=json.load(open('turing_ig.json'))['zd']; ee=[r['eps'] for r in TI]
ax[1].semilogx(ee,[r['M'] for r in TI],'o-',label=r'$M_\varepsilon(200)$'); ax[1].semilogx(ee,[79+r['top'] for r in TI],'s-',label=r'$79+\tau_\varepsilon(200)$'); ax[1].semilogx(ee,[r['sum'] for r in TI],'k^-',label=r'$M_\varepsilon+\tau_\varepsilon$')
ax[1].set_xlabel(r'$\varepsilon$'); ax[1].legend(fontsize=8); ax[1].set_title('The zero density identity returns the integer 79')
fig.tight_layout(); fig.savefig(F+'turingthorin.pdf'); plt.close()
# 7 Hankel minors and derivative signs
H=json.load(open('hankel.json'))
fig,ax=plt.subplots(1,2,figsize=(9.5,3.5))
for key,lab in [('0.01','t = 0.01'),('0.05','t = 0.05')]:
    ax[0].semilogy([2,3,4,5],[float(v) for v in H[key]['even']],'o-',label='even, '+lab); ax[0].semilogy([2,3,4,5],[float(v) for v in H[key]['odd']],'s--',label='odd, '+lab)
ax[0].set_xlabel('block size N'); ax[0].set_title('Normalised Hankel minors of the heat trace'); ax[0].legend(fontsize=7)
tt=np.linspace(0.005,0.1,200)
for k in range(0,6):
    ck=[np.sum(g**(2*k)*np.exp(-g**2*x)) for x in tt]
    ax[1].semilogy(tt,np.array(ck)/np.array(ck)[0],lw=1.3,label=r'$c_%d(t)/c_%d(0.005)$'%(k,k))
ax[1].set_xlabel('t'); ax[1].legend(fontsize=7); ax[1].set_title(r'$(-1)^kW^{(k)}$ are positive and log convex')
fig.tight_layout(); fig.savefig(F+'hankel.pdf'); plt.close()
# 9 cell leaks
S=json.load(open('spacing.json')); lk=np.array(S['leaks']); gw=g[g<200]
fig,ax=plt.subplots(figsize=(7.5,3.3))
ax.bar(np.arange(1,len(lk)+1),lk,color=np.where(lk>0,'#1b4f72','#c0392b'))
ax.set_xlabel('zero index k'); ax.set_ylabel(r'$(\mu_\varepsilon(C_k)-1)/\varepsilon$'); ax.set_title('First order cell errors below height 200: the imbalance of the spacings')
fig.tight_layout(); fig.savefig(F+'cells.pdf'); plt.close()
# 11 Li test functions
th=np.linspace(0.5,200,2000)
fig,ax=plt.subplots(figsize=(7.5,3.4))
for n in [1,10,50,100,200]: ax.semilogy(th,4*np.sin(n*np.arctan(1/(2*th)))**2+1e-12,lw=1.2,label='n = %d'%n)
for z in g[g<200]: ax.axvline(z,color='gray',lw=0.3,alpha=0.5)
ax.set_ylim(1e-6,5); ax.set_xlabel(r'$\theta=\sqrt{z}$'); ax.legend(fontsize=7,ncol=5); ax.set_title(r'Li test functions $4\sin^2(n\arctan(1/2\theta))$ against the zeros')
fig.tight_layout(); fig.savefig(F+'litest.pdf'); plt.close()
# 12 benchmark vs zeta cell TV
sb=json.load(open('sinh_bench.json')); an=json.load(open('analysis.json'))['rows']
fig,ax=plt.subplots(figsize=(6.5,3.5))
ax.loglog([r['eps'] for r in sb],[r['celltv']/50 for r in sb],'o-',label='sinh clock, 50 equally spaced atoms')
ax.loglog([r['eps'] for r in an],[r['tv_cells']/79 for r in an],'s-',label='zeta zeros, 79 atoms')
ax.set_xlabel(r'$\varepsilon$'); ax.set_ylabel('cell distance per atom'); ax.legend(fontsize=8); ax.set_title('Irregular spacing costs a factor of three')
fig.tight_layout(); fig.savefig(F+'bench.pdf'); plt.close()
print('done')
