import numpy as np, json, mpmath as mp, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':0.25,'axes.spines.top':False,'axes.spines.right':False})
A=json.load(open('alpha1.json')); b1=json.load(open('levy.json'))['b1']
fig,ax=plt.subplots(1,2,figsize=(9.8,3.7))
x=np.linspace(0.02,4.2,600); ac=np.exp(-3*x)/(1-np.exp(-2*x))-1
ax[0].plot(x,ac,color='#1b4f72',lw=1.5,label=r'$\frac{e^{-3x}}{1-e^{-2x}}-1$ (gamma factor minus pole)')
lam={}
for p in [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59]:
    pk=p
    while pk<=60: lam[pk]=np.log(p); pk*=p
for n_,l in lam.items(): ax[0].vlines(np.log(n_),0,l/n_,color='#c0392b',lw=1.5)
ax[0].plot([],[],color='#c0392b',label=r'prime atoms $\Lambda(n)/n$ at $\log n$')
ax[0].axhline(0,color='k',lw=0.6); ax[0].set_ylim(-1.05,1.2); ax[0].set_xlabel('x'); ax[0].legend(fontsize=7,loc='upper right')
ax[0].set_title(r'Signed canonical measure $x\,\Pi_1(dx)$ of $\log\xi(1+s)$',fontsize=10)
th=np.array(A['grid']['th']); r=np.array(A['grid']['rho'])
Ad=np.array([float((mp.re(mp.digamma(mp.mpc(1.5,t/2)))-mp.digamma(mp.mpf(1.5)))/2) for t in th])
P=np.pi*r-b1-Ad; m=th<=60
ax[1].plot(th[m],np.pi*r[m],color='k',lw=1.4,label=r'$\pi\rho_1(\theta)$')
ax[1].plot(th[m],Ad[m],color='#1b4f72',lw=1.2,label=r'gamma factor $A(\theta)$')
ax[1].plot(th[m],P[m],color='#c0392b',lw=1.0,label=r'primes minus pole $P(\theta)$')
ax[1].axhline(b1,color='#27ae60',lw=1.2,label=r'drift $b_1=\sum_\rho\rho^{-1}$ (stable part)')
ax[1].set_xlabel(r'$\theta$'); ax[1].legend(fontsize=7); ax[1].set_title('Thorin density from the sine squared lemma',fontsize=10)
fig.tight_layout(); fig.savefig('fig_levy.pdf'); print('ok', P.min(), P.max())
