import json, numpy as np, mpmath as mp, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':0.25,'axes.spines.top':False,'axes.spines.right':False})
mp.mp.dps=20
S=json.load(open('salem_riesz.json')); E=json.load(open('riesz_explicit.json')); Q=json.load(open('riesz_thorin.json'))
g=np.load('g_1_400.npy')[:60]
coef=[complex(mp.gamma(1-mp.mpc(0.5,t)/2)/(2*mp.zeta(mp.mpc(0.5,t),derivative=1))) for t in g]
triv=E['triv']
pred=lambda x: sum(2*(c*x**(0.5j*t)).real for c,t in zip(coef,g))+sum(r*x**(-k) for k,r in triv)/x**0.25
xs=np.array(S['riesz']['x']); R=np.array(S['riesz']['R']); m=(xs>=1e3)&(xs<=1e10)
xf=np.geomspace(1e3,1e10,3000)
fig,ax=plt.subplots(1,3,figsize=(12,3.5))
ax[0].semilogx(xf,[pred(x) for x in xf],color='#c0392b',lw=0.8,label='explicit formula, 60 zeros')
ax[0].semilogx(xs[m],R[m]/xs[m]**0.25,'o',ms=3,color='#1b4f72',label=r'M\"obius series')
ax[0].set_xlabel('x'); ax[0].set_title(r'$R(x)/x^{1/4}$',fontsize=10); ax[0].legend(fontsize=7)
C=Q['coef']; k=np.array([c['k'] for c in C])
ax[1].semilogy(k,[1-c['exact'] for c in C],'o-',color='#1b4f72',label=r'$1-1/\zeta(2k)$')
ax[1].semilogy(k,4.0**(-k),':',color='gray',label=r'$4^{-k}$')
ax[1].errorbar(k,[max(1-c['mc'],1e-6) for c in C],yerr=[c['se'] for c in C],fmt='s',ms=4,color='#c0392b',label='via the clock (Monte Carlo)')
ax[1].set_xlabel('k'); ax[1].set_title(r'Riesz coefficients from the clock',fontsize=10); ax[1].legend(fontsize=7)
X=np.array([r['x'] for r in Q['canc']])
ax[2].loglog(X,[r['raw'] for r in Q['canc']],color='k',lw=1.4,label='power series')
for Tv,c in [('0.01','#c0392b'),('0.023','#1b4f72'),('0.05','#27ae60')]:
    ax[2].loglog(X,[r['T'+Tv] for r in Q['canc']],color=c,lw=1.2,label=r'clock series, $T=%s$'%Tv)
ax[2].set_xlabel('x'); ax[2].set_ylabel(r'$\log_{10}$ of the largest term'); ax[2].legend(fontsize=7); ax[2].set_title('Cancellation needed',fontsize=10)
fig.tight_layout(); fig.savefig('fig_riesz.pdf'); print('ok')
