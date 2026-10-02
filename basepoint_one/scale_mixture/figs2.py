import numpy as np, json, mpmath as mp, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':0.25,'axes.spines.top':False,'axes.spines.right':False})
S=json.load(open('strip.json')); C=json.load(open('comp.json'))
mp.mp.dps=15
xi=lambda s: s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
dl=lambda s: mp.diff(lambda z: mp.log(xi(z)), s)
cols={'0.5':'#27ae60','0.6':'#c0392b','0.75':'#d68910','1.0':'#1b4f72'}
fig,ax=plt.subplots(1,3,figsize=(11.5,3.5))
th=np.linspace(0,40,801)
for al in ['0.6','0.75','1.0']:
    r=[float(mp.re(dl(mp.mpc(float(al),t)))/mp.pi) if t>0 else float(mp.re(dl(mp.mpf(al)))/mp.pi) for t in th]
    ax[0].semilogy(th,r,color=cols[al],lw=1.0,label=r'$\alpha=%s$'%al)
    ax[0].axhline(S['dens'][al]['b']/np.pi,color=cols[al],lw=0.7,ls=':')
ax[0].set_xlabel(r'$\theta$'); ax[0].legend(fontsize=7); ax[0].set_title(r'$\rho_\alpha$ in the strip, Cauchy floors $b_\alpha/\pi$',fontsize=9)
for al in ['0.6','0.75','1.0']:
    d=S['dens'][al]; y=np.array(d['y']); f=np.array(d['f']); m=y>=1
    if al=='0.5': ax[1].semilogy(y[y<=6],np.maximum(f[y<=6],1e-12),color=cols[al],lw=1.2,label=r'centre, $f$ itself')
    else: ax[1].semilogx(y[m],y[m]**2*f[m]/(d['b']/np.pi),color=cols[al],lw=1.2,label=r'$\alpha=%s$'%al)
ax[1].axhline(1,color='gray',lw=0.7,ls=':'); ax[1].set_xlabel('y'); ax[1].set_title(r'$\pi y^2 f_{Y_\alpha}(y)/b_\alpha\to1$',fontsize=9); ax[1].legend(fontsize=7)
ax[1].set_ylim(0.95,1.2); ax[1].set_xlim(1,50)
K=C['K']; y=np.array([q['y'] for q in K])
ax[2].loglog(y,[q['K1'] for q in K],'o-',color='#1b4f72',label=r'$K_1(y)$')
ax[2].loglog(y,[q['K1_asym'] for q in K],':',color='#c0392b',label=r'$b_1/(\pi y)$')
kc=np.array([q['Kcentre'] for q in K]); mm=kc>1e-12
ax[2].loglog(y[mm],kc[mm],'s--',color='#27ae60',label=r'$K_{1/2}(y)=\sum e^{-\gamma y}$')
ax[2].set_xlabel('y'); ax[2].legend(fontsize=7); ax[2].set_title(r'Poisson traces',fontsize=9)
fig.tight_layout(); fig.savefig('fig_strip.pdf'); print('ok')
