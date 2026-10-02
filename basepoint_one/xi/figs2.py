import numpy as np, json, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':0.25,'axes.spines.top':False,'axes.spines.right':False})
I=json.load(open('improve.json')); g=np.concatenate([np.load(f) for f in ['g_1_400.npy','g_401_1000.npy','g_1001_1700.npy']])
x=np.array(I['talbot']['x']); f=np.array(I['talbot']['f']); sd=np.array(I['subdens'])
fig,ax=plt.subplots(1,2,figsize=(9.8,3.6))
m=f>1e-12
ax[0].loglog(x[m],f[m],color='#1b4f72',lw=1.6,label='Talbot inversion of $\\varphi_1$')
ms=sd>1e-12; ax[0].loglog(x[ms][::2],sd[ms][::2],'o',ms=4,color='#c0392b',label='subordination formula')
xx=np.logspace(0,3,20); ax[0].loglog(xx,0.006515*xx**-1.5,'--',color='gray',label=r'$0.00652\,x^{-3/2}$')
ax[0].set_xlabel('x'); ax[0].set_title(r'Density of $X_1$',fontsize=10); ax[0].legend(fontsize=7)
th=np.array(I['alphas']['th']); cols={'0.6':'#c0392b','0.75':'#d68910','1.0':'#1b4f72'}
for a,r in I['alphas']['rows'].items(): ax[1].semilogy(th,r['rho'],color=cols[a],lw=1.1,label=r'$\alpha=%s$'%a)
for z in g[g<60]: ax[1].axvline(z,color='gray',lw=0.4,alpha=0.6)
ax[1].set_xlabel(r'$\theta$'); ax[1].set_title(r'Thorin densities $\rho_\alpha$ of $\xi(\alpha)/\xi(\alpha+\sqrt{s})$',fontsize=10); ax[1].legend(fontsize=7)
fig.tight_layout(); fig.savefig('fig_density2.pdf'); print('ok')
