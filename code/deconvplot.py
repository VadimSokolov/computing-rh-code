import numpy as np, json, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.optimize import nnls
plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':0.25,'axes.spines.top':False,'axes.spines.right':False})
g=np.concatenate([np.load(f) for f in ['g_1_400.npy','g_401_1000.npy','g_1001_1700.npy']])
d=np.load('eps_0.1.npz'); th=d['theta']; F=d['F']/np.pi; eps=0.1
P=lambda a,x: a/(np.pi*(a*a+x*x))
m=(th>=15)&(th<=55); th=th[m][::2]; y=F[m][::2]
far=g[(g<5)|(g>65)]; bg=np.array([np.sum(P(eps,t-far)+P(eps,t+far)) for t in th])
xs=np.arange(5,65,0.02); A=P(eps,th[:,None]-xs[None,:])+P(eps,th[:,None]+xs[None,:])
w,res=nnls(A,y-bg,maxiter=20000)
D=json.load(open('deconv.json'))
fig,ax=plt.subplots(1,2,figsize=(9.5,3.5))
ml,sl,bl=ax[0].stem(xs[w>1e-4],w[w>1e-4]); plt.setp(sl,color='#1b4f72',lw=1); plt.setp(ml,color='#1b4f72',ms=3); plt.setp(bl,visible=False)
for z in g[(g>5)&(g<65)]: ax[0].axvline(z,color='#c0392b',ls='--',lw=0.7)
ax[0].axvspan(15,55,color='gray',alpha=0.08); ax[0].set_xlabel(r'candidate atom location'); ax[0].set_title('Nonnegative deconvolution of $\\rho_{0.1}$ (dashed: zeros)')
dl=[0.005,0.01,0.02,0.05]; r0=D['0.0']['res']; ex=[np.sqrt(max(D[str(x)]['res']**2-r0**2,1e-12)) for x in dl]
ax[1].loglog(dl,ex,'o-',lw=1.8,label='excess residual'); ax[1].loglog(dl,1.385*(np.array(dl)/0.05)**2,'--',color='gray',label=r'$\propto\delta^2$')
ax[1].axhline(r0,color='#c0392b',lw=0.8,ls=':',label='off grid baseline')
ax[1].set_xlabel(r'displacement $\delta$ of the fifth zero'); ax[1].legend(fontsize=8); ax[1].set_title('Posterior predictive check')
fig.tight_layout(); fig.savefig('book/fig/deconv.pdf'); print('ok')
