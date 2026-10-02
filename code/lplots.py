exec(open('lfun.py').read().split('res={}')[0])
mp.mp.dps=40
import json, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':0.25,'axes.spines.top':False,'axes.spines.right':False})
S=json.load(open('lsin2.json')); Lz=json.load(open('lfun.json'))
setups={'zeta':(Phi,0.5,3.0),'chi12':(K12,2,5.0),'Delta':(KD,1,4.0)}
lab={'zeta':r'$\zeta$','chi12':r'$L(s,\chi_{12})$','Delta':r'$L(s,\Delta)$'}
col={'zeta':'#1b4f72','chi12':'#c0392b','Delta':'#27ae60'}
fig,ax=plt.subplots(1,2,figsize=(9.8,3.7))
for n in ['zeta','chi12','Delta']:
    R=np.array(S[n]['rows']); ax[0].plot(R[:,0],R[:,5],color=col[n],lw=1.3,label='prime piece, '+lab[n])
R=np.array(S['zeta']['rows']); ax[0].plot(R[:,0],R[:,4],'--',color=col['zeta'],lw=1.2,label=r'pole piece, $\zeta$')
ax[0].axhline(0,color='k',lw=0.6); ax[0].set_xlabel(r'$\theta$'); ax[0].set_title(r'Sine squared pieces at $\sigma=2$',fontsize=10); ax[0].legend(fontsize=7)
eps=mp.mpf('0.1'); th=np.arange(0.5,40,0.05)
out={}
for n in ['chi12','Delta']:
    K,scale,U=setups[n]; h=mp.mpf(1)/80; us,Ks=nodes(K,h,U); uKs=[u/scale*k for u,k in zip(us,Ks)]
    r=[float(mp.re(Lam(us,uKs,h,eps+1j*mp.mpf(t),scale)/Lam(us,Ks,h,eps+1j*mp.mpf(t),scale))/mp.pi) for t in th]
    out[n]=r; ax[1].semilogy(th,r,color=col[n],lw=1.1,label=r'$\rho_{0.1}$, '+lab[n])
    for z in Lz[n]['zeros']: ax[1].axvline(z,color=col[n],lw=0.3,alpha=0.5)
ax[1].set_xlim(0,40); ax[1].set_xlabel(r'$\theta$'); ax[1].set_title(r'Thorin densities at $\varepsilon=0.1$',fontsize=10); ax[1].legend(fontsize=7)
fig.tight_layout(); fig.savefig('book/fig/lfunc.pdf')
print('min rho chi12 %.4f Delta %.4f'%(min(out['chi12']),min(out['Delta'])))
json.dump(out,open('lrho.json','w'))
