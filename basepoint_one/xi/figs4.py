import numpy as np, json, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':0.25,'axes.spines.top':False,'axes.spines.right':False})
B=json.load(open('below.json'))
fig,ax=plt.subplots(1,2,figsize=(9.8,3.7))
cols={1.0:'#1b4f72',0.9:'#27ae60',0.75:'#d68910',0.6:'#c0392b'}
for al,c in cols.items():
    r=[x for x in B if x['alpha']==al and abs(x['theta']-14.134725)<1e-6]
    X=np.array([x['X'] for x in r]); e=np.abs([x['err'] for x in r])
    ax[0].loglog(X,e,'o-',color=c,label=r'$\alpha=%s$'%al)
    ax[0].loglog(X,e[-1]*(X/X[-1])**(0.5-al),':',color=c,lw=0.9)
ax[0].set_xlabel('prime cutoff X'); ax[0].set_ylabel(r'error at $\theta=\gamma_1$'); ax[0].legend(fontsize=7)
ax[0].set_title(r'Rate $X^{1/2-\alpha}$ (dotted): convergence slows as $\alpha\downarrow\frac{1}{2}$',fontsize=10)
r=[x for x in B if x['alpha']==0.5 and abs(x['theta']-14.134725)<1e-6]
X=np.array([x['X'] for x in r]); t=np.array([x['total'] for x in r])
ax[1].semilogx(X,t,'o-',color='#7d3c98',label=r'compensated sum at $\alpha=\frac{1}{2}$, $\theta=\gamma_1$')
ax[1].semilogx(X,np.log(X)-np.log(X[0])+t[0],'--',color='gray',label=r'slope $\log X$')
ax[1].set_xlabel('prime cutoff X'); ax[1].legend(fontsize=7); ax[1].set_title('At the centre the sum diverges at a zero: the atom forming',fontsize=10)
fig.tight_layout(); fig.savefig('fig_below.pdf'); print('ok')
