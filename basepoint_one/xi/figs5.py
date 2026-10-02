import numpy as np, json, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':0.25,'axes.spines.top':False,'axes.spines.right':False})
R=json.load(open('rs1.json')); D=json.load(open('rs1d.json'))
t=np.array(R['t']); b=np.array(R['b']); bt=np.array(R['bt'])
fig,ax=plt.subplots(1,2,figsize=(9.8,3.7))
for (a,c_,col,l) in [(0.006,0.011,'#d5f5e3','small time'),(0.011,0.05,'#fdebd0','bulk'),(0.05,100,'#fadbd8','right tail')]:
    ax[0].axvspan(a,c_,color=col,alpha=0.7,lw=0)
ax[0].semilogx(t,b,color='#1b4f72',lw=2,label='exact boundary (Volterra)')
ax[0].semilogx(t,bt,'--',color='#c0392b',lw=1.3,label='tangent approximation (no solution after 0.091)')
ax[0].axhline(b[-1],color='gray',lw=0.8,ls=':',label='flat level %.4f'%b[-1])
ax[0].axhline(0.016331,color='#27ae60',lw=0.8,ls=':',label=r'$E\,b(X_1)=\sum\rho^{-1}/\sqrt{2}=0.0163$')
ax[0].set_xlabel('t'); ax[0].set_ylabel('b(t)'); ax[0].legend(fontsize=6.5,loc='upper right'); ax[0].set_title(r'Brownian boundary of $X_1$ in three zones',fontsize=10)
G=np.array(D['G']); ax[1].semilogx(G,np.array(D['S']),color='#1b4f72',lw=1.4,label='Monte Carlo, $1.5\\times10^5$ paths')
ax[1].semilogx(G,np.array(D['ex']),'--',color='#c0392b',lw=1.2,label='exact survival')
tt=np.logspace(-1,0,20); ax[1].semilogx(tt,0.01303/np.sqrt(tt),':',color='gray',label=r'$0.01303\,t^{-1/2}$')
ax[1].set_xlabel('t'); ax[1].legend(fontsize=7); ax[1].set_title('Survival across the exact boundary (KS 0.0019)',fontsize=10)
fig.tight_layout(); fig.savefig('fig_rs1.pdf'); print('ok')
