import numpy as np, json, mpmath as mp, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':0.25,'axes.spines.top':False,'axes.spines.right':False})
B=json.load(open('base.json')); A=json.load(open('arith.json')); H=json.load(open('hit.json')); g=np.load('g12.npy')
Z=json.load(open('../xi/alpha1.json'))
th=np.array(B['th']); r=np.array(B['rho']); M=np.array(B['M'])
fig,ax=plt.subplots(1,2,figsize=(9.8,3.6))
ax[0].plot(th,r,color='#c0392b',lw=1.3,label=r'$L(s,\chi_{12})$')
zt=np.array(Z['grid']['th']); zr=np.array(Z['grid']['rho']); m=zt<=60
ax[0].plot(zt[m],zr[m],color='#1b4f72',lw=0.9,alpha=0.7,label=r'$\zeta$')
for z in g[g<60]: ax[0].axvline(z,color='#c0392b',lw=0.4,alpha=0.5)
ax[0].set_xlabel(r'$\theta$'); ax[0].legend(fontsize=8); ax[0].set_title(r'Thorin densities at the basepoint one',fontsize=10)
ax[1].plot(th,M,color='#c0392b',lw=1.5,label='Thorin mass'); ax[1].step(th,[np.sum(g<t) for t in th],where='post',color='k',lw=0.9,label='N(T)')
ax[1].set_xlabel('T'); ax[1].legend(fontsize=8); ax[1].set_title(r'Mass against the zero count of $L(s,\chi_{12})$',fontsize=10)
fig.tight_layout(); fig.savefig('fig_eta_density.pdf'); plt.close()
chi=[0,1,0,0,0,-1,0,-1,0,0,0,1]
fig,ax=plt.subplots(1,2,figsize=(9.8,3.6))
x=np.linspace(0.02,4.2,500); ax[0].plot(x,np.exp(-x)/(1-np.exp(-2*x)),color='#1b4f72',lw=1.5,label=r'gamma density $e^{-x}/(1-e^{-2x})$')
for p in [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59]:
    pk=p;k=1
    while pk<=60:
        v=chi[p%12]**k*np.log(p)/pk
        if v!=0: ax[0].vlines(np.log(pk),0,v,color='#27ae60' if v>0 else '#c0392b',lw=1.6)
        pk*=p;k+=1
ax[0].axhline(0,color='k',lw=0.6); ax[0].set_ylim(-0.6,1.5); ax[0].set_xlabel('x'); ax[0].legend(fontsize=7)
ax[0].set_title(r'Canonical measure of $\log\Lambda(1+s,\chi_{12})$: no pole, signed atoms',fontsize=9)
b1=B['b']; m=th<=40
Ad=np.array([float((mp.re(mp.digamma(mp.mpc(1,t)/2))-mp.digamma(mp.mpf(0.5)))/2) for t in th[m]])
P=np.pi*r[m]-b1-Ad
ax[1].plot(th[m],np.pi*r[m],color='k',lw=1.3,label=r'$\pi\rho_1$'); ax[1].plot(th[m],Ad,color='#1b4f72',lw=1.1,label=r'gamma $A(\theta)$')
ax[1].plot(th[m],P,color='#c0392b',lw=1.0,label=r'primes $P(\theta)$'); ax[1].axhline(b1,color='#27ae60',lw=1.1,label=r'drift $b=\sum\rho^{-1}$')
ax[1].set_xlabel(r'$\theta$'); ax[1].legend(fontsize=7); ax[1].set_title('Thorin density from the sine squared lemma',fontsize=10)
fig.tight_layout(); fig.savefig('fig_eta_levy.pdf'); plt.close()
S2=A['sin2']
fig,ax=plt.subplots(1,2,figsize=(9.8,3.6)); cols={1.0:'#1b4f72',0.9:'#27ae60',0.75:'#d68910',0.6:'#c0392b'}
for al,c in cols.items():
    rr=[q for q in S2 if q['alpha']==al and abs(q['theta']-3.8046276)<1e-6]
    X=[q['X'] for q in rr]; ax[0].loglog(X,[abs(q['err']) for q in rr],'o-',color=c,label=r'$\alpha=%s$'%al)
    ax[0].loglog(X,[abs(q['corr']) for q in rr],'s--',color=c,lw=0.8)
ax[0].set_xlabel('prime cutoff X'); ax[0].set_ylabel(r'error at $\theta=\gamma_1(\chi_{12})$'); ax[0].legend(fontsize=7); ax[0].set_title('Solid: primes alone. Dashed: after the zero correction',fontsize=9)
rr=[q for q in S2 if q['alpha']==0.5 and abs(q['theta']-3.8046276)<1e-6]
ax[1].semilogx([q['X'] for q in rr],[q['total'] for q in rr],'o-',color='#7d3c98',label='sum at the centre, first ordinate')
ax[1].set_xlabel('prime cutoff X'); ax[1].legend(fontsize=7); ax[1].set_title(r'At the centre the atom forms: growth like $\log X$',fontsize=10)
fig.tight_layout(); fig.savefig('fig_eta_below.pdf'); plt.close()
t=np.array(H['t']); b=np.array(H['b']); bt=np.array(H['bt'])
fig,ax=plt.subplots(1,2,figsize=(9.8,3.6))
for (a0,a1,col) in [(0.03,0.1,'#d5f5e3'),(0.1,0.4,'#fdebd0'),(0.4,1000,'#fadbd8')]: ax[0].axvspan(a0,a1,color=col,alpha=0.7,lw=0)
ax[0].semilogx(t,b,color='#c0392b',lw=2,label='exact boundary'); ax[0].semilogx(t,bt,'--',color='k',lw=1.1,label='tangent approximation (fails at 1.83)')
ax[0].axhline(b[-1],color='gray',ls=':',lw=0.8,label='flat level %.4f'%b[-1]); ax[0].set_ylim(-1.0,1.1)
ax[0].set_xlabel('t'); ax[0].legend(fontsize=7); ax[0].set_title(r'Brownian boundary of $X_1$ for $L(s,\chi_{12})$',fontsize=10)
mc=H['mc']; ax[1].semilogx(mc['G'],mc['S'],color='#c0392b',lw=1.4,label='Monte Carlo'); ax[1].semilogx(mc['G'],mc['ex'],'--',color='k',lw=1,label='exact survival')
ax[1].set_xlabel('t'); ax[1].legend(fontsize=8); ax[1].set_title('Survival across the exact boundary (KS 0.0014)',fontsize=10)
fig.tight_layout(); fig.savefig('fig_eta_hit.pdf'); plt.close(); print('ok')
