import numpy as np, json, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':0.25,'axes.spines.top':False,'axes.spines.right':False})
C=json.load(open('comp.json'))
fig,ax=plt.subplots(1,2,figsize=(9.8,3.6))
x=np.array(C['bdlp']['x']); d=np.array(C['bdlp']['mdW'])
ax[0].loglog(x,d,color='#1b4f72',lw=1.6,label=r"$-W'(x)=\sum\gamma^2e^{-\gamma^2x}$ (zeros)")
wx=[0.005,0.01,0.02,0.05,0.1]; wv=[165.1022811964,33.72924560562,3.741068142123,0.009166117863113,4.205189085735e-7]
ax[0].loglog(wx,wv,'o',color='#c0392b',label='explicit formula (primes, gamma, pole)')
ax[0].set_xlabel('x'); ax[0].set_title('Levy density of the driving process of the centre clock',fontsize=10); ax[0].legend(fontsize=7)
T=C['tail']; xs=np.array([q['x'] for q in T])
ax[1].loglog(xs,[q['S'] for q in T],'o-',color='#1b4f72',label=r'$P(X_1>x)$')
ax[1].loglog(xs,[q['nu'] for q in T],'s--',color='#c0392b',label=r'Levy tail $\nu((x,\infty))$')
ax[1].loglog(xs,[q['asy'] for q in T],':',color='gray',label=r'$0.01303\,x^{-1/2}$')
ax[1].set_xlabel('x'); ax[1].legend(fontsize=7); ax[1].set_title('Subexponential tail: the law follows its Levy tail',fontsize=10)
fig.tight_layout(); fig.savefig('fig_bdlp_tail.pdf'); plt.close()
U=C['uni']; fig,ax=plt.subplots(figsize=(6.5,3.4))
ax.semilogx(U['centre']['t'],U['centre']['f'],color='#1b4f72',lw=1.4,label='centre clock (mode %.3f)'%U['centre']['mode'])
xx=np.array(U['x1']['x']); ff=np.array(U['x1']['f']); m=ff>1e-8
ax.semilogx(xx[m],ff[m],color='#c0392b',lw=1.4,label=r'$X_1$ (mode %.3f)'%U['x1']['mode'])
ax.set_xlim(5e-3,1); ax.set_xlabel('x'); ax.set_ylabel('density'); ax.legend(fontsize=8); ax.set_title("Unimodal, as Yamazato's theorem requires",fontsize=10)
fig.tight_layout(); fig.savefig('fig_unimodal.pdf'); print('ok')
