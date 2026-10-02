import numpy as np, json, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':0.25,'axes.spines.top':False,'axes.spines.right':False})
rng=np.random.default_rng(4)
fig,ax=plt.subplots(2,2,figsize=(10,6.4))
# (a) BES(3) paths to level 1
dt=2e-4
for k,c in enumerate(['#1b4f72','#c0392b','#27ae60','#d68910']):
    x=np.zeros(3); r=[0.0]; t=[0.0]
    while r[-1]<1:
        x=x+np.sqrt(dt)*rng.standard_normal(3); r.append(np.linalg.norm(x)); t.append(t[-1]+dt)
    ax[0,0].plot(t,r,color=c,lw=0.7)
ax[0,0].axhline(1,color='k',lw=0.8,ls=':'); ax[0,0].set_xlabel('t'); ax[0,0].set_title(r'BES(3) paths to level 1: $E\,T^{(3)}=1/3$ (Williams)',fontsize=9)
# (b) Riemann OU path
O=json.load(open('ou.json')); y=np.array(O['path']); tt=np.arange(len(y))*0.02
ax[0,1].plot(tt,y,color='#7d3c98',lw=0.8); ax[0,1].axhline(O['p1'],color='k',lw=0.7,ls=':',label=r'mean $\sum\gamma^{-2}$')
ax[0,1].set_xlabel('t'); ax[0,1].legend(fontsize=7); ax[0,1].set_title(r"Riemann OU process: $dY=-Y\,dt+dZ$, jumps $\sim -W'$",fontsize=9)
# (c) boundaries across alpha
V=json.load(open('volterra.json')); R=json.load(open('rs1.json'))
ax[1,0].semilogx(V['t'],V['b'],color='#27ae60',lw=1.5,label=r'$\alpha=\frac{1}{2}$ (centre)')
for al,c in [('0.6','#c0392b'),('0.75','#d68910')]:
    B=json.load(open('bnd_%s.json'%al)); ax[1,0].semilogx(B['t'],B['b'],color=c,lw=1.5,label=r'$\alpha=%s$'%al)
ax[1,0].semilogx(R['t'],R['b'],color='#1b4f72',lw=1.5,label=r'$\alpha=1$')
ax[1,0].set_ylim(-1.2,0.6); ax[1,0].set_xlim(6e-3,100); ax[1,0].set_xlabel('t'); ax[1,0].legend(fontsize=7)
ax[1,0].set_title('Exact Brownian boundaries: linear at the centre, flat above it',fontsize=9)
# (d) Brownian paths against the alpha=1 boundary
tb=np.array(R['t']); bb=np.array(R['b']); grid=np.arange(0,0.3,1e-4); bnd=np.interp(grid,tb,bb,left=bb[0])
ax[1,1].plot(grid,bnd,color='#1b4f72',lw=2,label='boundary, $\\alpha=1$')
for k in range(6):
    w=np.concatenate([[0],np.cumsum(np.sqrt(1e-4)*rng.standard_normal(len(grid)-1))])
    hit=np.argmax(w>=bnd) if np.any(w>=bnd) else len(grid)-1
    ax[1,1].plot(grid[:hit+1],w[:hit+1],lw=0.6); ax[1,1].plot(grid[hit],w[hit],'o',ms=3,color='k')
ax[1,1].set_xlabel('t'); ax[1,1].legend(fontsize=7); ax[1,1].set_title(r'Brownian paths crossing the boundary of $X_1$',fontsize=9)
fig.tight_layout(); fig.savefig('fig_paths.pdf'); print('ok')
