import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from flint import acb, arb
import tilted
plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':0.25,'axes.spines.top':False,'axes.spines.right':False})
S=tilted.setup(0.01,256,0.0006,3.9); J=S['J']; h=S['h']; eta=S['eta']
nodes=[acb((k-J)*h,eta) for k in range(2*J+1)]; P0=list(S['P'].coeffs())
def Z(t,theta):
    s=acb(0,arb(str(theta))); tt=arb(str(t))/4; A=acb(0)
    for k,v in enumerate(nodes): A+=P0[k]*(tt*v*v+s*v).exp()
    return float(A.real*(arb.pi()*arb(str(theta))/4).exp())
th=np.linspace(7005.02,7005.14,61)
fig,ax=plt.subplots(figsize=(7,3.6))
for i,t in enumerate([0.0,-0.0004,-0.0007,-0.001]):
    v=np.array([Z(t,x) for x in th]); ax.plot(th,v/np.max(np.abs(v)),lw=1.6,label='t = %g'%t)
ax.axhline(0,color='k',lw=0.6); ax.set_xlabel(r'$\theta$'); ax.set_ylabel('normalised $H_t$')
ax.set_title('The Lehmer pair near 7005 under the de Bruijn Newman flow'); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig('book/fig/lehmerflow.pdf'); print('ok')
