import numpy as np, json, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
C=np.load('plotcurves.npz'); Z=C['zeros']
eps=[1.0,0.5,0.25,0.1,0.05,0.02,0.01,0.005]
tg=np.linspace(10,60,2500)
img=np.array([np.interp(tg,C['th_%s'%e],np.log10(C['F_%s'%e]/np.pi)) for e in eps])
fig,ax=plt.subplots(figsize=(8,3.6))
im=ax.pcolormesh(tg,np.arange(len(eps)),img,shading='auto',cmap='magma')
ax.set_yticks(range(len(eps))); ax.set_yticklabels([str(e) for e in eps]); ax.set_ylabel(r'resolution $\varepsilon$')

ax.set_xlabel(r'$\theta=\sqrt{z}$'); cb=fig.colorbar(im,ax=ax); cb.set_label(r'$\log_{10}\rho_\varepsilon(\theta)$')
ax.set_title('The Cauchy flow of Thorin densities sharpening onto the zeros')
fig.tight_layout(); fig.savefig('book/fig/flow.pdf'); plt.close(); print('flow ok')
