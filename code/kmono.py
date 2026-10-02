import numpy as np, json, matplotlib
from scipy.stats import gamma as G
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':0.25,'axes.spines.top':False,'axes.spines.right':False})
g=np.concatenate([np.load(f) for f in ['g_1_400.npy','g_401_1000.npy','g_1001_1700.npy']])
th=np.linspace(10,60,4000)
def Uk_theta(k):
    # atom a=gamma^2 spread to z=a k/G, G~Gamma(k); density in theta=sqrt(z)
    out=np.zeros_like(th); z=th**2
    for gm in g:
        a=gm*gm; x=a*k/z
        out+=np.exp(G.logpdf(x,k)+np.log(a*k/z**2)+np.log(2*th))
    return out
def pois(eps):
    return np.sum(eps/np.pi/(eps**2+(th[:,None]-g[None,:])**2)+eps/np.pi/(eps**2+(th[:,None]+g[None,:])**2),axis=1)
fig,ax=plt.subplots(2,1,figsize=(8,5.6),sharex=True)
res={}
for k,c in [(10**4,'#1b4f72'),(10**5,'#27ae60'),(10**6,'#c0392b')]:
    u=Uk_theta(k); ax[0].semilogy(th,u,color=c,lw=1.1,label='k = $10^{%d}$'%int(np.log10(k)))
    res[str(k)]=dict(min=float(u.min()))
for eps,c in [(30/(2*np.sqrt(1e4)),'#1b4f72'),(30/(2*np.sqrt(1e5)),'#27ae60'),(30/(2*np.sqrt(1e6)),'#c0392b')]:
    p=pois(eps); ax[1].semilogy(th,p,color=c,lw=1.1,label=r'$\varepsilon=%.3f$'%eps)
for a in ax:
    for z in g[(g>10)&(g<60)]: a.axvline(z,color='gray',lw=0.3,alpha=0.5)
    a.legend(fontsize=7,loc='lower right')
ax[0].set_title(r'Williamson measures $U_k$ (Gamma$(k)$ smoothing), density in $\theta=\sqrt{z}$',fontsize=10)
ax[1].set_title(r'Poisson smoothing $\rho_\varepsilon$ with $\varepsilon=30/(2\sqrt{k})$ (widths matched at $\theta=30$)',fontsize=10)
ax[1].set_xlabel(r'$\theta$'); ax[0].set_ylim(1e-6,50); ax[1].set_ylim(1e-3,50)
fig.tight_layout(); fig.savefig('book/fig/kmono.pdf')
# widths at first and 20th zero for k=1e5 (sd in theta)
for k in [1e4,1e5,1e6]:
    rel=np.sqrt(k)*np.sqrt(1/(k-2)) # approx
    res['relwidth_%d'%k]=float(1/(2*np.sqrt(k)))
print(res)
json.dump(res,open('kmono.json','w'))
