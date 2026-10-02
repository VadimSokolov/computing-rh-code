import mpmath as mp, numpy as np, time
mp.mp.dps=18
chi=[0,1,0,0,0,-1,0,-1,0,0,0,1]
def Z(t):
    s=mp.mpf(0.5)+1j*mp.mpf(t)
    v=(12/mp.pi)**(s/2)*mp.gamma(s/2)*mp.dirichlet(s,chi)
    return mp.re(v)*mp.exp(mp.pi*t/4)   # remove the exponential decay
T=np.arange(0.05,340,0.06); t0=time.time()
v=np.array([float(Z(t)) for t in T]); print('scan',time.time()-t0)
idx=np.where(np.sign(v[:-1])*np.sign(v[1:])<0)[0]
zs=[float(mp.findroot(Z,(T[i],T[i+1]),solver='anderson')) for i in idx]
np.save('g12.npy',np.array(zs)); print(len(zs),zs[:5],zs[-1],'time',time.time()-t0)
