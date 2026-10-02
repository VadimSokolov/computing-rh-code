import numpy as np, mpmath as mp, json
mp.mp.dps=18
xi=lambda s: mp.mpf(0.5) if s==1 else s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
dl=lambda s: mp.diff(lambda z: mp.log(xi(z)), s)
U=np.linspace(0,45,45001); S=json.load(open('strip.json')); dens={}
for al in [0.5,0.6,0.75,1.0]:
    xa=xi(mp.mpf(al)); b=float(mp.re(dl(mp.mpf(al)))) if al>0.5 else 0.0
    psi=np.array([float(xa/xi(mp.mpf(al)+mp.mpf(u))) for u in U])
    ys=np.concatenate([np.linspace(0,1,51),np.geomspace(1.05,50,40)])
    f=np.array([np.trapezoid(np.cos(U*y)*psi,U)/np.pi for y in ys])
    dens[str(al)]=dict(y=ys.tolist(),f=f.tolist(),b=b)
    idx=[np.argmin(abs(ys-q)) for q in [5,10,20,50]]
    print('alpha',al,'f(0) %.5f'%f[0],'y^2 f(y) at 5,10,20,50:',[round(ys[i]**2*f[i],6) for i in idx],'b/pi %.6f'%(b/np.pi),flush=True)
S['dens']=dens; json.dump(S,open('strip.json','w'))
