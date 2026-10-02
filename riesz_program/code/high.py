import numpy as np, mpmath as mp, json, time
mp.mp.dps=25
def dxi(s):  # xi'/xi via zeta and gamma factor
    return 1/s+1/(s-1)-mp.log(mp.pi)/2+mp.digamma(s/2)/2+mp.zeta(s,derivative=1)/mp.zeta(s)
P=lambda a,x: a/(np.pi*(a*a+x*x))
out={}
al=0.6
for gam0 in [300.0, 1.0e5]:
    t0=time.time()
    th=np.linspace(gam0-4,gam0+4,801)
    bg=np.array([float(mp.re(dxi(mp.mpc(al,t)))/mp.pi) for t in th])
    print('height',gam0,'background computed in %.1fs; min %.4f mean %.4f  (1/2pi)log(t/2pi)=%.4f'%(time.time()-t0,bg.min(),bg.mean(),np.log(gam0/(2*np.pi))/(2*np.pi)),flush=True)
    rows=[]
    for beta in [0.62,0.65,0.7,0.75,0.8,0.9,0.95]:
        w=beta-al
        r=bg+P(-w,th-gam0)+P(al-(1-beta),th-gam0)
        dth=th[1]-th[0]; neg=float(np.sum(np.clip(-r,0,None))*dth)
        rows.append(dict(beta=beta,min=float(r.min()),neg=neg,threshold=float(1/(np.pi*bg.min()))))
        print('   beta %.2f width %.2f  min rho %.3f  negative mass %.4f  (no negativity guaranteed if width >= 1/(pi min bg) = %.3f)'%(beta,w,r.min(),neg,1/(np.pi*bg.min())))
    out[str(gam0)]=dict(th=th.tolist(),bg=bg.tolist(),rows=rows)
json.dump(out,open('high.json','w'))
