import numpy as np, mpmath as mp, json
mp.mp.dps=25
g=np.concatenate([np.load(f) for f in ['g_1_400.npy','g_401_1000.npy','g_1001_1700.npy']])
out={}
# Riesz sensitivity: amplitude of the term of each zero in R(x)/x^{1/4}
amp=[]
for t in g[:120]:
    rho=mp.mpc(0.5,t); a=abs(mp.gamma(1-rho/2)/(2*mp.zeta(rho,derivative=1)))
    amp.append(float(a))
amp=np.array(amp); fit=np.polyfit(g[:120],np.log(amp),1)
print('Riesz amplitude: slope of log amplitude per unit height',fit[0],' (pi/4 would be',np.pi/4,')')
print('amplitude of zero 1, 10, 50, 100:',amp[0],amp[9],amp[49],amp[99])
out['riesz_amp']=dict(g=g[:120].tolist(),amp=amp.tolist(),slope=float(fit[0]))
# margin statistics: unfolded zeros vs CUE eigenvalues
theta=lambda t: float(mp.siegeltheta(t))
u=np.array([theta(t)/np.pi+1 for t in g])   # unfolded ordinates, mean spacing 1
P=lambda k,x: k/(np.pi*(k*k+x*x))
def minima(pts,kap,lo,hi):
    mids=(pts[:-1]+pts[1:])/2; mids=mids[(mids>lo)&(mids<hi)]
    return np.array([np.sum(P(kap,m-pts)) for m in mids])
rng=np.random.default_rng(0)
def cue(N):
    Z=(rng.standard_normal((N,N))+1j*rng.standard_normal((N,N)))/np.sqrt(2)
    Q,R=np.linalg.qr(Z); d=np.diag(R); Q=Q*(d/abs(d))
    ang=np.sort(np.angle(np.linalg.eigvals(Q)))%(2*np.pi)
    return np.sort(ang)*N/(2*np.pi)
res={}
for kap in [0.1,0.25,0.5,1.0]:
    mz=minima(u,kap,50,u[-1]-50)
    mc=np.concatenate([minima(np.concatenate([c-400,c,c+400]),kap,50,350) for c in [cue(400) for _ in range(12)]])
    qs=[0.001,0.01,0.05,0.5]
    res[str(kap)]=dict(zeta=[float(np.quantile(mz,q)) for q in qs],cue=[float(np.quantile(mc,q)) for q in qs],nz=len(mz),nc=len(mc),zmin=float(mz.min()),cmin=float(mc.min()))
    print('kappa',kap,'zeta quantiles',np.round(res[str(kap)]['zeta'],4),'min',round(mz.min(),4),'| CUE',np.round(res[str(kap)]['cue'],4),'min',round(mc.min(),4))
out['margins']=res
json.dump(out,open('sens.json','w'))
