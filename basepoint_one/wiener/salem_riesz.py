import numpy as np, mpmath as mp, json, time
mp.mp.dps=20
out={}
# Salem kernel transform: K_sigma^(t) = Gamma(s) eta(s), s = sigma + i t ; Polya analogue: xi(sigma+it)
ts=np.linspace(0.05,60,2400)
sal={}
for sg in [0.5,0.6,0.75,1.0]:
    v=[float(abs(mp.gamma(mp.mpc(sg,t))*mp.altzeta(mp.mpc(sg,t)))) for t in ts]
    # normalise by the size of Gamma to remove the exponential decay
    w=[float(abs(mp.altzeta(mp.mpc(sg,t)))) for t in ts]
    sal[str(sg)]=dict(t=ts.tolist(),K=v,eta=w)
    loc=[ (ts[i],w[i]) for i in range(1,len(w)-1) if w[i]<w[i-1] and w[i]<w[i+1] and w[i]<0.2]
    print('sigma',sg,'local minima of |eta(sigma+it)| below 0.2 on (0,60]:',[(round(a,3),round(b,5)) for a,b in loc][:8])
out['salem']=sal
# Riesz function via Mobius sieve: R(x) = x * [ sum_{n<=N} mu(n) n^-2 (e^{-x/n^2}-1) + 6/pi^2 ]
N=10**7; t0=time.time()
mu=np.ones(N+1,dtype=np.int8); isp=np.ones(N+1,bool); isp[:2]=False
for p in range(2,N+1):
    if isp[p]:
        if p*p<=N: isp[p*p::p]=False
        mu[p::p]*=-1
        if p*p<=N: mu[p*p::p*p]=0
print('mobius sieve',time.time()-t0)
n=np.arange(1,N+1,dtype=np.float64); m=mu[1:].astype(np.float64); w=m/n**2
xs=np.geomspace(10,1e12,56); R=[]
for x in xs:
    val=x*(np.sum(w*np.expm1(-x/n**2))+6/np.pi**2)
    R.append(float(val))
R=np.array(R); print('R(x)/x^(1/4) range',(R/xs**0.25).min(),(R/xs**0.25).max())
for x,r in list(zip(xs,R))[::7]: print('x %.2e Riesz %.4f  /x^1/4 %.4f'%(x,r,r/x**0.25))
out['riesz']=dict(x=xs.tolist(),R=R.tolist())
json.dump(out,open('salem_riesz.json','w'))
