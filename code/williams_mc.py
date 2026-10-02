import numpy as np, mpmath as mp, json, time
rng=np.random.default_rng(7); n=1_000_000; K=2000
w=2/(np.pi**2*np.arange(1,K+1)**2)
t0=time.time(); S=np.zeros(n)
for k in range(K): S+=w[k]*rng.gamma(2.0,1.0,n)
S+=4/(np.pi**2*K)      # tail replaced by its mean
T=0.5*np.log(np.pi*S/2); print('sample sec',time.time()-t0)
mp.mp.dps=20
def Xi(z): s=mp.mpf(0.5)+1j*mp.mpf(z); return mp.re(s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s))
rows=[]
for z in [0,2,5,8,10,14.134725,20,30]:
    y=np.exp(T/2)*np.cos(z*T); est=y.mean()/2; se=y.std()/np.sqrt(n)/2; ex=float(Xi(z))
    need=(y.std()/2/abs(ex))**2 if ex!=0 else float('inf')
    rows.append(dict(z=z,mc=est,se=se,exact=ex,n_for_rel_se_1=need)); print(z,est,se,ex,'%.2e'%need)
json.dump(rows,open('williams_mc.json','w'))
