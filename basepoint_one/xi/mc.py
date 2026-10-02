import numpy as np, mpmath as mp, json
g=np.concatenate([np.load(f) for f in ['g_1_400.npy','g_401_1000.npy','g_1001_1700.npy']])
p1=0.0231049931154189707888; p2=0.0000371725992852696862
lam=g**2+0.25
K=200
rest=lam[K:]
# zeros beyond 1700: tail of sum 1/gamma^2 and 1/gamma^4 (shift by 1/4 is negligible there)
tm=np.sum(1/rest)+(p1-np.sum(1/g**2)); tv=np.sum(1/rest**2)+(p2-np.sum(1/g**4))
shape=tm*tm/tv; scale=tv/tm
rng=np.random.default_rng(7); n=1_000_000
Tt=np.zeros(n)
for c in range(0,n,100000):
    Tt[c:c+100000]=(rng.exponential(size=(100000,K))/lam[:K]).sum(1)+rng.gamma(shape,scale,100000)
Z=rng.standard_normal(n); H=Tt**2/(2*Z**2); X=Tt+H
xi=lambda s: s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
rows=[]
for s in [0.5,1,10,100,1000,10000]:
    emp=np.mean(np.exp(-s*X)); se=np.std(np.exp(-s*X))/np.sqrt(n)
    ex=float(mp.mpf(0.5)/xi(1+mp.sqrt(s)))
    rows.append(dict(s=s,emp=float(emp),se=float(se),exact=ex)); print('s',s,'MC %.6f +- %.6f'%(emp,se),'exact %.6f'%ex)
print('E Ttilde MC',Tt.mean(),' sum 1/rho',0.0230957089661)
tails=[]
for x in [1e2,1e4,1e6]:
    pr=np.mean(X>x); tails.append(dict(x=x,p=float(pr),scaled=float(pr*np.sqrt(x)))); print('x',x,'P(X>x) sqrt(x)',pr*np.sqrt(x),' predicted',0.0230957089661/np.sqrt(np.pi))
json.dump(dict(lt=rows,tails=tails,ET=float(Tt.mean())),open('mc.json','w'))
np.save('X_sample.npy',np.sort(X)[::50])
xs=np.logspace(-1.5,6,60); Xsort=np.sort(X)
surv=1-np.searchsorted(Xsort,xs)/n
json.dump(dict(xs=xs.tolist(),surv=surv.tolist()),open('surv.json','w'))
