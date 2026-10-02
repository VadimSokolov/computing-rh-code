import numpy as np, mpmath as mp, json, time
mp.mp.dps=20
xi=lambda s: s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
g=np.concatenate([np.load(f) for f in ['g_1_400.npy','g_401_1000.npy','g_1001_1700.npy']])
p1=0.0231049931154189707888; K=len(g); d=p1-np.sum(1/g**2)
rng=np.random.default_rng(11); Tmax=4000.0; D=0.02; t0=time.time()
nj=rng.poisson(K*Tmax); tau=rng.uniform(0,Tmax,nj); J=rng.exponential(1/g[rng.integers(0,K,nj)]**2)
cell=np.minimum((tau//D).astype(np.int64),int(Tmax/D)-1); w=J*np.exp(-((cell+1)*D-tau))
inc=np.bincount(cell,weights=w,minlength=int(Tmax/D))
e=np.exp(-D); Y=np.zeros(len(inc)+1)
for i in range(len(inc)): Y[i+1]=Y[i]*e+d*(1-e)+inc[i]
Y=Y[int(20/D):]; print('sim',time.time()-t0,'jumps',nj)
out={'mean':float(Y.mean()),'p1':p1}
lt=[]
for s in [10,50,200]:
    ex=float(xi(mp.mpf(0.5))/xi(0.5+mp.sqrt(s))); emp=float(np.mean(np.exp(-s*Y)))
    # batch means standard error
    B=np.array_split(np.exp(-s*Y),100); se=float(np.std([b.mean() for b in B])/10)
    lt.append(dict(s=s,emp=emp,se=se,exact=ex)); print('s',s,emp,'+-',se,'exact',ex)
ac=[]
x=Y-Y.mean()
for lag in [0.5,1,2]:
    L=int(lag/D); r=float(np.sum(x[:-L]*x[L:])/np.sum(x*x)); ac.append(dict(lag=lag,acf=r,exact=float(np.exp(-lag)))); print('acf',lag,r,np.exp(-lag))
out['lt']=lt; out['acf']=ac; out['path']=Y[:int(10/D)].tolist(); print('mean',Y.mean(),p1)
json.dump(out,open('ou.json','w'))
