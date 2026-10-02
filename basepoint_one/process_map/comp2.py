import numpy as np, mpmath as mp, json, time
mp.mp.dps=20
xi=lambda s: s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s) if s!=1 else mp.mpf(0.5)
g=np.concatenate([np.load(f) for f in ['g_1_400.npy','g_401_1000.npy','g_1001_1700.npy']])
p1=0.0231049931154189707888; p2=0.0000371725992852696862
rng=np.random.default_rng(7); out={}
# 1. Riemann OU process: dY=-Y dt+dZ, Z = sum_k compound Poisson(rate 1, Exp(gamma_k^2)) + drift for missing zeros
K=len(g); rates=np.ones(K); d=p1-np.sum(1/g**2)
Tmax=60.0; t0=time.time()
njump=rng.poisson(K*Tmax); times=np.sort(rng.uniform(0,Tmax,njump)); which=rng.integers(0,K,njump)
sizes=rng.exponential(1/g[which]**2)
grid=np.arange(0,Tmax,0.01); Y=np.zeros(len(grid)); y=0.0; tprev=0.0; j=0
for i,tg in enumerate(grid):
    while j<njump and times[j]<=tg:
        dtt=times[j]-tprev; y=y*np.exp(-dtt)+d*(1-np.exp(-dtt)); y+=sizes[j]; tprev=times[j]; j+=1
    dtt=tg-tprev; Y[i]=y*np.exp(-dtt)+d*(1-np.exp(-dtt))
print('OU sim time',time.time()-t0,'jumps',njump)
st=Y[grid>5]; lt=[]
for s in [10,50,200]:
    emp=float(np.mean(np.exp(-s*st))); ex=float(xi(mp.mpf(0.5))/xi(0.5+mp.sqrt(s))); lt.append(dict(s=s,emp=emp,exact=ex)); print('OU stationary LT s',s,emp,ex)
ac=[]
for lag in [0.5,1.0,2.0]:
    L=int(lag/0.01); x=st-st.mean(); r=float(np.sum(x[:-L]*x[L:])/np.sum(x*x)); ac.append(dict(lag=lag,acf=r,exact=float(np.exp(-lag)))); print('acf',lag,r,np.exp(-lag))
out['ou']=dict(lt=lt,acf=ac,mean=float(st.mean()),path_t=grid[(grid>=5)&(grid<=15)][::2].tolist(),path_y=Y[(grid>=5)&(grid<=15)][::2].tolist())
print('OU mean',st.mean(),'p1',p1)
# 2. Wiener gamma with exact cell means of 1/z
A=json.load(open('alpha1.json')); th=np.array(A['grid']['th']); r=np.array(A['grid']['rho']); M=np.array(A['grid']['M'])
edges=np.concatenate([np.geomspace(1e-12,0.005,300),np.arange(0.01,M[-1],0.005)])
dts=np.diff(edges); inv=np.zeros(len(dts))
for i in range(len(dts)):
    sub=np.geomspace(edges[i],edges[i+1],41) if i<299 else np.linspace(edges[i],edges[i+1],41)
    zz=np.interp(sub,M,th)**2; inv[i]=np.trapezoid(1/zz,sub)/dts[i]
tailmean=float(mp.quad(lambda x: (mp.log(x/(2*mp.pi))/(2*mp.pi))/x**2,[100,mp.inf]))
N=200000; X=np.zeros(N)
for c in range(0,N,10000): X[c:c+10000]=(rng.gamma(dts,size=(10000,len(dts)))*inv).sum(1)+tailmean
wg=[]
for s in [0.5,1,10,100]:
    emp=np.mean(np.exp(-s*X)); se=np.std(np.exp(-s*X))/np.sqrt(N); ex=float(xi(mp.mpf(1))/xi(1+mp.sqrt(s))); wg.append(dict(s=s,emp=float(emp),se=float(se),exact=ex)); print('WG s',s,'%.5f +- %.5f exact %.5f'%(emp,se,ex))
out['wg']=wg
# 3. subordinated Brownian motion across alpha
bm=[]
for al in [0.6,0.75,1.0]:
    a=al-0.5; lam=g**2+a*a; K2=200
    tm=np.sum(1/lam[K2:])+(p1-np.sum(1/g**2)); tv=np.sum(1/lam[K2:]**2)+(p2-np.sum(1/g**4))
    n=300000; Tt=(rng.exponential(size=(n,K2))/lam[:K2]).sum(1)+rng.gamma(tm*tm/tv,tv/tm,n)
    c=np.sqrt(2)*a*Tt; Xa=Tt+c**2/rng.standard_normal(n)**2
    Yb=np.sqrt(Xa)*rng.standard_normal(n)
    ba=float(mp.diff(lambda z: mp.log(xi(z)),mp.mpf(al)))
    row=dict(alpha=al,b=ba,pred=np.sqrt(2)*ba/np.pi)
    for u in [1.0,3.0]:
        row['cf%g'%u]=float(np.mean(np.cos(u*Yb))); row['cfx%g'%u]=float(xi(mp.mpf(al))/xi(al+mp.mpf(u)/mp.sqrt(2)))
    for y in [10,100]: row['tail%g'%y]=float(np.mean(np.abs(Yb)>y)*y)
    bm.append(row); print(row)
out['bm']=bm
json.dump(out,open('comp2.json','w'))
