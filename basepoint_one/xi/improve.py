import numpy as np, mpmath as mp, json, time
mp.mp.dps=30
g=np.concatenate([np.load(f) for f in ['g_1_400.npy','g_401_1000.npy','g_1001_1700.npy']])
xi=lambda s: s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
out={}
# 1. Talbot inversion of phi_1 : density and CDF
phi=lambda s: mp.mpf(0.5)/xi(1+mp.sqrt(s))
xs=np.logspace(-2.3,3,40); dens=[]; t0=time.time()
for x in xs:
    f=mp.invertlaplace(phi,mp.mpf(x),method='talbot')
    dens.append(float(f))
print('talbot time',time.time()-t0)
out['talbot']=dict(x=xs.tolist(),f=dens)
# 2. density by subordination formula, averaging over Ttilde draws
p1=0.0231049931154189707888; p2=0.0000371725992852696862
lam=g**2+0.25; K=200; rest=lam[K:]
tm=np.sum(1/rest)+(p1-np.sum(1/g**2)); tv=np.sum(1/rest**2)+(p2-np.sum(1/g**4))
rng=np.random.default_rng(11); n=400000
Tt=(rng.exponential(size=(n,K))/lam[:K]).sum(1)+rng.gamma(tm*tm/tv,tv/tm,n)
sub=[]
for x in xs:
    d=x-Tt; m=d>0; c=Tt[m]/np.sqrt(2)
    val=np.sum(c/np.sqrt(2*np.pi*d[m]**3)*np.exp(-c*c/(2*d[m])))/n
    sub.append(float(val))
out['subdens']=sub
for x,a,b in list(zip(xs,dens,sub))[::5]: print('x %.4f talbot %.6e subordination %.6e'%(x,a,b))
# 3. Thorin densities and tail constants for several basepoints
th=np.linspace(0,60,1201); rows={}
for al in [0.6,0.75,1.0]:
    r=[float(mp.re(mp.diff(lambda s: mp.log(xi(s)), mp.mpc(al,t)))/mp.pi) if t>0 else float(mp.diff(lambda s: mp.log(xi(s)),mp.mpf(al))/mp.pi) for t in th]
    c=float(mp.diff(lambda s: mp.log(xi(s)),mp.mpf(al)))
    rows[str(al)]=dict(rho=r,xiprime=c,tail=c/np.sqrt(np.pi),minrho=min(r))
    print('alpha',al,'xi\'/xi(alpha)=%.8f tail const %.6f min rho %.5f'%(c,c/np.sqrt(np.pi),min(r)))
out['alphas']=dict(th=th.tolist(),rows=rows)
json.dump(out,open('improve.json','w'))
