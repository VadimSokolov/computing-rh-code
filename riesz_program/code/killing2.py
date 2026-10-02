import numpy as np, mpmath as mp, json
mp.mp.dps=30
g=np.concatenate([np.load(f) for f in ['g_1_400.npy','g_401_1000.npy','g_1001_1700.npy']])
p1=0.0231049931154189707888; p2=0.0000371725992852696862
b1=float(1+mp.euler/2-mp.log(4*mp.pi)/2)
out=json.load(open('killing.json'))
rng=np.random.default_rng(21); n=1000000; K=200
lam=g**2+0.25
tm=np.sum(1/lam[K:])+(p1-np.sum(1/g**2)); tv=np.sum(1/lam[K:]**2)+(p2-np.sum(1/g**4))
Tt=np.zeros(n)
for c in range(0,n,100000): Tt[c:c+100000]=(rng.exponential(size=(100000,K))/lam[:K]).sum(1)+rng.gamma(tm*tm/tv,tv/tm,100000)
X1=Tt+(Tt/np.sqrt(2))**2/rng.standard_normal(n)**2
# all killing probabilities (Rao Blackwellised) up to k=30 and the full Riesz series at small x
C=lambda k: mp.mpf(2*k*(2*k-1))*mp.factorial(k-1)/mp.pi**k
pk=[float(np.mean(np.exp(-(2*k-1)**2*X1))) for k in range(1,31)]
small=[]
for x in [0.5,1.0,2.0,3.0]:
    Rmc=float(-sum((-mp.mpf(x))**k/mp.factorial(k-1)*C(k)*pk[k-1] for k in range(1,31)))
    Rex=float(-sum((-mp.mpf(x))**k/(mp.factorial(k-1)*mp.zeta(2*k)) for k in range(1,80)))
    small.append(dict(x=x,mc=Rmc,exact=Rex)); print('x',x,'Riesz as alternating sum of killing probabilities',Rmc,' exact',Rex)
out['small']=small
# half Cauchy: moments of log(C1 C2) and zeta(2n)
Y=np.log(np.abs(rng.standard_cauchy(n)))+np.log(np.abs(rng.standard_cauchy(n)))
hc=[]
for nn in [1,2,3,4]:
    m=2*nn-2
    exact_mom=mp.taylor(lambda s: mp.sec(mp.pi*s/2)**2,0,m)[m]*mp.factorial(m)   # E[Y^m]
    z_from=mp.pi**2/(8*(2*nn-1)*(1-mp.mpf(2)**(-2*nn)))*exact_mom/mp.factorial(m)
    mc=np.mean(Y**m); se=np.std(Y**m)/np.sqrt(n)
    dual=exact_mom/mp.factorial(m)*C(nn)**-1/mp.zeta(2*nn)*mp.zeta(2*nn)  # placeholder
    prod=float(exact_mom/mp.factorial(m))*pk[nn-1]
    target=float(8*(2*nn-1)*(1-mp.mpf(2)**(-2*nn))/(mp.pi**2*C(nn)))
    hc.append(dict(n=nn,moment_exact=float(exact_mom),moment_mc=float(mc),se=float(se),zeta_from_moment=float(z_from),zeta=float(mp.zeta(2*nn)),dual_mc=prod,dual_target=target))
    print('n',nn,'E[(log C1C2)^%d] exact %.6f MC %.4f+-%.4f | zeta(2n) from moment %.10f zeta(2n) %.10f | duality product (MC) %.6f target %.6f'%(m,exact_mom,mc,se,z_from,mp.zeta(2*nn),prod,target))
out['halfcauchy']=hc
# Keating Snaith / Barnes G: CUE characteristic polynomial moments
def cue(N):
    Z=(rng.standard_normal((N,N))+1j*rng.standard_normal((N,N)))/np.sqrt(2)
    Q,R=np.linalg.qr(Z); d=np.diag(R); return Q*(d/abs(d))
ks=[]
for N in [10,20]:
    M=20000
    vals=np.array([abs(np.prod(1-np.linalg.eigvals(cue(N))))**2 for _ in range(M)])
    # Bourgade Hughes Nikeghbali Yor decomposition: Z_N = prod_{j=1}^N (1 + e^{i theta_j} sqrt(B_{1,j-1})), B_{1,0}=1
    th=rng.uniform(0,2*np.pi,(200000,N)); B=np.ones((200000,N))
    for j in range(2,N+1): B[:,j-1]=rng.beta(1,j-1,200000)
    zb=np.abs(np.prod(1+np.exp(1j*th)*np.sqrt(B),axis=1))**2
    for k in [1,2]:
        ex=float(mp.fprod([mp.gamma(j)*mp.gamma(j+2*k)/mp.gamma(j+k)**2 for j in range(1,N+1)]))
        mcv=np.mean(vals**k); se=np.std(vals**k)/np.sqrt(M); mb=np.mean(zb**k); seb=np.std(zb**k)/np.sqrt(len(zb))
        barnes=float(mp.barnesg(1+k)**2/mp.barnesg(1+2*k))
        ks.append(dict(N=N,k=k,exact=ex,cue=float(mcv),cue_se=float(se),beta=float(mb),beta_se=float(seb),ratio=ex/N**(k*k),barnes=barnes))
        print('N',N,'k',k,'E|Z_N|^{2k} exact %.4f  CUE MC %.4f+-%.4f  beta decomposition %.4f+-%.4f  | exact/N^{k^2} %.4f  G(1+k)^2/G(1+2k) %.4f'%(ex,mcv,se,mb,seb,ex/N**(k*k),barnes))
out['ks']=ks
json.dump(out,open('killing.json','w'))
