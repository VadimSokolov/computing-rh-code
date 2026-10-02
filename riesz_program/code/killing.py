import numpy as np, mpmath as mp, json, time
mp.mp.dps=25
xi=lambda s: mp.mpf(0.5) if s==1 else s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
g=np.concatenate([np.load(f) for f in ['g_1_400.npy','g_401_1000.npy','g_1001_1700.npy']])
p1=0.0231049931154189707888; p2=0.0000371725992852696862
b1=float(1+mp.euler/2-mp.log(4*mp.pi)/2)
out={}
# 1. simulate X1 = Ttilde + H_{Ttilde/sqrt2}, Ttilde = centre clock tilted by e^{-T/4} (RH based sampler)
rng=np.random.default_rng(12); n=1000000; K=200
lam=g**2+0.25
tm=np.sum(1/lam[K:])+(p1-np.sum(1/g**2)); tv=np.sum(1/lam[K:]**2)+(p2-np.sum(1/g**4))
Tt=np.zeros(n)
for c in range(0,n,100000): Tt[c:c+100000]=(rng.exponential(size=(100000,K))/lam[:K]).sum(1)+rng.gamma(tm*tm/tv,tv/tm,100000)
X1=Tt+(Tt/np.sqrt(2))**2/rng.standard_normal(n)**2
cb=b1/np.sqrt(2)
rows=[]
for k in range(1,7):
    rate=(2*k-1)**2
    E=rng.exponential(1/rate,n)
    pk=np.mean(X1<E); se=np.sqrt(pk*(1-pk)/n)
    pk_rb=np.mean(np.exp(-rate*X1)); se_rb=np.std(np.exp(-rate*X1))/np.sqrt(n)   # Rao Blackwellised
    Ck=2*k*(2*k-1)*float(mp.factorial(k-1))/float(mp.pi)**k
    target=float(1/mp.zeta(2*k))/Ck
    # Brownian drift factor: P(H_c < E) with c = b1/sqrt2, H_c = c^2/Z^2
    Hc=cb**2/rng.standard_normal(n)**2; E2=rng.exponential(1/rate,n)
    q=np.mean(Hc<E2); qse=np.sqrt(q*(1-q)/n)
    rows.append(dict(k=k,rate=rate,p=float(pk),se=float(se),p_rb=float(pk_rb),se_rb=float(se_rb),target=target,brown=float(q),brown_se=float(qse),stable=float(np.exp(-b1*(2*k-1))),C=Ck))
    print('k %d rate %3d  P(X1<E) %.6f+-%.6f  RB %.6f+-%.6f  target 1/(C_k zeta(2k)) %.6f | P(H_c<E) %.5f+-%.5f  e^{-b1(2k-1)} %.5f'%(k,rate,pk,se,pk_rb,se_rb,target,q,qse,np.exp(-b1*(2*k-1))))
out['kill']=rows
# Riesz function at small x as an alternating sum of killing probabilities
small=[]
for x in [0.5,1.0,2.0]:
    Rk=-sum((-x)**r['k']/float(mp.factorial(r['k']-1))*r['C']*r['p_rb'] for r in rows)
    Rex=float(-sum((-mp.mpf(x))**k/(mp.factorial(k-1)*mp.zeta(2*k)) for k in range(1,60)))
    # truncation at k<=6 from the exact coefficients, to separate Monte Carlo error from truncation
    R6=float(-sum((-mp.mpf(x))**k/(mp.factorial(k-1)*mp.zeta(2*k)) for k in range(1,7)))
    small.append(dict(x=x,mc=float(Rk),exact6=R6,exact=Rex)); print('x',x,'Riesz from killing probabilities (k<=6)',Rk,' exact k<=6',R6,' exact',Rex)
out['small']=small
json.dump(out,open('killing.json','w'))
