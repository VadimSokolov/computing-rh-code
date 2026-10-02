import numpy as np, mpmath as mp, json, time
mp.mp.dps=20
xi=lambda s: s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s) if s!=1 else mp.mpf(0.5)
out={}
rng=np.random.default_rng(1)
# (a) Biane Pitman Yor / Williams: S2 = sum 2/(pi^2 n^2) G_n, G_n ~ Gamma(2): two BES(3) hitting times of 1
n=400000; K=400
nn=np.arange(1,K+1); w=2/(np.pi**2*nn**2)
S2=np.zeros(n)
for c in range(0,n,50000): S2[c:c+50000]=(rng.gamma(2.0,size=(50000,K))*w).sum(1)
S2+=2*(1/3-np.sum(w))   # tail of the series replaced by its mean (E T = 1/3 for BES(3))
print('E S2 MC',S2.mean(),' exact 2/3')
rows=[]
for s in [0.5,1.0,1.5,2.0,3.0]:
    m=np.mean((np.pi*S2/2)**(s/2)); se=np.std((np.pi*S2/2)**(s/2))/np.sqrt(n)
    rows.append(dict(s=s,mc=float(m),se=float(se),exact=float(2*xi(mp.mpf(s)))))
    print('s',s,'E(pi S2/2)^(s/2) = %.5f +- %.5f   2 xi(s) = %.5f'%(m,se,2*xi(mp.mpf(s))))
out['bpy']=rows
# (c) Wiener gamma representation of X1: X1 = int dgamma_t / Minv(t), M(z)=U1((0,z]) = int_0^sqrt z rho1
A=json.load(open('alpha1.json')); th=np.array(A['grid']['th']); r=np.array(A['grid']['rho']); M=np.array(A['grid']['M'])
edges=np.concatenate([np.geomspace(1e-12,0.005,300),np.arange(0.01,M[-1],0.005)])
tt=(edges[1:]+edges[:-1])/2; dts=np.diff(edges)
thq=np.interp(tt,M,th); z=thq**2   # Minv(t)=theta(t)^2
b=0.0230957089661
# Thorin mass above theta=100: Weyl density (1/2pi)log(theta/2pi); tail enters through its mean
tailmean=float(mp.quad(lambda x: (mp.log(x/(2*mp.pi))/(2*mp.pi))/x**2,[100,mp.inf]))
N=200000; X=np.zeros(N)
for c in range(0,N,10000):
    X[c:c+10000]=(rng.gamma(dts,size=(10000,len(tt)))/z).sum(1)+tailmean
lt=[]
for s in [0.5,1,10,100]:
    emp=np.mean(np.exp(-s*X)); se=np.std(np.exp(-s*X))/np.sqrt(N); ex=float(xi(mp.mpf(1))/xi(1+mp.sqrt(s)))
    lt.append(dict(s=s,emp=float(emp),se=float(se),exact=ex)); print('Wiener gamma: s',s,'MC %.5f +- %.5f exact %.5f'%(emp,se,ex))
out['wg']=dict(lt=lt,tailmean=tailmean)
# (e) Brownian motion subordinated to X1: Y = sqrt(X1) N; char function xi(1)/xi(1+|u|/sqrt2); Cauchy type tail
Y=np.sqrt(X)*rng.standard_normal(N)
cf=[]
for u in [0.5,1,2,5]:
    emp=np.mean(np.cos(u*Y)); ex=float(xi(mp.mpf(1))/xi(1+mp.mpf(u)/mp.sqrt(2))); cf.append(dict(u=u,emp=float(emp),exact=ex)); print('char fn u',u,'MC %.5f exact %.5f'%(emp,ex))
tl=[]
for y in [1,10,100]:
    p=np.mean(np.abs(Y)>y); tl.append(dict(y=y,scaled=float(p*y))); print('y',y,'y P(|Y|>y)',p*y,' predicted',np.sqrt(2)*b/np.pi)
out['bm']=dict(cf=cf,tail=tl,pred=float(np.sqrt(2)*b/np.pi))
json.dump(out,open('comp.json','w'))
