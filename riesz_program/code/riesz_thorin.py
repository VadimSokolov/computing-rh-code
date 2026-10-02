import numpy as np, mpmath as mp, json
mp.mp.dps=30
xi=lambda s: s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
g=np.concatenate([np.load(f) for f in ['g_1_400.npy','g_401_1000.npy','g_1001_1700.npy']])
p1=0.0231049931154189707888; p2=0.0000371725992852696862
out={}
# coefficients via the clock
rng=np.random.default_rng(2); n=400000; K=200
tm=p1-np.sum(1/g[:K]**2); tv=p2-np.sum(1/g[:K]**4)
T=(rng.exponential(size=(n,K))/g[:K]**2).sum(1)+rng.gamma(tm*tm/tv,tv/tm,n)
x05=float(xi(mp.mpf(0.5))); coef=[]
for k in range(1,7):
    ex=float(1/mp.zeta(2*k)); pref=2*k*(2*k-1)*float(mp.factorial(k-1))/(2*float(mp.pi)**k*x05)
    v=np.exp(-(2*k-0.5)**2*T); mc=pref*v.mean(); se=pref*v.std()/np.sqrt(n)
    coef.append(dict(k=k,exact=ex,mc=float(mc),se=float(se))); print('k',k,ex,mc,'+-',se)
out['coef']=coef
# Riesz function at small x: Mobius series vs averaged theta series over the clock
def Rmob(x,N=200000):
    import sympy
    return None
mu=np.ones(2000001,dtype=np.int8); isp=np.ones(2000001,bool); isp[:2]=False
for p in range(2,2000001):
    if isp[p]:
        if p*p<=2000000: isp[p*p::p]=False
        mu[p::p]*=-1
        if p*p<=2000000: mu[p*p::p*p]=0
nn=np.arange(1,2000001,dtype=float); w=mu[1:]/nn**2
small=[]
kk=np.arange(1,120)
for x in [0.5,1.0,2.0,5.0]:
    Rm=float(x*(np.sum(w*np.expm1(-x/nn**2))+6/np.pi**2))
    Rp=float(-sum((-mp.mpf(x))**k/(mp.factorial(k-1)*mp.zeta(2*k)) for k in range(1,200)))
    y=x/np.pi
    vals=np.zeros(n)
    for c in range(0,n,50000):
        Tc=T[c:c+50000][:,None]
        vals[c:c+50000]=np.sum(kk*(2*kk-1)*(-y)**kk*np.exp(-(2*kk-0.5)**2*Tc),axis=1)
    est=-vals.mean()/x05; se=vals.std()/np.sqrt(n)/x05
    small.append(dict(x=x,mobius=Rm,power=Rp,clock=float(est),se=float(se))); print('x',x,'Mobius',Rm,'power',Rp,'clock average',est,'+-',se)
out['small']=small
# cancellation sizes
xs=np.geomspace(10,1e8,40); canc=[]
for x in xs:
    raw=float(mp.log10(mp.mpf(x)**int(x)/mp.factorial(max(int(x)-1,0)))) if x<2e5 else float(x/np.log(10))
    k=np.arange(1,20000); y=x/np.pi
    row=dict(x=float(x),raw=raw)
    for Tv in [0.01,0.023,0.05]:
        lt=k*np.log(y)+np.log(k*(2*k-1))-(2*k-0.5)**2*Tv; row['T%s'%Tv]=float(lt.max()/np.log(10))
    canc.append(row)
out['canc']=canc
print([ (round(r['x']),round(r['raw'],1),round(r['T0.023'],1)) for r in canc[::8]])
json.dump(out,open('riesz_thorin.json','w'))
