import numpy as np, mpmath as mp, json
mp.mp.dps=30
g=np.concatenate([np.load(f) for f in ['g_1_400.npy','g_401_1000.npy','g_1001_1700.npy']])
def xi(s): return s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
X0=xi(mp.mpf(0.5))
def Gex(w): return xi(mp.mpf(0.5)+mp.sqrt(w))/X0
rho=lambda x: (mp.re(mp.digamma(mp.mpf(1)/4+1j*x/2))/2-mp.log(mp.pi)/2)/mp.pi   # theta'(x)/pi
def Tstar(N):
    return mp.findroot(lambda T: mp.siegeltheta(T)/mp.pi+1-(N+0.5) if False else mp.siegeltheta(T)/mp.pi+1-N, (g[N-1]+g[N])/2)
rows=[]
for N in [1,20,79,300,1700]:
    T=float(Tstar(N)) if N<1700 else float(mp.findroot(lambda T: mp.siegeltheta(T)/mp.pi+1-1700,2197.9))
    for w in [1,100,1000,10000]:
        w=mp.mpf(w)
        part=mp.fprod([1+w/mp.mpf(float(x))**2 for x in g[:N]])
        tail=mp.quad(lambda x: mp.log(1+w/x**2)*rho(x),[T,10*T,100*T,mp.inf])
        ex=Gex(w)
        e1=abs(part/ex-1); e2=abs(part*mp.exp(tail)/ex-1)
        rows.append(dict(N=N,w=float(w),plain=float(e1),tail=float(e2)))
        print(N,float(w),mp.nstr(e1,3),mp.nstr(e2,3),flush=True)
json.dump(rows,open('kenttail.json','w'))
