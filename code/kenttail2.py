import numpy as np, mpmath as mp, json
mp.mp.dps=40
g=np.concatenate([np.load(f) for f in ['g_1_400.npy','g_401_1000.npy','g_1001_1700.npy']])
def xi(s): return s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
X0=xi(mp.mpf(0.5))
# exact power sums p_m = sum gamma^{-2m} from Taylor coefficients of f(s)=xi(1/2+s)/xi(1/2)
c=mp.taylor(lambda s: xi(mp.mpf(0.5)+s)/X0, 0, 6)
c2,c4,c6=c[2],c[4],c[6]
p1=c2; p2=c2**2-2*c4
print('p1 = sum 1/gamma^2 =',mp.nstr(p1,20),' p2 =',mp.nstr(p2,20))
rho=lambda x: (mp.re(mp.digamma(mp.mpf(1)/4+1j*x/2))/2-mp.log(mp.pi)/2)/mp.pi
def Gex(w): return xi(mp.mpf(0.5)+mp.sqrt(w))/X0
def tailint(f,T): return mp.quad(lambda x: f(x)*rho(x),[T,2*T,10*T,100*T,mp.inf])
rows=[]
for N in [79,300,1700]:
    gs=[mp.mpf(float(x)) for x in g[:N]]
    tau1=p1-mp.fsum([1/x**2 for x in gs]); tau2=p2-mp.fsum([1/x**4 for x in gs])
    T0=mp.findroot(lambda T: mp.siegeltheta(T)/mp.pi+1-N,(g[N-1]+g[N])/2 if N<1700 else 2197.9)
    # moment matched tail: density lam*rho on [T,inf) with int x^-2 = tau1 and int x^-4 = tau2
    def eqs(T,lam): return [lam*tailint(lambda x:x**-2,T)-tau1, lam*tailint(lambda x:x**-4,T)-tau2]
    T1,lam=mp.findroot(eqs,(T0,mp.mpf(1)))
    for w in [1,1e2,1e4,1e5,1e6,1e7]:
        w=mp.mpf(w); ex=Gex(w); part=mp.fprod([1+w/x**2 for x in gs])
        weyl=part*mp.exp(tailint(lambda x: mp.log(1+w/x**2),T0))
        mm=part*mp.exp(lam*tailint(lambda x: mp.log(1+w/x**2),T1))
        e=[abs(ex/part-1),abs(weyl/ex-1),abs(mm/ex-1)]
        rows.append(dict(N=N,w=float(w),plain=float(e[0]),weyl=float(e[1]),mm=float(e[2]),T0=float(T0),T1=float(T1),lam=float(lam)))
        print(N,float(w),*[mp.nstr(x,3) for x in e],' T0 %.3f T1 %.3f lam %.6f'%(T0,T1,lam),flush=True)
json.dump(dict(p1=str(p1),p2=str(p2),rows=rows),open('kenttail2.json','w'))
