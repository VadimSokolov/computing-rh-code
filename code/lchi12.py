import numpy as np, mpmath as mp, json, time
from flint import arb, acb, acb_poly, ctx
chi=[0,1,0,0,0,-1,0,-1,0,0,0,1]
def setupL(prec=320,hinv=200,umax=3.5):
    ctx.prec=prec; h=arb(1)/hinv; J=int(umax*hinv); pi=arb.pi(); K=[]
    for j in range(J+1):
        u=j*h; x=(2*u).exp(); s=arb(0)
        for n in range(1,400):
            c=chi[n%12]
            if c==0: continue
            t=c*(-pi*n*n*x/12).exp(); s+=t
            if abs(float(t.mid()))<2.0**(-prec-20) and n>12: break
        K.append(2*(u/2).exp()*s)
    cA=[h*K[0]/2]+[h*K[j] for j in range(1,J+1)]
    cB=[arb(0)]+[h*j*h*K[j] for j in range(1,J+1)]
    return dict(h=h,PA=acb_poly(cA),PB=acb_poly(cB))
def Lpair(S,s):
    s=acb(s); E=(s*S['h']).exp(); Ei=1/E
    return S['PA'](E)+S['PA'](Ei), S['PB'](E)-S['PB'](Ei)
S=setupL()
mp.mp.dps=30
def Lam(s): return (12/mp.pi)**(s/2)*mp.gamma(s/2)*mp.dirichlet(s,chi)
for t in [0,5,20]:
    A,_=Lpair(S,acb(0.1,t)); ref=Lam(mp.mpf('0.6')+1j*t)
    print('check',t,complex(A),complex(ref))
# zeros to height 100: sign changes of Lambda(1/2+it) (real)
th=np.arange(0,100.0001,0.02); v=[Lpair(S,acb(0,float(t)))[0].real for t in th]
roots=[]
for i in range(len(th)-1):
    if v[i]*v[i+1]<0:
        t=arb((th[i]+th[i+1])/2)
        for it in range(40):
            A,B=Lpair(S,acb(0,t)); dt=A.real/(-B.imag); t=(t-dt).mid()
            if abs(float(dt.mid()))<1e-40: break
        roots.append(t)
print('zeros below 100:',len(roots),[float(r.mid()) for r in roots[:6]])
# argument count and zero density at T=100
T=100.0
def argpath(pts):
    tot=0.0; prev=None
    for p in pts:
        a=Lpair(S,acb(*p))[0]
        if prev is not None: tot+=float((a/prev).arg())
        prev=a
    return tot
res=[]
for eps in [0.5,0.25,0.1,0.05,0.02]:
    M=argpath([(eps,float(t)) for t in np.arange(0,T+1e-9,eps/8)])/np.pi
    tau=argpath([(float(x),T) for x in np.linspace(eps,0,2001)])/np.pi
    res.append(dict(eps=eps,M=M,tau=tau,sum=M+tau)); print(eps,M,tau,M+tau,flush=True)
# density at eps 0.1 for plot
grid=np.arange(0,100.0001,0.0125); rho=[float((lambda ab:(ab[1]/ab[0]).real)(Lpair(S,acb(0.1,float(t)))))/np.pi for t in grid]
json.dump(dict(roots=[r.mid().str(30,radius=False) for r in roots],zd=res,grid=list(grid),rho=rho),open('lchi12.json','w'))
