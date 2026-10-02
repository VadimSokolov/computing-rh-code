import numpy as np, mpmath as mp, json
from scipy.integrate import quad
mp.mp.dps=20
xi=lambda s: s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
xiv=lambda s: mp.mpf(0.5) if s==1 else xi(s)
dl=lambda s: mp.diff(lambda z: mp.log(xi(z)), s)
g=np.concatenate([np.load(f) for f in ['g_1_400.npy','g_401_1000.npy','g_1001_1700.npy']])
p1=0.0231049931154189707888; p2=0.0000371725992852696862
out={}
# Thorin densities on a grid for alpha = 1 and 0.75 (direct from xi'/xi)
th=np.concatenate([np.linspace(0,60,2401),np.linspace(60.05,200,1400)])
rho={}
for al in [1.0,0.75]:
    rho[al]=np.array([float(mp.re(dl(mp.mpc(al,t)))/mp.pi) if t>0 else float(mp.re(dl(mp.mpf(al)))/mp.pi) for t in th])
    print('alpha',al,'rho computed; min on [0,200]',rho[al].min(),'at',th[rho[al].argmin()])
weyl=lambda t: mp.log(t/(2*mp.pi))/(2*mp.pi)
# (b) Thorin form of psi_alpha(u)=xi(alpha)/xi(alpha+|u|)
rows=[]
for al in [1.0,0.75]:
    r=rho[al]
    for u in [0.5,1,2,5,20]:
        I=np.trapezoid(np.log1p(u*u/th[1:]**2)*r[1:],th[1:])+r[0]*th[1]*0  # grid from th[1]
        # small theta piece [0,th1]: log(1+u^2/theta^2) integrable singularity
        h=th[1]; I+= r[0]*float(mp.quad(lambda x: mp.log(1+u*u/x**2),[0,h]))
        I+= float(mp.quad(lambda x: mp.log(1+u*u/x**2)*weyl(x),[200,mp.inf]))
        thor=np.exp(-I); ex=float(xiv(mp.mpf(al))/xiv(al+mp.mpf(u)))
        rows.append(dict(alpha=al,u=u,thorin=float(thor),exact=ex)); print('alpha',al,'u',u,'Thorin form %.7f  xi ratio %.7f'%(thor,ex))
out['thorin']=rows
# (c) Levy density kernel K_alpha(y)=int e^{-theta y} rho_alpha, and the centre Poisson trace sum e^{-gamma y}
Ky=[]
for y in [0.05,0.1,0.5,1,2,5,10]:
    K1=np.trapezoid(np.exp(-th*y)*rho[1.0],th)+float(mp.quad(lambda x: mp.e**(-x*y)*weyl(x),[200,mp.inf]))
    Kc=float(np.sum(np.exp(-g*y)))
    Ky.append(dict(y=y,K1=K1,K1_asym=0.0230957089661/(np.pi*y),Kcentre=Kc)); print('y',y,'K_1(y)',K1,' b/(pi y)',0.0230957089661/(np.pi*y),' centre sum e^{-gamma y}',Kc)
out['K']=Ky
# (a),(d) Monte Carlo of Y_alpha = sqrt(2 X_alpha) Z
rng=np.random.default_rng(5); mc=[]
for al in [0.5,0.75,1.0]:
    a=al-0.5; lam=g**2+a*a; K2=200; n=400000
    tm=np.sum(1/lam[K2:])+(p1-np.sum(1/g**2)); tv=np.sum(1/lam[K2:]**2)+(p2-np.sum(1/g**4))
    Tt=(rng.exponential(size=(n,K2))/lam[:K2]).sum(1)+rng.gamma(tm*tm/tv,tv/tm,n)
    X=Tt+(np.sqrt(2)*a*Tt)**2/rng.standard_normal(n)**2 if a>0 else Tt
    Y=np.sqrt(2*X)*rng.standard_normal(n)
    row=dict(alpha=al)
    for u in [1.0,5.0,20.0]:
        row['mc%g'%u]=float(np.mean(np.cos(u*Y))); row['se%g'%u]=float(np.std(np.cos(u*Y))/np.sqrt(n)); row['ex%g'%u]=float(xiv(mp.mpf(al))/xiv(al+mp.mpf(u)))
    ba=float(mp.re(dl(mp.mpf(al)))) if al>0.5 else 0.0
    row['var']=float(np.var(Y)) if al==0.5 else None
    row['tailpred']=2*ba/np.pi
    row['tail']=[float(np.mean(np.abs(Y)>yy)*yy) for yy in [10,100]]
    hist,edges=np.histogram(Y,bins=np.linspace(-1,1,201),density=True)
    row['hist']=hist.tolist(); row['edges']=edges.tolist()
    mc.append(row); print('alpha',al,{k:v for k,v in row.items() if k not in('hist','edges')})
out['mc']=mc
out['th']=th[::10].tolist(); out['rho1']=rho[1.0][::10].tolist(); out['rho075']=rho[0.75][::10].tolist()
json.dump(out,open('comp.json','w'))
