import numpy as np, mpmath as mp, json, time
mp.mp.dps=20
xi=lambda s: mp.mpf(0.5) if s==1 else s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
dl=lambda s: mp.diff(lambda z: mp.log(xi(z)), s)
g=np.concatenate([np.load(f) for f in ['g_1_400.npy','g_401_1000.npy','g_1001_1700.npy']])
out={}
# (a) density of Y_alpha by Fourier inversion f(y)=(1/pi) int_0^inf cos(uy) xi(a)/xi(a+u) du
t0=time.time()
U=np.linspace(0,80,16001)
dens={}
for al in [0.5,0.75,1.0]:
    xa=xi(mp.mpf(al))
    psi=np.array([float(xa/xi(mp.mpf(al)+mp.mpf(u))) for u in U])
    ys=np.concatenate([np.linspace(0,1,101),np.geomspace(1.05,300,60)])
    f=[float(np.trapezoid(np.cos(U*y)*psi,U)/np.pi) for y in ys]
    ba=float(mp.re(dl(mp.mpf(al)))) if al>0.5 else 0.0
    dens[al]=dict(y=ys.tolist(),f=f,b=ba)
    print('alpha',al,'f(0)',f[0],'y^2 f(y) at y=10,100,300:',[round(ys[i]**2*f[i],6) for i in [np.argmin(abs(ys-10)),np.argmin(abs(ys-100)),len(ys)-1]],'b/pi',ba/np.pi,'time',time.time()-t0)
out['dens']=dens
# (b) Poisson trace as stable one half subordinate of the heat trace: K(y)=int y/sqrt(4 pi t^3) e^{-y^2/4t} W(t) dt
def Wz(t): return float(np.sum(np.exp(-g*g*t)))
Kc=[]
for y in [0.2,0.5,1.0]:
    k=lambda t: y/np.sqrt(4*np.pi*t**3)*np.exp(-y*y/(4*t))*Wz(t)
    from scipy.integrate import quad
    val=sum(quad(k,a,c,limit=400)[0] for a,c in [(1e-6,0.005),(0.005,0.05),(0.05,1),(1,50)])
    direct=float(np.sum(np.exp(-g*y)))
    Kc.append(dict(y=y,subordinated=val,direct=direct)); print('y',y,'int kernel*W',val,'sum e^{-gamma y}',direct)
out['Ksub']=Kc
# (c) rho_1 - b/pi >= 0 further: coarse scan theta in [200,600] plus local minima
th=np.linspace(0,600,24001); b=0.0230957089661
P=lambda x: 0.5/(np.pi*(0.25+x*x))
rz=np.array([np.sum(P(t-g)+P(t+g)) for t in th])   # zero form with beta=1/2, 1700 zeros (valid well below 2198)
off=float(np.pi*rz[0]); print('zero-sum rho1(0)*pi',off,'true b',b)
m=np.argmin(rz[th>1]); print('min over (1,600] of rho1 (zero form):',rz[th>1][m],'at',th[th>1][m],' b/pi',b/np.pi)
out['minrho']=dict(minval=float(rz[th>1].min()),at=float(th[th>1][m]),bpi=b/np.pi)
json.dump(out,open('improve.json','w'))
