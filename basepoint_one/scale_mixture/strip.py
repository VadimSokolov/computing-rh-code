import numpy as np, mpmath as mp, json
mp.mp.dps=20
xi=lambda s: mp.mpf(0.5) if s==1 else s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
dl=lambda s: mp.diff(lambda z: mp.log(xi(z)), s)
g=np.concatenate([np.load(f) for f in ['g_1_400.npy','g_401_1000.npy','g_1001_1700.npy']])
out={}
# (1) densities of Y_alpha by Fourier inversion with the Cauchy corner removed analytically
U=np.concatenate([np.linspace(0,2,20001),np.linspace(2.0005,60,23000)])
dens={}
for al in [0.5,0.6,0.75,1.0]:
    xa=xi(mp.mpf(al)); b=float(mp.re(dl(mp.mpf(al)))) if al>0.5 else 0.0
    psi=np.array([float(xa/xi(mp.mpf(al)+mp.mpf(u))) for u in U]); rem=psi-np.exp(-b*U)
    ys=np.concatenate([np.linspace(0,1,51),np.geomspace(1.05,200,50)])
    f=[float(np.trapezoid(np.cos(U*y)*rem,U)/np.pi + b/(np.pi*(b*b+y*y))) for y in ys]
    dens[str(al)]=dict(y=ys.tolist(),f=f,b=b)
    idx=[np.argmin(abs(ys-q)) for q in [10,50,200]]
    print('alpha',al,'b',round(b,6),'f(0)',round(f[0],5),'y^2 f(y) at 10,50,200:',[round(ys[i]**2*f[i],6) for i in idx],'b/pi',round(b/np.pi,6))
out['dens']=dens
# (2) Thorin form and the product over verified zeros in the strip
weyl=lambda t: mp.log(t/(2*mp.pi))/(2*mp.pi)
rows=[]
for al in [0.6,0.75]:
    a=al-0.5; xa=xi(mp.mpf(al))
    th=np.linspace(0,200,8001)
    r=np.array([float(mp.re(dl(mp.mpc(al,t)))/mp.pi) if t>0 else float(mp.re(dl(mp.mpf(al)))/mp.pi) for t in th])
    print('alpha',al,'min rho on [0,200]',r.min(),'at',th[r.argmin()])
    for u in [1.0,5.0,20.0]:
        I=np.trapezoid(np.log1p(u*u/th[1:]**2)*r[1:],th[1:])+r[0]*float(mp.quad(lambda x: mp.log(1+u*u/x**2),[0,th[1]]))+float(mp.quad(lambda x: mp.log(1+u*u/x**2)*weyl(x),[200,mp.inf]))
        prodv=np.prod((a*a+g*g)/((a+u)**2+g*g))
        # zeros above 2198: smooth density, pair factor log((a+u)^2+t^2)-log(a^2+t^2)
        T=g[-1]+0.5
        tail=float(mp.quad(lambda t: (mp.log((a+u)**2+t*t)-mp.log(a*a+t*t))*weyl(t),[T,mp.inf]))
        ex=float(xa/xi(mp.mpf(al)+mp.mpf(u)))
        rows.append(dict(alpha=al,u=u,thorin=float(np.exp(-I)),product=float(prodv*np.exp(-tail)),exact=ex))
        print(' u',u,'Thorin %.7f  product over zeros %.7f  exact %.7f'%(np.exp(-I),prodv*np.exp(-tail),ex))
    out['rho%s'%al]=dict(min=float(r.min()))
out['rows']=rows
# (3) remainder bound from zeros above the verification height H
H=3e12
for th0 in [1e3,1e6,1e9]:
    bound=float(mp.log(H/(2*mp.pi))/(2*mp.pi)*0.5/mp.pi/(H-th0))
    print('theta',th0,'bound on |negative part| of rho from zeros above H:',bound)
out['H']=H
json.dump(out,open('strip.json','w'))
