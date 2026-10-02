import mpmath as mp, numpy as np, json
mp.mp.dps=25
chi=[0,1,0,0,0,-1,0,-1,0,0,0,1]
Lam=lambda s: (12/mp.pi)**(s/2)*mp.gamma(s/2)*mp.dirichlet(s,chi)
dlog=lambda s: mp.log(12/mp.pi)/2+mp.digamma(s/2)/2+mp.dirichlet(s,chi,derivative=1)/mp.dirichlet(s,chi)
g=np.load('g12.npy')
b=mp.re(dlog(mp.mpf(1)))
print('b = Lambda\'/Lambda(1) =',mp.nstr(b,15))
print('check sum over zeros 2*sum 0.5/(0.25+g^2) (297 zeros):',np.sum(1/(0.25+g**2)))
c=mp.taylor(lambda s: Lam(mp.mpf(0.5)+s)/Lam(mp.mpf(0.5)),0,4)
p1=c[2]; p2=c[2]**2-2*c[4]
print('p1 sum 1/g^2 =',mp.nstr(p1,15),' from zeros',np.sum(1/g**2),' p2',mp.nstr(p2,15),np.sum(1/g**4))
out=dict(b=float(b),p1=float(p1),p2=float(p2),tail=float(b/mp.sqrt(mp.pi)),Lam05=float(Lam(mp.mpf(0.5))),Lam1=float(Lam(mp.mpf(1))))
th=np.linspace(0,60,1201); r=[float(mp.re(dlog(mp.mpc(1,t)))/mp.pi) for t in th]
M=np.concatenate([[0],np.cumsum((np.array(r[1:])+np.array(r[:-1]))/2*np.diff(th))])
out.update(th=th.tolist(),rho=r,M=M.tolist())
print('rho(0)',r[0],'min',min(r),'M(60)',M[-1],'N(60)',int(np.sum(g<60)))
for t in [0,3.8046276,5,10,30,60]:
    ex=float(mp.re(dlog(mp.mpc(1,t)))/mp.pi); zs=np.sum(0.5/(np.pi*(0.25+(t-g)**2))+0.5/(np.pi*(0.25+(t+g)**2)))
    print('theta',t,'rho',ex,'zeros(297)',zs)
json.dump(out,open('base.json','w'))
