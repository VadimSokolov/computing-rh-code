import numpy as np, mpmath as mp, json, time
mp.mp.dps=25
chi=[0,1,0,0,0,-1,0,-1,0,0,0,1]
c12=lambda n: chi[int(n)%12]
Lam=lambda s: (12/mp.pi)**(s/2)*mp.gamma(s/2)*mp.dirichlet(s,chi)
dl=lambda s: mp.diff(lambda x: mp.log(Lam(x)), s)
g=np.load('g12.npy')
out={}
b=mp.log(12/mp.pi)/2+mp.digamma(mp.mpf(0.5))/2+mp.dirichlet(1,chi,derivative=1)/mp.dirichlet(1,chi)
print('b = Lambda\'/Lambda(1) =',mp.nstr(b,15),' check',mp.nstr(dl(mp.mpf(1)),15))
out['b']=float(b)
# power sums at the centre from Taylor coefficients
mp.mp.dps=40
L0=Lam(mp.mpf(0.5)); c=mp.taylor(lambda s: Lam(mp.mpf(0.5)+s)/L0,0,4)
p1=c[2]; p2=c[2]**2-2*c[4]; print('p1=sum gamma^-2',mp.nstr(p1,15),' partial 297',np.sum(g**-2.0),' p2',mp.nstr(p2,12))
out['p1']=float(p1); out['p2']=float(p2)
mp.mp.dps=25
# density at alpha=1, table and grid
def rho(al,th): return float(mp.re(dl(mp.mpc(al,th))))/np.pi
tab=[]
for th in [0,2,3.804627633,5,10,30]:
    d=rho(1,th) if th>0 else float(b)/np.pi
    zs=np.sum(0.5/(np.pi*(0.25+(th-g)**2))+0.5/(np.pi*(0.25+(th+g)**2)))
    tab.append((th,d,float(zs))); print('theta',th,'rho1',d,'zero sum (297 zeros)',zs)
out['table']=tab
th=np.linspace(0,60,1201); r=np.array([rho(1,t) if t>0 else float(b)/np.pi for t in th])
M=np.concatenate([[0],np.cumsum((r[1:]+r[:-1])/2*np.diff(th))])
out['grid']=dict(th=th.tolist(),rho=r.tolist(),M=M.tolist())
print('min rho1 on [0,60]',r.min(),'at',th[r.argmin()],' mass(60)',M[-1],' N(60)',int(np.sum(g<60)))
rows={}
for al in [0.75,0.6]:
    rr=[rho(al,t) if t>0 else float(dl(mp.mpf(al)))/np.pi for t in th]; rows[str(al)]=rr
    print('alpha',al,'Lambda\'/Lambda =',float(dl(mp.mpf(al))),'min rho',min(rr))
out['grid_alpha']=rows
out['balpha']={str(al):float(dl(mp.mpf(al))) for al in [0.6,0.75,0.9,1.0]}
# no-go growth
out['nogo']=[(t,float(abs(Lam(mp.mpf(1))/Lam(mp.mpc(1,t))))) for t in [10,20,40,80]]
print('nogo',out['nogo'])
json.dump(out,open('eta1.json','w'))
