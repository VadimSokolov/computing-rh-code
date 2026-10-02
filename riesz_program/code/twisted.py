import numpy as np, mpmath as mp, json, time
mp.mp.dps=25
g=np.concatenate([np.load(f) for f in ['g_1_400.npy','g_401_1000.npy']])[:300]
out={}
T=50.0
amp=[float(abs(mp.gamma(1-mp.mpc(0.5,t)/2+1j*T)/(2*mp.zeta(mp.mpc(0.5,t),derivative=1)))) for t in g]
amp0=[float(abs(mp.gamma(1-mp.mpc(0.5,t)/2)/(2*mp.zeta(mp.mpc(0.5,t),derivative=1)))) for t in g]
k=int(np.argmax(amp)); print('twisted amplitude peaks at gamma =',g[k],'(2T =',2*T,')  value',amp[k],' untwisted amplitude there',amp0[k])
out['amp']=dict(g=g.tolist(),tw=amp,plain=amp0,T=T)
# twisted Riesz function from the Mobius series
N=10**7
mu=np.ones(N+1,dtype=np.int8); isp=np.ones(N+1,bool); isp[:2]=False
for p in range(2,N+1):
    if isp[p]:
        if p*p<=N: isp[p*p::p]=False
        mu[p::p]*=-1
        if p*p<=N: mu[p*p::p*p]=0
n=np.arange(1,N+1,dtype=np.float64); m=mu[1:].astype(np.float64); w=m/n**2
izT=complex(1/mp.zeta(2+2j*T))
xs=np.geomspace(1e4,1e10,40); vals=[]
for x in xs:
    ph=np.exp(1j*T*np.log(x/n**2))
    s=np.sum(w*(ph*np.exp(-x/n**2)-np.exp(1j*T*np.log(x))*n**(-2j*T)))+np.exp(1j*T*np.log(x))*izT
    vals.append(complex(x*s))
vals=np.array(vals)
# explicit formula with zeros up to 300 (both signs of gamma)
co=[(complex(mp.gamma(1-mp.mpc(0.5,t)/2+1j*T)/(2*mp.zeta(mp.mpc(0.5,t),derivative=1))),t) for t in g]
cn=[(complex(mp.gamma(1-mp.mpc(0.5,-t)/2+1j*T)/(2*mp.zeta(mp.mpc(0.5,-t),derivative=1))),-t) for t in g]
triv=[(complex(mp.gamma(1+kk+1j*T)/(2*mp.zeta(-2*kk,derivative=1))),kk) for kk in range(1,4)]
def pred(x): return sum(c*x**(0.25+0.5j*t) for c,t in co+cn)+sum(c*x**(-kk) for c,kk in triv)
pr=np.array([pred(x) for x in xs])
dev=np.abs(vals-pr)/xs**0.25
print('max |R_T - explicit|/x^(1/4):',dev.max(),'  signal |R_T|/x^(1/4) max',np.abs(vals/xs**0.25).max())
contrib_near=[sum(c*x**(0.25+0.5j*t) for c,t in co if abs(t-2*T)<15) for x in xs]
share=np.abs(np.array(contrib_near))/np.abs(pr)
print('median share of |R_T| from zeros with |gamma-2T|<15:',np.median(share))
out['series']=dict(x=xs.tolist(),re=vals.real.tolist(),im=vals.imag.tolist(),pre=pr.real.tolist(),pim=pr.imag.tolist(),dev=float(dev.max()),share=float(np.median(share)))
json.dump(out,open('twisted.json','w'))
