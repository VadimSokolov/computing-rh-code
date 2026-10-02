import time, json, numpy as np, mpmath as mp
from flint import acb, arb
import tilted
t0=time.time(); S=tilted.setup(0.01,256,0.0006,3.9); print('setup',time.time()-t0,'nodes',2*S['J']+1,flush=True)
mp.mp.dps=30
# Hardy style normalisation: Z(theta) ~ Xi(theta) * exp(pi theta/4)
th=np.round(np.arange(7004.95,7005.2001,0.005),4); out={'theta':list(th)}
t1=time.time(); Z=[]
for t in th:
    A,_=tilted.xi_pair_direct(S,acb(0,arb(str(t))))
    Z.append(float((A.real*(arb.pi()*arb(str(t))/4).exp())))
print('Xi grid',time.time()-t1, 'per point',(time.time()-t1)/len(th),flush=True)
out['Z']=Z
sc=[ (th[i],th[i+1]) for i in range(len(th)-1) if Z[i]*Z[i+1]<0]; print('sign changes',sc)
for eps in [0.05,0.02,0.01]:
    F=[]
    for t in th:
        A,B=tilted.xi_pair_direct(S,acb(arb(str(eps)),arb(str(t))))
        F.append(float((B/A).real)/np.pi)
    out['rho_%s'%eps]=F; print(eps,'max rho',max(F),'min',min(F),flush=True)
t2=time.time(); zs=[float(mp.siegelz(t)) for t in th[::10]]; print('mpmath siegelz per point',(time.time()-t2)/len(zs))
out['siegel_sample']=zs
json.dump(out,open('lehmer.json','w'))
