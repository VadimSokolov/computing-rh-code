import time, numpy as np, mpmath as mp, json
from flint import acb, arb
import tilted
S=tilted.setup(0.1,200,0.004,4.2)
T=1000.0
def xi(eps,t):
    A,_=tilted.xi_pair_direct(S,acb(arb(str(eps)),arb(repr(float(t))))); return A
t0=time.time()
# argument along Re s = 1 (eps = 1/2) from 0 to T; step 0.1
acc=0.0; prev=xi(0.5,0.0); mx=0
for t in np.arange(0.1,T+1e-9,0.1):
    v=xi(0.5,t); d=float((v/prev).arg()); mx=max(mx,abs(d)); acc+=d; prev=v
M=acc/np.pi
# horizontal segment 1 -> 1/2 at height T
acc2=0.0; prev=xi(0.5,T)
for x in np.linspace(0.5,0,2001)[1:]:
    v=xi(float(x),T); acc2+=float((v/prev).arg()); prev=v
tau=acc2/np.pi; tp=time.time()-t0
print('M_1/2(1000)',M,'tau',tau,'N',M+tau,'max step',mx,'sec',tp,flush=True)
mp.mp.dps=15; t1=time.time(); n=mp.nzeros(T); tm=time.time()-t1
print('mpmath nzeros(1000) (Riemann Siegel plus Turing):',n,'sec',tm)
json.dump(dict(M=M,tau=tau,N=M+tau,maxstep=mx,sec=tp,mp_n=int(n),mp_sec=tm),open('turing1000.json','w'))
