import numpy as np, json, time
from scipy.stats import norm
from scipy.optimize import brentq
R=json.load(open('rs1.json')); t=np.array(R['t']); b=np.array(R['b']); F=np.array(R['F']); i0=R['i0']
f=np.gradient(F,t); t0=time.time()
bt=b.copy()
for k in range(i0+1,len(t)):
    tk=t[k]; dt=tk-t[k-1]; fk=max(f[k],1e-300)
    # backward Euler on t b' = b - sqrt(2 pi) t^{3/2} f e^{b^2/2t}, in log form for the exponential
    def H(x):
        return tk*(x-bt[k-1])/dt-(x-np.exp(np.log(np.sqrt(2*np.pi)*tk**1.5*fk)+x*x/(2*tk)))
    xs=bt[k-1]+np.linspace(-0.5,0.5,201)
    hv=np.array([H(x) for x in xs]); ok=[q for q in range(200) if np.isfinite(hv[q]) and np.isfinite(hv[q+1]) and hv[q]*hv[q+1]<=0]
    if not ok: print('tangent stop',k,tk); bt[k:]=np.nan; break
    q=min(ok,key=lambda q: abs(xs[q]-bt[k-1])); bt[k]=brentq(H,xs[q],xs[q+1])
print('time',time.time()-t0)
for q in [0.012,0.02,0.05,0.1,1,10,100]:
    i=np.argmin(abs(t-q)); print('t %.3f exact %.5f tangent %.5f'%(t[i],b[i],bt[i]))
print('tangent tail level vs sqrt(2pi)*0.006515',np.sqrt(2*np.pi)*0.0230957089661/(2*np.sqrt(np.pi)))
R['bt']=bt.tolist(); json.dump(R,open('rs1.json','w'))
