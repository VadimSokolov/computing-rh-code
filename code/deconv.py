import numpy as np, json, time
from scipy.optimize import nnls
g=np.concatenate([np.load(f) for f in ['g_1_400.npy','g_401_1000.npy','g_1001_1700.npy']])
d=np.load('eps_0.1.npz'); th=d['theta']; F=d['F']/np.pi
eps=0.1
P=lambda a,x: a/(np.pi*(a*a+x*x))
m=(th>=15)&(th<=55); th=th[m][::2]; y=F[m][::2]
far=g[(g<5)|(g>65)]
bg=np.array([np.sum(P(eps,t-far)+P(eps,t+far)) for t in th])
xs=np.arange(5,65,0.02)
A=P(eps,th[:,None]-xs[None,:])+P(eps,th[:,None]+xs[None,:])
print(A.shape)
out={}
k=int(np.argmin(abs(g-32.935)))   # fifth zero region
for delta in [0.0,0.005,0.01,0.02,0.05,0.1]:
    pert=0.5*(P(eps+delta,th-g[k])+P(eps-delta,th-g[k]))-P(eps,th-g[k])
    t0=time.time(); w,res=nnls(A,y-bg+pert,maxiter=20000)
    big=xs[w>0.05]; mass=w.sum()
    out[str(delta)]=dict(res=res,mass=float(mass),sec=time.time()-t0)
    print(delta,'residual',res,'mass',mass,'sec',round(time.time()-t0,1),flush=True)
# atom recovery at delta=0
w,res=nnls(A,y-bg,maxiter=20000)
cl=[]; cur=[]
for x,wi in zip(xs,w):
    if wi>1e-3: cur.append((x,wi))
    elif cur: cl.append(cur); cur=[]
if cur: cl.append(cur)
cent=[(sum(x*wi for x,wi in c)/sum(wi for x,wi in c), sum(wi for x,wi in c)) for c in cl]
print('recovered clusters',[(round(a,4),round(b,4)) for a,b in cent])
out['clusters']=cent; out['true']=[float(x) for x in g[(g>5)&(g<65)]]
json.dump(out,open('deconv.json','w'))
