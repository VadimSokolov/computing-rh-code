import numpy as np, json, time
R=json.load(open('rs1.json')); t=np.array(R['t']); b=np.array(R['b']); S=np.array(R['S'])
rng=np.random.default_rng(5); n=150000
grid=np.concatenate([np.arange(0,0.05,2e-5),np.arange(0.05,1.0+1e-9,2e-4)])
bnd=np.interp(grid,t,b)
x=np.zeros(n); alive=np.ones(n,bool); Sv=[]; gaps={}; t0=time.time()
for k in range(1,len(grid)):
    dt=grid[k]-grid[k-1]; idx=np.nonzero(alive)[0]
    xn=x[idx]+np.sqrt(dt)*rng.standard_normal(idx.size)
    d0=bnd[k-1]-x[idx]; d1=bnd[k]-xn
    cr=(xn>=bnd[k])|(rng.random(idx.size)<np.exp(-2*np.clip(d0,0,None)*np.clip(d1,0,None)/dt))
    alive[idx[cr]]=False; x[idx]=xn; Sv.append(alive.mean())
    for q in [0.25,0.5,1.0]:
        if abs(grid[k]-q)<1e-9: gaps[q]=float(np.sum(bnd[k]-x[alive])/n)
G=grid[1:]; Sv=np.array(Sv); ex=np.interp(G,t,S)
print('sec',time.time()-t0,'KS',np.max(np.abs(Sv-ex)))
for q in [0.015,0.02,0.03,0.05,0.1,0.5,1.0]:
    i=np.argmin(abs(G-q)); print('t %.3f MC S %.5f exact S %.5f'%(G[i],Sv[i],ex[i]))
print('E[(b(t)-B_t) 1{alive}]:',gaps,' predicted sqrt(pi/2)*0.013030 =',np.sqrt(np.pi/2)*0.0230957089661/np.sqrt(np.pi))
json.dump(dict(G=G[::25].tolist(),S=Sv[::25].tolist(),ex=ex[::25].tolist(),ks=float(np.max(np.abs(Sv-ex))),gaps=gaps),open('rs1d.json','w'))
