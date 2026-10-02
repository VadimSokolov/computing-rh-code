import numpy as np, json
g=np.concatenate([np.load(f) for f in ['g_1_400.npy','g_401_1000.npy','g_1001_1700.npy']])
A=json.load(open('almost.json')); tj=A['synthetic']['tj']
M=json.load(open('more.json')) if False else {}
P=lambda a,x: a/(np.pi*(a*a+x*x))
out={}
g1=g[0]; tail=1.4797e-12
eps=4*tail/(2*np.pi)
out['window']=[]
for k in [1,2,3,4,6,8]:
    a=10.0**(-k); mids=np.concatenate([(g[:-1]+g[1:])/2,[0.0]])
    r=np.array([np.sum(P(a,t-g)+P(a,t+g)) for t in mids])
    lb=2*a/(np.pi*(a*a+g1*g1))
    out['window'].append(dict(alpha=0.5+a,min_mid=float(r.min()),lb=float(lb),eps=eps))
    print('alpha=1/2+1e-%d  min %.3e  lower bound %.3e  eps %.1e'%(k,r.min(),lb,eps))
print('threshold a* where lower bound = eps:',eps*np.pi*g1**2/2)
out['astar']=float(eps*np.pi*g1**2/2)
curve=[]
w=np.linspace(-6,6,12001); dw=w[1]-w[0]
for al in [0.52,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95]:
    bg=[]
    for t in tj:
        x=t+w; bg.append(np.sum(P(al-0.5,x[:,None]-g[None,:])+P(al-0.5,x[:,None]+g[None,:]),axis=1))
    for beta in [0.6,0.75,0.9]:
        if al>=beta: curve.append(dict(beta=beta,alpha=al,neg=0.0)); continue
        neg=0.0
        for t,b in zip(tj,bg):
            x=t+w; r=b.copy()
            for s in tj:
                r+=P(al-beta,x-s)+P(al-(1-beta),x-s)+P(al-beta,x+s)+P(al-(1-beta),x+s)
            neg+=np.sum(np.clip(-r,0,None))*dw
        curve.append(dict(beta=beta,alpha=al,neg=float(neg)))
        print('beta',beta,'alpha',al,'negative mass',round(neg,4),flush=True)
out['curve']=curve; out['tail']=tail; out['maxgap']=float(np.diff(g).max())
json.dump(out,open('more.json','w'))
