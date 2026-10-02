import numpy as np, json
g=np.concatenate([np.load(f) for f in ['g_1_400.npy','g_401_1000.npy','g_1001_1700.npy']])
out={}
# 1. zero density exponents e(sigma) with N(sigma,T) <= T^{e(sigma)+o(1)}
def ingham(s): return 3*(1-s)/(2-s)
def huxley(s): return 12*(1-s)/5
def gm(s): return 15*(1-s)/(3+5*s) if s>=0.7 else np.inf
def best(s): return min(ingham(s),huxley(s) if s>=0.75 else np.inf,gm(s))
rows=[]
for a in [0.55,0.6,0.65,0.7,0.75,0.8,0.9]:
    b=best(a); rows.append(dict(alpha=a,ingham=ingham(a),gm=gm(a) if a>=0.7 else None,best=b,saving=1-b,dh=2*(1-a)))
    print(a,round(ingham(a),4),gm(a) if a>=0.7 else '-',round(b,4),'saving',round(1-b,4),'DH',round(2*(1-a),3))
out['exp']=rows
# 2. synthetic experiment: add hypothetical off-line zeros (beta=0.75) at heights t_j with count ~ c t^e
th=np.linspace(0,2100,42001)
P=lambda a,x: a/(np.pi*(a*a+x*x))
def rho(alpha,off):
    r=np.zeros_like(th)
    for gg in g:
        r+=P(alpha-0.5,th-gg)+P(alpha-0.5,th+gg)
    for (b,t) in off:   # quadruple: b+it, 1-b+it (and conjugates via +theta)
        r+=P(alpha-b,th-t)+P(alpha-(1-b),th-t)+P(alpha-b,th+t)+P(alpha-(1-b),th+t)
    return r
e=0.6; c=1.0
tj=[]; k=1
while True:
    t=(k/c)**(1/e)*20   # hypothetical off-line zero heights: count grows like t^e
    if t>2000: break
    tj.append(t); k+=1
print('synthetic off line zeros',len(tj),[round(x,1) for x in tj[:6]])
off=[(0.75,t) for t in tj]
res=[]
for alpha in [0.6,0.7,0.8]:
    r=rho(alpha,off); r0=rho(alpha,[])
    dth=th[1]-th[0]
    for T in [250,500,1000,2000]:
        m=th<=T
        neg=float(np.sum(np.clip(-r[m],0,None))*dth); pos=float(np.sum(np.clip(r[m],0,None))*dth)
        cnt=sum(1 for t in tj if t<=2*T) if alpha<0.75 else 0
        res.append(dict(alpha=alpha,T=T,neg=neg,pos=pos,frac=neg/pos,bound=cnt,Npos=int(np.sum(g<=T))))
        print('alpha',alpha,'T',T,'Jordan negative mass %.3f'%neg,'positive %.1f'%pos,'fraction %.2e'%(neg/pos),'decomposition bound N_off(2T)=',cnt)
    out['curve_%s'%alpha]=dict(th=th[::10].tolist(),r=r[::10].tolist(),r0=r0[::10].tolist())
out['synthetic']=dict(tj=tj,res=res,e=e)
json.dump(out,open('almost.json','w'))
