import numpy as np, json, mpmath as mp
from scipy.special import exp1
A=json.load(open('../basepoint_one/xi/alpha1.json')); R=json.load(open('../basepoint_one/xi/rs1.json')); I=json.load(open('../basepoint_one/xi/improve.json'))
g=np.concatenate([np.load(f) for f in ['g_1_400.npy','g_401_1000.npy','g_1001_1700.npy']])
th=np.array(A['grid']['th']); r=np.array(A['grid']['rho']); b=0.0230957089661
out={}
xs=[0.02,0.05,0.1,0.5,1,5,10,50,100]; t=np.array(R['t']); S=np.array(R['S'])
rows=[]
for x in xs:
    from scipy.integrate import quad
    rho=lambda q: np.interp(q,th,r)
    nu=sum(quad(lambda q: rho(q)*exp1(x*q*q),a,c,limit=400)[0] for a,c in [(0,0.05),(0.05,1),(1,10),(10,100)])
    Sx=float(np.interp(x,t,S)); asy=b/np.sqrt(np.pi*x)
    rows.append(dict(x=x,nu=float(nu),S=Sx,asy=asy,ratio=Sx/nu)); print('x %6.2f nu((x,inf)) %.6e  P(X1>x) %.6e  ratio %.4f  asympt %.6e'%(x,nu,Sx,Sx/nu,asy))
out['tail']=rows
# BDLP Levy density -W'(x) and its CM checks from zeros
xg=np.geomspace(0.002,0.2,200); dW=[float(np.sum(g*g*np.exp(-g*g*x))) for x in xg]
out['bdlp']=dict(x=xg.tolist(),mdW=dW)
# unimodality: X1 density (Talbot) and centre clock density
xt=np.array(I['talbot']['x']); ft=np.array(I['talbot']['f']); i=np.argmax(ft); print('X1 density mode near',xt[i],'value',ft[i],'monotone before/after:',np.all(np.diff(ft[:i+1])>=-1e-12),np.all(np.diff(ft[i:])<=1e-12))
CG=json.load(open('../data/clockgrid.json'))['rows']; HS=json.load(open('../data/hazard_spec.json'))['rows']
tc=np.array([q['t'] for q in CG]+[q[0] for q in HS]); fc=np.array([float(mp.mpf(q['f'])) for q in CG]+[float(mp.mpf(q[1])) for q in HS])
o=np.argsort(tc); tc,fc=tc[o],fc[o]; j=np.argmax(fc); print('centre clock mode near',tc[j],'f',fc[j],'unimodal:',np.all(np.diff(fc[:j+1])>=-1e-9),np.all(np.diff(fc[j:])<=1e-9))
out['uni']=dict(x1=dict(x=xt.tolist(),f=ft.tolist(),mode=float(xt[i])),centre=dict(t=tc.tolist(),f=fc.tolist(),mode=float(tc[j])))
json.dump(out,open('comp.json','w'))
