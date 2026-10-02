import mpmath as mp, json, numpy as np
from clock import *
from flint import acb
out={}
# IG mixture no go: |L(u)| on the parabola u=(-q^2+2icq)/2 must be <= 1
ig=[]
for c in [0,5,10,20,24.93,30,50,100]:
    qs=np.linspace(0.05,60,600); first=None; mx=0
    for q in qs:
        u=mp.mpc(-q*q/2, c*q) if c>0 else mp.mpf(-q*q/2)
        v=abs(Lap(u))
        if v>1 and first is None: first=q
        mx=max(mx,float(v))
    ig.append(dict(c=c,qstar=first,max=mx)); print(c, first, mx, flush=True)
out['ig']=ig
# Turing plus Thorin: N(T) - M_eps(T) - top/pi = 2 nu(1/2+eps, T)
from polya_flint import xi_pair as xp
Theta=200; rows=[]
for eps in [1.0,0.5,0.25,0.1,0.05,0.02,0.01,0.005]:
    d=np.load(f'eps_{eps}.npz'); M=d['arg'][-1]/np.pi
    xs=np.linspace(eps,0,4001); prev=None; acc=0.0
    for x in xs:
        v=complex(xp(S, acb(float(x),Theta))[0])
        if prev is not None: acc+=np.angle(v/prev)
        prev=v
    top=acc/np.pi   # change of arg/pi going from 1/2+eps+iT to 1/2+iT
    rows.append(dict(eps=eps,M=M,top=top,sum=M+top)); print(eps, M, top, M+top, flush=True)
out['zd']=rows
json.dump(out,open('turing_ig.json','w'),indent=1)
