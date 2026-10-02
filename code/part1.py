import numpy as np, mpmath as mp, json, time
from flint import acb, arb
from polya_flint import setup, xi_pair
S = setup(420, 200, 3)
Theta = 200
t0=time.time()
# line zeros: sign changes of Xi(theta) = xi(1/2 + i theta) on step 0.05, Newton refinement
th = np.arange(0, Theta+1e-9, 0.05)
Xi = [xi_pair(S, acb(0, float(t)))[0].real for t in th]
sg = [1 if x > 0 else -1 for x in Xi]   # check sign certainty via balls
unc = sum(1 for x in Xi if x.contains(0))
roots=[]
for k in range(len(th)-1):
    if sg[k] != sg[k+1]:
        t = arb((th[k]+th[k+1])/2)
        for it in range(40):
            A,B = xi_pair(S, acb(0, t))
            f = A.real; fp = -B.imag   # d/dtheta xi(1/2+i theta) = i xi', real part = -Im xi'
            dt = f/fp; t = (t - dt).mid()
            if abs(float(dt.mid())) < 1e-45: break
        roots.append(t)
print('sign changes', len(roots), 'uncertain signs', unc, time.time()-t0, flush=True)
# argument principle: s-1/2 along 1.5 -> 1.5 + i Theta -> i Theta
def argsum(pts):
    tot=0.0; prev=None; mx=0
    for p in pts:
        v = complex(xi_pair(S, acb(*p))[0])
        if prev is not None:
            d = np.angle(v/prev); mx=max(mx,abs(d)); tot += d
        prev = v
    return tot, mx
a1,m1 = argsum([(1.5, float(t)) for t in np.arange(0, Theta+1e-9, 0.01)])
a2,m2 = argsum([(float(x), Theta) for x in np.linspace(1.5, 0, 1501)])
N_arg = (a1+a2)/np.pi
print('N from argument principle', N_arg, 'max step', m1, m2, flush=True)
mp.mp.dps = 50
ref = [mp.zetazero(n).imag for n in range(1, len(roots)+2)]
diffs = [abs(mp.mpf(r.mid().str(50, radius=False)) - ref[i]) for i,r in enumerate(roots)]
out = dict(Theta=Theta, n_sign=len(roots), uncertain=unc, N_arg=N_arg, maxstep=[m1,m2],
           roots=[r.mid().str(45, radius=False) for r in roots],
           ref=[mp.nstr(g, 45) for g in ref], maxdiff=mp.nstr(max(diffs),5), next_zero=mp.nstr(ref[-1],20))
json.dump(out, open('part1.json','w'), indent=1)
print('max |polya root - zetazero|', out['maxdiff'], 'next', out['next_zero'], time.time()-t0)
