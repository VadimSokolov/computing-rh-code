import numpy as np, time, sys
from flint import acb
from polya_flint import setup, xi_pair
S = setup(420, 200, 3)
Theta = 200
for eps in [float(e) for e in sys.argv[1:]]:
    t0=time.time()
    dth = eps/(4 if eps<0.009 else 8); K = int(round(Theta/dth))
    F = np.zeros(K+1); A = np.zeros(K+1); prev=None; acc=0.0; mx=0
    for k in range(K+1):
        a,b = xi_pair(S, acb(eps, k*dth))
        F[k] = float((b/a).real)
        v = complex(a)
        if prev is not None:
            d = np.angle(v/prev); mx=max(mx,abs(d)); acc += d
        prev = v; A[k] = acc
    np.savez(f'eps_{eps}.npz', theta=np.arange(K+1)*dth, F=F, arg=A, maxstep=mx)
    print(eps, K, 'mass', A[-1]/np.pi, 'maxstep', mx, 'minF', F.min(), time.time()-t0, flush=True)
