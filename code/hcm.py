import mpmath as mp, json
from flint import acb, arb
from polya_flint import setup, xi_pair
S = setup(420, 200, 3)
mp.mp.dps = 60
def xi_half(x):   # xi(1/2 + x), x real, from Polya
    a,_ = xi_pair(S, acb(arb(mp.nstr(x, 70))))
    return mp.mpf(a.real.mid().str(70, radius=False))
out = {}
for alpha in [0.5, 1.0]:
    sh = mp.mpf(alpha) - mp.mpf(0.5)
    base = xi_half(sh)
    G = lambda x: xi_half(sh + mp.sqrt(x))/base
    def H(w, u):
        v = (w + mp.sqrt(w*w-4))/2
        return 1/(G(u*v)*G(u/v))
    rows = []
    for u in [1, 25, 400]:
        for w in [2.05, 3, 5, 10, 20, 50]:
            h0 = H(mp.mpf(w), u)
            ds = [float((-1)**k*mp.diff(lambda x: H(x,u), mp.mpf(w), k, h=mp.mpf('1e-8'))/h0) for k in range(1,7)]
            rows.append(dict(u=u, w=w, H=float(h0), d=ds))
            print(alpha, u, w, float(h0), ['%.3e'%d for d in ds], flush=True)
    out[str(alpha)] = rows
json.dump(out, open('hcm.json','w'), indent=1)
