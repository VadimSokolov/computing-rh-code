# RECONSTRUCTED SCRIPT, not the authors' original: no script in the archive writes data/hankel.json, and this
# one is meant to reproduce that file (reconstructed in September 2026 under misc/repro/reconstruct/hankel).
# data/hankel.json holds Table tab:ch12:hankel of Chapter ch:heat (ch/ch12.tex) and the left panel of Figure
# fig:ch12:hankel, which code/plots3.py (section 7) draws from it. The book describes the numbers as the
# normalised leading Hankel minors "recomputed from the first 1700 zeros at sixty digits":
#   c_k(t) = (-1)^k W^(k)(t) = sum_j a_j^k exp(-a_j t), a_j = gamma_j^2, one zero of each pair,
#   even block det[c_{i+j}]_{i,j<N} / prod_{i<N} c_{2i}, odd block det[c_{i+j+1}]_{i,j<N} / prod_{i<N} c_{2i+1},
#   for N = 2, 3, 4, 5, and the log convexity ratio c_0 c_2 / c_1^2, at t = 0.01 and t = 0.05.
# Conventions follow code/arith.py and code/plots3.py: mpmath with mp.dps = 60, the 1700 float64 ordinates in
# g_1_400.npy, g_401_1000.npy and g_1001_1700.npy (misc/repro/zeros.py writes them), t handed to mpmath as a
# float. The output has the layout of the archived file: minors as mp.nstr(x, 4), the ratio as mp.nstr(x, 10).
# Run: python3 hankel.py   (writes hankel.json and prints every value to 15 digits)
import numpy as np, mpmath as mp, json
mp.mp.dps=60
g=np.concatenate([np.load(f) for f in ['g_1_400.npy','g_401_1000.npy','g_1001_1700.npy']])
a=[mp.mpf(float(x))**2 for x in g]
out={}
for t in [0.01,0.05]:
    T=mp.mpf(t); e=[mp.exp(-x*T) for x in a]
    c=[mp.fsum(x**k*y for x,y in zip(a,e)) for k in range(10)]
    ev=[]; od=[]
    for N in range(2,6):
        He=mp.matrix([[c[i+j] for j in range(N)] for i in range(N)])
        Ho=mp.matrix([[c[i+j+1] for j in range(N)] for i in range(N)])
        ev.append(mp.det(He)/mp.fprod(c[2*i] for i in range(N)))
        od.append(mp.det(Ho)/mp.fprod(c[2*i+1] for i in range(N)))
    r=c[0]*c[2]/c[1]**2
    out[str(t)]=dict(ratio=mp.nstr(r,10),even=[mp.nstr(v,4) for v in ev],odd=[mp.nstr(v,4) for v in od])
    print('t',t,'zeros',len(a),'ratio',mp.nstr(r,15),'even',[mp.nstr(v,15) for v in ev],'odd',[mp.nstr(v,15) for v in od],flush=True)
json.dump(out,open('hankel.json','w'))
print(json.dumps(out))
