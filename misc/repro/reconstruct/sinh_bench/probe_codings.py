# Exploration for the reconstruction of data/sinh_bench.json, not an archive script.
# The paired coding of sinh_bench.py reproduces six of the ten archived values bit for bit (celltv and mass at
# eps=0.5, 0.05, 0.005) and misses the other four by at most 1.1e-12; probe_arctan.py showed that numpy's arctan
# dispatch is not the cause. This tries other codings of the same cumulative sums, which differ only in rounding.
# Usage: python3 probe_codings.py ARCHIVED.json
import numpy as np, json, sys, math, time
book = json.load(open(sys.argv[1]))
EPS = [0.5, 0.1, 0.05, 0.01, 0.005]; K = 400000; N = 50
b = np.concatenate([[0.0], np.arange(N+1)+0.5]); nat = np.array([0]+[1]*N); k = np.arange(1, K+1, dtype=float)
def fin(M): cm = np.diff(M); return float(np.sum(np.abs(cm-nat))), float(M[-1]-M[0])
codings = {
 'pair (sinh_bench.py)': lambda e: fin(np.array([np.sum(np.arctan((t-k)/e)+np.arctan((t+k)/e)) for t in b])/np.pi),
 'pair, M[-1] not M[-1]-M[0]': lambda e: (lambda M: (float(np.sum(np.abs(np.diff(M)-nat))), float(M[-1])))(np.array([np.sum(np.arctan((t-k)/e)+np.arctan((t+k)/e)) for t in b])/np.pi),
 'pair, /pi per term': lambda e: fin(np.array([np.sum(np.arctan((t-k)/e)/np.pi+np.arctan((t+k)/e)/np.pi) for t in b])),
 'pair, 2-D axis=1': lambda e: fin(np.sum(np.arctan((b[:, None]-k)/e)+np.arctan((b[:, None]+k)/e), axis=1)/np.pi),
 'pair, 2-D axis=0': lambda e: fin(np.sum(np.arctan((b[None, :]-k[:, None])/e)+np.arctan((b[None, :]+k[:, None])/e), axis=0)/np.pi),
 'pair, sequential sum': lambda e: fin(np.array([np.cumsum(np.arctan((t-k)/e)+np.arctan((t+k)/e))[-1] for t in b])/np.pi),
 'pair, math.fsum': lambda e: fin(np.array([math.fsum(np.arctan((t-k)/e)+np.arctan((t+k)/e)) for t in b])/np.pi),
 'atoms, reflections, math.fsum': lambda e: fin(np.array([math.fsum(np.concatenate([np.arctan((t-k)/e), np.arctan((t+k)/e)])) for t in b])/np.pi),
 'combined atan2(2te, e2+k2-t2)': lambda e: fin(np.array([np.sum(np.arctan2(2*t*e, e*e+k*k-t*t)) for t in b])/np.pi),
 'atom masses, two sums': lambda e: fin(np.array([np.sum(np.arctan((t-k)/e)+np.arctan(k/e)) + np.sum(np.arctan((t+k)/e)-np.arctan(k/e)) for t in b])/np.pi),
 'pair, arctan2(t-k, e)': lambda e: fin(np.array([np.sum(np.arctan2(t-k, e)+np.arctan2(t+k, e)) for t in b])/np.pi),
 'pair, x/e as x*(1/e)': lambda e: fin(np.array([np.sum(np.arctan((t-k)*(1/e))+np.arctan((t+k)*(1/e))) for t in b])/np.pi),
}
out = {}
for name, f in codings.items():
    t0 = time.time(); vals = []; bits = []
    for i, e in enumerate(EPS):
        tv, ms = f(e)
        for q, v in [('celltv', tv), ('mass', ms)]:
            a = book[i][q]; vals.append(v); bits.append('=' if v == a else '%.1e' % ((v-a)/a))
    nb = bits.count('=')
    out[name] = dict(values=vals, vs_archived=bits, bitwise_equal=nb)
    print('%-32s %2d of 10 bitwise  %s  (%.1f s)' % (name, nb, ' '.join(bits), time.time()-t0), flush=True)
json.dump(out, open('probe_codings.json', 'w'), indent=1)
