# Second diagnostic companion of the reconstructed spacing.py (not the authors' code). The first order leaks of
# spacing.py agree with the archived data/spacing.json to about one float64 unit in the sums L(b), but only a few
# are identical. This script searches the natural float64 recipes (zero list, source of the cell boundaries,
# summation, final scaling) for one that reproduces the 79 archived leaks and coef exactly.
# Run: python3 variants2.py ARCHIVED.json DIR_WITH_30_DIGIT_NPY PART1.json   (writes variants2.json)
import json, math, os, sys
import numpy as np, mpmath as mp
from multiprocessing import Pool

FILES = ['g_1_400.npy', 'g_401_1000.npy', 'g_1001_1700.npy']
NZ, THETA = 1700, 200.0


def zz15(n):
    mp.mp.dps = 15
    return float(mp.zetazero(n).imag)


def sums(g, b):
    """L(c) = sum_j [1/(c-g_j)+1/(c+g_j)] at the boundaries c = b, by several summation recipes."""
    out = {}
    out['np.sum(1/(c-g)+1/(c+g))'] = np.array([np.sum(1 / (c - g) + 1 / (c + g)) for c in b])
    out['np.sum(1/(c-g))+np.sum(1/(c+g))'] = np.array([np.sum(1 / (c - g)) + np.sum(1 / (c + g)) for c in b])
    out['np.sum(2c/(c^2-g^2))'] = np.array([np.sum(2 * c / (c * c - g * g)) for c in b])
    out['2D sum over axis 0'] = np.sum(1 / (b[None, :] - g[:, None]) + 1 / (b[None, :] + g[:, None]), axis=0)
    out['2D sum over axis 1'] = np.sum(1 / (b[:, None] - g[None, :]) + 1 / (b[:, None] + g[None, :]), axis=1)
    out['python sum, sequential'] = np.array([sum(1 / (c - x) + 1 / (c + x) for x in g.tolist()) for c in b.tolist()])
    out['math.fsum'] = np.array([math.fsum(1 / (c - x) + 1 / (c + x) for x in g.tolist()) for c in b.tolist()])
    gl, bl = g.astype(np.longdouble), b.astype(np.longdouble)
    out['longdouble np.sum'] = np.array([np.sum(1 / (c - gl) + 1 / (c + gl)) for c in bl])
    mp.mp.dps = 30
    G = [mp.mpf(float(x)) for x in g]
    out['30 digits, float64 boundaries'] = [mp.fsum(1 / (mp.mpf(float(c)) - x) + 1 / (mp.mpf(float(c)) + x) for x in G)
                                            for c in b]
    return out


def finals(L):
    if isinstance(L, list):   # 30 digit sums
        return {'-(L_k-L_k-1)/pi at 30 digits': np.array([float(-(L[k + 1] - L[k]) / mp.pi) for k in range(79)])}
    L = np.asarray(L)
    f = {'-np.diff(L)/np.pi': -np.diff(L) / np.pi, 'np.diff(-L/np.pi)': np.diff(-L / np.pi),
         '-np.diff(L)*(1/np.pi)': -np.diff(L) * (1 / np.pi), '-np.diff(L/np.pi)': -np.diff(L / np.pi)}
    return {k: np.asarray(v, dtype=float) for k, v in f.items()}


if __name__ == '__main__':
    book = json.load(open(sys.argv[1])); bl = np.array(book['leaks'])
    g30 = np.concatenate([np.load(os.path.join(sys.argv[2], f)) for f in FILES])
    P1 = json.load(open(sys.argv[3]))
    ref = [mp.mpf(s) for s in P1['ref']] if isinstance(P1['ref'][0], str) else [mp.mpf(x) for x in P1['ref']]
    with Pool(int(os.environ.get('SLURM_CPUS_PER_TASK', os.cpu_count()))) as p:
        g15 = np.array(p.map(zz15, range(1, NZ + 1), chunksize=10))
    np.save('g15.npy', g15)
    near = np.array([np.longdouble(str(s)) for s in P1['ref']])
    rows = []
    for zname, g in [('dps30', g30), ('dps15', g15)]:
        gw = g[g < THETA]
        mp.mp.dps = 60
        B = {'float64 midpoints (g_k+g_k+1)/2': np.concatenate([[0.0], (gw[:-1] + gw[1:]) / 2, [THETA]]),
             'longdouble midpoints of part1 ref, as analyse.py': np.concatenate(
                 [[0.0], [float((near[k] + near[k + 1]) / 2) for k in range(78)], [THETA]]),
             'rounded exact midpoints of part1 ref': np.concatenate(
                 [[0.0], [float((ref[k] + ref[k + 1]) / 2) for k in range(78)], [THETA]])}
        for bname, b in B.items():
            for sname, L in sums(g, b).items():
                for fname, lk in finals(L).items():
                    same = int(np.sum(lk == bl)); coef = float(np.sum(np.abs(lk)))
                    d = float(np.max(np.abs(lk - bl) / np.abs(bl)))
                    rows.append(dict(zeros=zname, boundaries=bname, sum=sname, final=fname, identical=same,
                                     max_rel=d, coef_identical=coef == book['coef']))
    rows.sort(key=lambda r: (-r['identical'], r['max_rel']))
    json.dump(rows, open('variants2.json', 'w'), indent=1)
    print('%d recipes; best first: identical leaks of 79, max rel diff, coef identical' % len(rows))
    for r in rows[:25]:
        print('%2d %.2e %-5s | %s | %s | %s | %s' % (r['identical'], r['max_rel'], r['coef_identical'], r['zeros'],
                                                   r['boundaries'], r['sum'], r['final']))
    diffb = {z: int(np.sum(np.concatenate([[0.0], (g[g < THETA][:-1] + g[g < THETA][1:]) / 2, [THETA]])
                           != np.concatenate([[0.0], [float((near[k] + near[k + 1]) / 2) for k in range(78)], [THETA]])))
             for z, g in [('dps30', g30), ('dps15', g15)]}
    print('boundaries differing between float64 midpoints and the analyse.py cells:', diffb)
