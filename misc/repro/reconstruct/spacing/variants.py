# Diagnostic companion of the reconstructed spacing.py (not the authors' code): which conventions reproduce the
# archived data/spacing.json, and how sensitive its numbers are to them. For ordinates rounded to float64 from
# mpmath.zetazero at 30, 20 and 15 digits, it recomputes mean, var, min and the leaks by several recipes and
# reports, for each, how many of the archived numbers it reproduces exactly and the largest relative difference.
# Run: python3 variants.py ARCHIVED.json DIR_WITH_30_DIGIT_NPY   (writes variants.json)
import json, os, sys
import numpy as np, mpmath as mp
from multiprocessing import Pool

FILES = ['g_1_400.npy', 'g_401_1000.npy', 'g_1001_1700.npy']
NZ, THETA = 1700, 200.0


def zz(arg):
    n, dps = arg
    mp.mp.dps = dps
    return float(mp.zetazero(n).imag)


def rel(a, b):
    a, b = np.atleast_1d(np.asarray(a, float)), np.atleast_1d(np.asarray(b, float))
    d = np.abs(a - b) / np.maximum(np.abs(a), np.abs(b))
    return int(np.sum(a == b)), float(d.max())


def spacings(g):
    out = {}
    mp.mp.dps = 30
    th30 = [mp.siegeltheta(t) for t in g]
    out['theta30, float(theta/pi), np.diff'] = np.diff(np.array([float(t / mp.pi) for t in th30]))
    out['theta30, float(theta)/np.pi, np.diff'] = np.diff(np.array([float(t) for t in th30]) / np.pi)
    out['theta30, exact differences'] = np.array([float((th30[k + 1] - th30[k]) / mp.pi) for k in range(NZ - 1)])
    mp.mp.dps = 15
    out['theta15, float(theta/pi), np.diff'] = np.diff(np.array([float(mp.siegeltheta(t) / mp.pi) for t in g]))
    mp.mp.dps = 30
    asym = lambda t: t / 2 * np.log(t / (2 * np.pi)) - t / 2 - np.pi / 8 + 1 / (48 * t) + 7 / (5760 * t ** 3)
    out['asymptotic theta to 1/t^3, float64'] = np.diff(asym(g) / np.pi)
    asym0 = lambda t: t / 2 * np.log(t / (2 * np.pi)) - t / 2 - np.pi / 8
    out['asymptotic theta without 1/t terms'] = np.diff(asym0(g) / np.pi)
    out['local unfolding by log(g_k/2pi)/2pi'] = np.diff(g) * np.log(g[:-1] / (2 * np.pi)) / (2 * np.pi)
    return out


def leaks(g):
    out = {}
    gw = g[g < THETA]
    b = np.concatenate([[0.0], (gw[:-1] + gw[1:]) / 2, [THETA]])
    L = np.array([np.sum(1 / (c - g) + 1 / (c + g)) for c in b])
    out['first order, differences of L(b), float64'] = -np.diff(L) / np.pi
    out['first order, four terms per cell, float64'] = np.array(
        [-np.sum(1 / (b[k + 1] - g) - 1 / (b[k] - g) + 1 / (b[k + 1] + g) - 1 / (b[k] + g)) / np.pi for k in range(79)])
    mp.mp.dps = 30
    G = [mp.mpf(float(t)) for t in g]
    B = [mp.mpf(0)] + [(G[k] + G[k + 1]) / 2 for k in range(78)] + [mp.mpf(THETA)]
    Lt = [mp.fsum(1 / (c - t) + 1 / (c + t) for t in G) for c in B]
    out['first order, 30 digits'] = np.array([float(-(Lt[k + 1] - Lt[k]) / mp.pi) for k in range(79)])
    for eps in [0.005, 0.001]:
        F = lambda y: np.arctan(y / eps) / np.pi
        out['finite eps=%g, arctan masses, float64' % eps] = np.array(
            [(np.sum(F(b[k + 1] - g) - F(b[k] - g) + F(b[k + 1] + g) - F(b[k] + g)) - 1) / eps for k in range(79)])
    b2 = b.copy(); b2[-1] = (g[78] + g[79]) / 2
    L2 = np.array([np.sum(1 / (c - g) + 1 / (c + g)) for c in b2])
    out['first order, last boundary at the midpoint 199.64'] = -np.diff(L2) / np.pi
    L3 = np.array([np.sum(1 / (c - g)) for c in b])
    out['first order, without the reflected ordinates'] = -np.diff(L3) / np.pi
    return out


if __name__ == '__main__':
    book = json.load(open(sys.argv[1]))
    g30 = np.concatenate([np.load(os.path.join(sys.argv[2], f)) for f in FILES])
    with Pool(int(os.environ.get('SLURM_CPUS_PER_TASK', os.cpu_count()))) as p:
        g20 = np.array(p.map(zz, [(n, 20) for n in range(1, NZ + 1)], chunksize=10))
        g15 = np.array(p.map(zz, [(n, 15) for n in range(1, NZ + 1)], chunksize=10))
    res = dict(zeros={'dps20 vs dps30': rel(g20, g30), 'dps15 vs dps30': rel(g15, g30)})
    print('ordinates identical in float64 / max rel diff:', res['zeros'], flush=True)
    stats = ['mean', 'var', 'min']
    for zname, g in [('dps30', g30), ('dps20', g20), ('dps15', g15)]:
        R = res[zname] = {}
        for name, s in spacings(g).items():
            v = dict(mean=np.mean(s), var=np.var(s), min=np.min(s))
            R['spacings: ' + name] = {k: rel(v[k], book[k]) for k in stats}
            R['spacings: ' + name]['var ddof=1'] = rel(np.var(s, ddof=1), book['var'])
        for name, lk in leaks(g).items():
            R['leaks: ' + name] = dict(leaks=rel(lk, book['leaks']), coef=rel(np.sum(np.abs(lk)), book['coef']))
        print('\nordinates', zname, '(numbers identical to the archived ones, max relative difference)')
        for k, v in R.items(): print('  %-58s %s' % (k, v), flush=True)
    json.dump(res, open('variants.json', 'w'), indent=1)
