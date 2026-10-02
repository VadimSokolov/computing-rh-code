# Third diagnostic companion of the reconstructed spacing.py (not the authors' code). variants2.py found that the
# sums A(c) = np.sum(1/(c-g)) and C(c) = np.sum(1/(c+g)) over the float64 ordinates, at the float64 midpoints,
# reproduce 61 of the 79 archived leaks exactly and the other 18 to one unit in the last place. This script tries
# the orderings of the final differences and of the sum giving coef.
# Run: python3 variants3.py ARCHIVED.json DIR_WITH_30_DIGIT_NPY   (writes variants3.json)
import json, math, os, sys
import numpy as np

FILES = ['g_1_400.npy', 'g_401_1000.npy', 'g_1001_1700.npy']
THETA, PI = 200.0, np.pi

if __name__ == '__main__':
    book = json.load(open(sys.argv[1])); bl = np.array(book['leaks'])
    g = np.concatenate([np.load(os.path.join(sys.argv[2], f)) for f in FILES])
    gw = g[g < THETA]
    b = np.concatenate([[0.0], (gw[:-1] + gw[1:]) / 2, [THETA]])
    A = np.array([np.sum(1 / (c - g)) for c in b]); C = np.array([np.sum(1 / (c + g)) for c in b])
    A0, A1, C0, C1 = A[:-1], A[1:], C[:-1], C[1:]
    L = A + C
    rec = {
        '-((A1+C1)-(A0+C0))/pi': -np.diff(L) / PI,
        '-((A1-A0)+(C1-C0))/pi': -((A1 - A0) + (C1 - C0)) / PI,
        '(((A0-A1)+C0)-C1)/pi': (((A0 - A1) + C0) - C1) / PI,
        '-(((A1-A0)+C1)-C0)/pi': -(((A1 - A0) + C1) - C0) / PI,
        '(((A0+C0)-A1)-C1)/pi': (((A0 + C0) - A1) - C1) / PI,
        '-(((A1+C1)-A0)-C0)/pi': -(((A1 + C1) - A0) - C0) / PI,
        '(A0-A1)/pi+(C0-C1)/pi': (A0 - A1) / PI + (C0 - C1) / PI,
        '-(A1-A0)/pi-(C1-C0)/pi': -(A1 - A0) / PI - (C1 - C0) / PI,
        '-(A1/pi+C1/pi)+(A0/pi+C0/pi)': -(A1 / PI + C1 / PI) + (A0 / PI + C0 / PI),
        '(A0+C0)/pi-(A1+C1)/pi': (A0 + C0) / PI - (A1 + C1) / PI,
        '-((A1-A0)+(C1-C0))*(1/pi)': -((A1 - A0) + (C1 - C0)) * (1 / PI),
        '-(1/pi)*((A1+C1)-(A0+C0))': -(1 / PI) * np.diff(L),
        '((A0+C0)-(A1+C1))/pi': ((A0 + C0) - (A1 + C1)) / PI,
        '-(np.diff(A)+np.diff(C))/pi': -(np.diff(A) + np.diff(C)) / PI,
        'np.diff(-A-C)/pi': np.diff(-A - C) / PI,
        '(np.diff(-A)+np.diff(-C))/pi': (np.diff(-A) + np.diff(-C)) / PI,
    }
    rows = []
    for name, lk in rec.items():
        lk = np.asarray(lk, float)
        cs = {'np.sum(np.abs)': float(np.sum(np.abs(lk))), 'sum(abs) sequential': float(sum(abs(v) for v in lk.tolist())),
              'math.fsum(abs)': math.fsum(abs(v) for v in lk.tolist())}
        rows.append(dict(recipe=name, identical=int(np.sum(lk == bl)),
                         max_rel=float(np.max(np.abs(lk - bl) / np.abs(bl))),
                         coef={k: v == book['coef'] for k, v in cs.items()}))
    rows.sort(key=lambda r: (-r['identical'], r['max_rel']))
    json.dump(rows, open('variants3.json', 'w'), indent=1)
    for r in rows: print('%2d %.2e %-36s coef identical: %s' % (r['identical'], r['max_rel'], r['recipe'], r['coef']))
