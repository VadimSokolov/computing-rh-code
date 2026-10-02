# Check for the item part6_almost_dip: compares every vertex of the four curves drawn in the middle and right panels of
# the authors' fig_almost.pdf (read by fig_extract.py into fig_extract.json) with the density of almost.py recomputed at
# the same heights, and tests whether the vertices lie on the grid of step 0.005 from 284.
# Run on Hopper: bash misc/tools/hopper_run.sh -c 1 -t 10 compare_fig.py fig_extract.json g_1_400.npy g_401_1000.npy g_1001_1700.npy
import json
import numpy as np

g = np.concatenate([np.load(f) for f in ['g_1_400.npy', 'g_401_1000.npy', 'g_1001_1700.npy']])
tj, k = [], 1
while True:
    t = (k / 1.0) ** (1 / 0.6) * 20
    if t > 2000:
        break
    tj.append(t); k += 1
b = 0.75
P = lambda a, x: a / (np.pi * (a * a + x * x))


def r0(th, alpha):
    th = np.asarray(th)[:, None]; a = alpha - 0.5
    return (P(a, th - g[None, :]) + P(a, th + g[None, :])).sum(axis=1)


def r(th, alpha):
    th = np.asarray(th); s = r0(th, alpha)
    for t in tj:
        s = s + P(alpha - b, th - t) + P(alpha - (1 - b), th - t) + P(alpha - b, th + t) + P(alpha - (1 - b), th + t)
    return s


F = json.load(open('fig_extract.json'))
res = {}
for panel, alpha in (('panel1', 0.6), ('panel2', 0.7)):
    for c in F[panel]['curves']:
        V = np.array(c['vertices']); x, y = V[:, 0], V[:, 1]
        gray = abs(c['rgb'][0] - 0.502) < 0.01
        m = r0(x, alpha) if gray else r(x, alpha)
        d = np.abs(m - y); i = int(np.argmax(d))
        off = np.abs((x - 284) / 0.005 - np.round((x - 284) / 0.005)) * 0.005
        key = '%s_%s' % (alpha, 'line_only' if gray else 'with_hypothetical')
        res[key] = dict(vertices=int(len(x)), max_abs_diff=float(d.max()), at_theta=float(x[i]), median_abs_diff=float(np.median(d)),
                        max_distance_from_step_0005_grid=float(off.max()), drawn_min=float(y.min()), drawn_min_theta=float(x[np.argmin(y)]))
        print(key, res[key], flush=True)
json.dump(res, open('compare_fig.json', 'w'), indent=1)
