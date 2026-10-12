# maxexp_flow.py: the trajectories of the real zeros of the de Bruijn Newman flow H^Delta_t of Section sec:ch16:dbn (ch/ch16.tex,
# Section 8.5.5), for Figure fig:maxexp:flow, and an independent recomputation of the collisions of Table tab:ch16:coll. Added
# for the book in October 2026; the first author's note on maxima of exponentials (misc/reviews/nick-maxexp/, Figure 4 there)
# plotted the zeros below height 33.
#   H^Delta_t(theta) = int_0^4 e^{t u^2} K_Delta(u) cos(theta u) du,  K_Delta(u) = e^{6u} Delta(i e^u) = e^{6u - 2 pi e^u} F(e^u),
# with F(y) = prod (1 - e^{-2 pi n y})^24 for y >= 1; K_Delta < 10^-130 beyond u = 4. Gauss Legendre on [0, 4], 8 panels of
# 96 nodes, kernel and weights at 70 digits, as in dbn_delta.py. On the grid theta = 0.5, 0.52, ..., 52 the cosines are kept
# in fixed point (integers scaled by 2^220, from the recurrence cos((k+1)hu) = 2 cos(hu) cos(khu) - cos((k-1)hu)), so that
# H^Delta_t on the whole grid costs one integer matrix product per t. For t = 0.5, 0.49, ..., -1.1 the zeros are the sign
# changes on the grid, located by cubic inverse interpolation. Between consecutive t the count of zeros below 52 is compared
# (a step with more than one event is cut into eight, up to three times): a drop by one is the highest zero leaving through
# 52, a drop by two a collision of the closest neighbouring pair, which is then
# located as the time at which the critical value of H^Delta_t between the pair changes sign (the critical point by Newton's
# method on dH/dtheta = 0, the time by bisection, both with direct evaluation at 60 digits). Zeros keep their order on the real
# line, so each colliding pair is labelled by its zeros at t = 0.
# Run from any folder: python3 maxexp_flow.py; writes maxexp_flow.json there.
import json, time
from operator import mul
import mpmath as mp

t0 = time.time()
mp.mp.dps = 70
out = {}
def put(k, v):
    out[k] = v
    print(k, v, round(time.time() - t0, 1), flush=True)
TINY = mp.mpf(10) ** -80
def K(u):
    y = mp.exp(u); p = mp.mpf(1); n = 1
    while True:
        q = mp.exp(-2 * mp.pi * n * y)
        if q < TINY: break
        p *= (1 - q) ** 24; n += 1
    return mp.exp(6 * u - 2 * mp.pi * y) * p
gl = mp.calculus.quadrature.GaussLegendre(mp.mp)
BASE = gl.calc_nodes(6, mp.mp.prec)   # 96 nodes on [-1, 1]
U, W = [], []
for k in range(8):
    lo = mp.mpf(k) / 2
    for x, w in BASE:
        u = lo + (x + 1) / 4
        U.append(u); W.append(w / 4 * K(u))
NJ = len(U)
put('nodes', NJ)
mp.mp.dps = 60
def H(t, th, d=0):
    # d-th derivative in theta of H^Delta_t at theta, direct evaluation
    s = mp.mpf(0)
    for u, w in zip(U, W):
        e = w * mp.exp(t * u * u) * u ** d
        s += e * [mp.cos, lambda z: -mp.sin(z), lambda z: -mp.cos(z)][d](th * u)
    return s
put('H_0(0) = xi(1/2, Delta)/2 (maxlaw.py: 0.00154487936039502720604/2)', mp.nstr(H(0, 0), 25))
# fixed point grid
SH = 220; SCALE = mp.mpf(2) ** SH
TH0 = mp.mpf('0.5'); HSTEP = mp.mpf('0.02'); NK = 2576   # theta_k = 0.5 + 0.02 k, k < NK, up to 52
c0 = [int(mp.nint(mp.cos(TH0 * u) * SCALE)) for u in U]
c1 = [int(mp.nint(mp.cos((TH0 + HSTEP) * u) * SCALE)) for u in U]
ch = [int(mp.nint(mp.cos(HSTEP * u) * SCALE)) for u in U]
C = [c0, c1]
for k in range(2, NK):
    a, b = C[-1], C[-2]
    C.append([((2 * c * x) >> SH) - y for c, x, y in zip(ch, a, b)])
chk = max(abs(C[NK - 1][j] - int(mp.nint(mp.cos((TH0 + HSTEP * (NK - 1)) * U[j]) * SCALE))) for j in range(0, NJ, 37))
put('largest error of the recurrence at theta = 52, in units of 2^-220', chk)
THG = [float(TH0 + HSTEP * k) for k in range(NK)]
def grid(t):
    wv = [int(mp.nint(w * mp.exp(t * u * u) * SCALE)) for u, w in zip(U, W)]
    return [sum(map(mul, row, wv)) for row in C]
def zeros_at(t):
    v = grid(t)
    zs = []
    for k in range(NK - 1):
        if (v[k] > 0) != (v[k + 1] > 0):
            # cubic inverse interpolation through four grid points around the sign change
            ks = [max(0, min(NK - 4, k - 1)) + i for i in range(4)]
            xs = [THG[i] for i in ks]; ys = [v[i] / 2.0 ** (2 * SH) for i in ks]
            sc = max(abs(y) for y in ys); ys = [y / sc for y in ys]
            z = 0.0
            for i in range(4):
                term = xs[i]
                for j in range(4):
                    if j != i: term *= (0 - ys[j]) / (ys[i] - ys[j])
                z += term
            if not (THG[k] <= z <= THG[k + 1]):   # fall back to linear interpolation
                z = THG[k] + (THG[k + 1] - THG[k]) * v[k] / (v[k] - v[k + 1])
            zs.append(z)
    return zs
Z0 = zeros_at(mp.mpf(0))
put('zeros at t = 0', [round(z, 7) for z in Z0])
def crit(t, a, b):
    # the critical point of H^Delta_t between a and b, by Newton on dH/dtheta, from the midpoint
    x = mp.mpf((a + b) / 2)
    for _ in range(40):
        f1 = H(t, x, 1); f2 = H(t, x, 2)
        dx = f1 / f2
        x -= dx
        if abs(dx) < mp.mpf(10) ** -30: break
    return x, H(t, x)
def collision(ta, tb, a, b):
    # ta: the pair (a, b) is real; tb: it is gone. Bisection on the sign of the critical value between them
    xa, ga = crit(ta, a, b)
    sgn = ga > 0
    lo, hi = mp.mpf(ta), mp.mpf(tb); x = xa
    for _ in range(20):   # the pair may only have fallen inside one grid cell at tb: move tb back until it is gone
        xb, gb = crit(hi, float(x) - 0.3, float(x) + 0.3)
        if (gb > 0) != sgn: break
        lo, x = hi, xb; hi -= mp.mpf(1) / 100
    for _ in range(45):
        mid = (lo + hi) / 2
        x2, g = crit(mid, float(x) - 0.3, float(x) + 0.3)
        if (g > 0) == sgn:
            lo, x = mid, x2
        else:
            hi = mid
    return (lo + hi) / 2, x
labels = [round(z, 4) for z in Z0]
traj = {'0': [round(z, 6) for z in Z0]}
events = []
def advance(ta, za, tb, depth=0):
    # the zeros at tb < ta from those at ta, with the events between; a step with more than one event is cut into 8
    global labels
    zb = zeros_at(tb)
    drop = len(za) - len(zb)
    if drop == 0:
        return zb
    if drop == 1 and za[-1] > 51.0:
        events.append({'event': 'the highest zero passes 52', 'between t': [float(ta), float(tb)], 'zero at t = 0': labels[-1]})
        labels = labels[:-1]
        return zb
    if drop == 2:
        gaps = [za[j + 1] - za[j] for j in range(len(za) - 1)]
        j = gaps.index(min(gaps))
        tc, xc = collision(ta, tb, za[j], za[j + 1])
        D = labels[j + 1] - labels[j]
        events.append({'event': 'collision', 'pair at t = 0': [labels[j], labels[j + 1]], '-D^2/8': round(-D * D / 8, 4),
                       't_*': mp.nstr(tc, 10), 'theta_*': mp.nstr(xc, 8), 'ratio t_*/(-D^2/8)': round(float(tc) / (-D * D / 8), 3)})
        print(events[-1], round(time.time() - t0, 1), flush=True)
        labels = labels[:j] + labels[j + 2:]
        return zb
    if depth < 3:
        t, z = ta, za
        for k in range(1, 9):
            t2 = ta + (tb - ta) * k / 8
            z = advance(t, z, t2, depth + 1)
            traj[mp.nstr(t2, 6)] = [round(x, 6) for x in z]
            t = t2
        return z
    raise RuntimeError('unexpected change of the count from %d to %d between t = %s and %s' % (len(za), len(zb), ta, tb))
prev_t = mp.mpf(0); prev = Z0
for i in range(1, 111):   # backward, t = -0.01, ..., -1.1
    t = -mp.mpf(i) / 100
    cur = advance(prev_t, prev, t)
    traj[mp.nstr(t, 4)] = [round(z, 6) for z in cur]
    prev_t, prev = t, cur
fwd = []
for i in range(1, 51):   # forward, t = 0.01, ..., 0.5
    t = mp.mpf(i) / 100
    cur = zeros_at(t)
    traj[mp.nstr(t, 4)] = [round(z, 6) for z in cur]
    fwd.append(len(cur))
put('events', events)
put('zeros below 52 forward, t = 0.01, ..., 0.5 (counts)', sorted(set(fwd)))
out['trajectories'] = traj
out['seconds'] = time.time() - t0
json.dump(out, open('maxexp_flow.json', 'w'), indent=1)
print('done', round(time.time() - t0, 1))
