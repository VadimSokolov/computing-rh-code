# biane_clock.py: the index side clock of Polya's function as a hitting time, for Section sec:biane:clock of ch/biane.tex
# (Theorem thm:biane:clock, Table tab:biane:clock and Figure fig:biane:clock). Added for the book in October 2026, from the
# first author's note on Biane's paper on Polya (misc/reviews/nick-biane/, Section 7 there), recomputed independently.
# X = log R, with R the radial part of planar Brownian motion killed at rate 1/2 and Doob transformed by K_0, run on the
# clock tau = (1/2) int du / R_u^2, solves dX = -2 e^X (K_1/K_0)(e^X) dtau + sqrt(2) dW. Started at x_0 = log(2 pi) + 2.5
# (R_0 = 2 pi e^{2.5}), its hitting time of b = log(2 pi) has Laplace transform
#   L(theta) = K_{sqrt theta}(R_0) K_0(2 pi) / (K_0(R_0) K_{sqrt theta}(2 pi)),
# and from R_0 = +oo the hitting time is tau_oo = sum_n eps_n / nu_n^2, eps_n standard exponential, nu_n the zeros of
# nu -> K_{i nu}(2 pi). The script simulates 100000 paths by the Euler scheme with step 4 x 10^{-6} (the step of the note),
# recording two hitting times on the same paths: the plain one (first step that ends below b, as in the note) and the one with
# the Brownian bridge correction (a step from x to x' both above b crosses b with probability exp(-(x-b)(x'-b)/dt), and the
# crossing is then dated to the middle of the step). The plain scheme overshoots the barrier and is late by about
# 0.5826 sqrt(2 dt) / |drift at b|; the bridge correction removes this to O(dt). Exact mean and standard deviation from the
# expansion K_{sqrt theta}(x) = K_0(x) + theta k_2(x)/2 + theta^2 k_4(x)/24 + ..., k_{2j}(x) = int_0^oo w^{2j} e^{-x cosh w} dw;
# exact L(theta) on 0 <= theta <= 400 from Arb, against the Monte Carlo averages of e^{-theta tau} with their standard
# errors; the exact densities of the finite start clock and of tau_oo by Talbot's inversion of the Laplace transform (fixed
# Talbot contour, M = 32 nodes, Arb at 704 bits, since Arb's K_mu(x) loses about 2x/log 2 bits to cancellation, 221 bits
# at x = R_0), checked by their total mass and mean on the grid.
# Run from any folder: python3 biane_clock.py; writes biane_clock.json there (about ten minutes on one core).
import json
import math
import time

import mpmath as mp
import numpy as np
from flint import acb, arb, ctx
from scipy.special import k0e, k1e

t0 = time.time()
out = {}


def rec(k, v):
    out[k] = v
    print(k, v, round(time.time() - t0, 1), flush=True)


def af(x):
    # the midpoint of an Arb ball as a Python float
    return float(x.mid().str(25, radius=False))


mp.mp.dps = 40
PI = mp.pi
b = 2 * PI
R0 = 2 * PI * mp.e ** mp.mpf('2.5')


# ---------------------------------------------------------------- exact moments
def kmom(j, x):
    # e^x k_{2j}(x) on a finite range (the integrand is below e^{-400} beyond W)
    W = mp.acosh(1 + 400 / x)
    return mp.quad(lambda w: w ** (2 * j) * mp.exp(-x * (mp.cosh(w) - 1)), [0, W / 16, W / 8, W / 4, W / 2, W])


def a12(x):
    k0, k2, k4 = kmom(0, x), kmom(1, x), kmom(2, x)
    return k2 / (2 * k0), k4 / (24 * k0) - k2 ** 2 / (8 * k0 ** 2)


a1b, a2b = a12(b)
a1r, a2r = a12(R0)
mean_ex, var_ex = a1b - a1r, 2 * (a2r - a2b)
mean_oo, var_oo = a1b, -2 * a2b
rec('exact: R_0, mean and sd from R_0, mean and sd of tau_oo',
    [float(R0), float(mean_ex), float(mp.sqrt(var_ex)), float(mean_oo), float(mp.sqrt(var_oo))])

# ---------------------------------------------------------------- exact Laplace transforms (Arb)
ctx.prec = 704
Rb, Bb = arb(mp.nstr(R0, 60)), arb(2) * arb.pi()


def LT(th, start=True):
    """L(theta) for complex theta (acb); the ratio is even in the order, so the branch of sqrt is irrelevant."""
    mu = acb(th).sqrt()
    num = acb(Bb).bessel_k(acb(0))
    den = acb(Bb).bessel_k(mu)
    if start:
        num = num * acb(Rb).bessel_k(mu)
        den = den * acb(Rb).bessel_k(acb(0))
    return num / den


THG = [float(x) for x in np.linspace(0, 400, 81)]
LTex = []
for th in THG:
    v = LT(acb(th))
    assert v.real.rad() < 1e-30 * abs(af(v.real)) + 1e-300
    LTex.append(af(v.real))
rec('exact L(theta) at theta = 0, 100, 200, 300, 400', [LTex[0], LTex[20], LTex[40], LTex[60], LTex[80]])


def talbot(t, start, M=32):
    """Fixed Talbot inversion (Abate and Valko 2004) of LT at time t > 0."""
    tt = arb(mp.nstr(t, 30))
    r = arb(2 * M) / (5 * tt)
    tot = (LT(acb(r), start) * (r * tt).exp()).real / 2
    for k in range(1, M):
        th = arb(k) * arb.pi() / M
        c = th.cos() / th.sin()
        s = acb(r * th * c, r * th)
        sig = th + (th * c - 1) * c
        tot += (acb(s * tt).exp() * LT(s, start) * acb(1, sig)).real
    return r / M * tot


TD = [float(x) for x in np.linspace(0.015, 0.20, 75)]
def dens_list(start):
    vals, worst = [], 0.0
    for t in TD:
        v = talbot(t, start)
        vals.append(af(v))
        worst = max(worst, float(v.rad().str(5, radius=False)))
    return vals, worst


dens_start, rad_start = dens_list(True)
dens_oo, rad_oo = dens_list(False)
rec('Talbot densities: largest radius of the Arb balls (finite start, tau_oo)', [rad_start, rad_oo])
h = TD[1] - TD[0]
mass = [h * (sum(d) - (d[0] + d[-1]) / 2) for d in (dens_start, dens_oo)]
mean_grid = [h * (sum(t * x for t, x in zip(TD, d)) - (TD[0] * d[0] + TD[-1] * d[-1]) / 2) for d in (dens_start, dens_oo)]
rec('Talbot densities on [0.015, 0.2]: total mass and mean (finite start, tau_oo)', [mass, mean_grid])
# check of the inversion at one point by a doubled number of nodes
rec('Talbot check at t = 0.07: M = 32 and M = 64 (finite start)',
    [af(talbot(0.07, True)), af(talbot(0.07, True, 64))])

# ---------------------------------------------------------------- Monte Carlo
rng = np.random.default_rng(20261009)
N, dt = 100000, 4e-6
lb = math.log(2 * math.pi)
x = np.full(N, lb + 2.5)
tau_plain = np.full(N, np.nan)
tau_bridge = np.full(N, np.nan)
active = np.arange(N)
sq = math.sqrt(2 * dt)
t = 0.0
steps = 0
while active.size:
    xa = x[active]
    r = np.exp(xa)
    xn = xa - 2 * r * (k1e(r) / k0e(r)) * dt + sq * rng.standard_normal(xa.size)
    u = rng.random(xa.size)
    below = xn <= lb
    with np.errstate(over='ignore', under='ignore'):
        pc = np.where(below, 1.0, np.exp(-np.maximum(xa - lb, 0) * np.maximum(xn - lb, 0) / dt))
    newb = np.isnan(tau_bridge[active]) & (below | (u < pc))
    tau_bridge[active[newb]] = t + dt / 2
    tau_plain[active[below]] = t + dt
    x[active] = xn
    active = active[~below]
    t += dt
    steps += 1
rec('Monte Carlo: paths, step, number of steps to the last hit', [N, dt, steps])
drift_b = 2 * 2 * math.pi * k1e(2 * math.pi) / k0e(2 * math.pi)
rec('drift at the barrier and the predicted lateness 0.5826 sqrt(2 dt)/|drift| of the plain scheme',
    [drift_b, 0.5826 * math.sqrt(2 * dt) / drift_b])
mc = {}
for name, tau in [('plain', tau_plain), ('bridge', tau_bridge)]:
    m, s = float(tau.mean()), float(tau.std(ddof=1))
    se = s / math.sqrt(N)
    mc[name] = {'mean': m, 'sd': s, 'se': se, 'z of the mean': (m - float(mean_ex)) / se,
                'quantiles 0.01 0.5 0.99': [float(q) for q in np.quantile(tau, [0.01, 0.5, 0.99])]}
    print(name, mc[name], flush=True)
out['mc'] = mc
# Laplace transform: MC averages of e^{-theta tau} with standard errors, against the exact values
rows = []
for th, ex in zip(THG, LTex):
    e = np.exp(-th * tau_bridge)
    m, se = float(e.mean()), float(e.std(ddof=1) / math.sqrt(N))
    rows.append([th, ex, m, se, (m - ex) / se if se > 0 else 0.0])
out['LT'] = rows
zs = [rw[4] for rw in rows[1:]]
rec('Laplace transform (bridge): largest |z| over theta in (0, 400], z at 100, 200, 400; MC/exact at 400',
    [max(abs(z) for z in zs), rows[20][4], rows[40][4], rows[80][4], rows[80][2] / rows[80][1]])
rows_p = []
for th, ex in zip(THG, LTex):
    e = np.exp(-th * tau_plain)
    m, se = float(e.mean()), float(e.std(ddof=1) / math.sqrt(N))
    rows_p.append([th, ex, m, se, (m - ex) / se if se > 0 else 0.0])
out['LT_plain'] = rows_p
rec('Laplace transform (plain): z at 100, 200, 400; MC/exact at 400', [rows_p[20][4], rows_p[40][4], rows_p[80][4],
                                                                        rows_p[80][2] / rows_p[80][1]])
edges = np.linspace(0.015, 0.20, 75)
hb, _ = np.histogram(tau_bridge, bins=edges, density=True)
out['fig'] = {'theta': THG, 'LT_exact': LTex, 'LT_mc': [rw[2] for rw in rows], 'LT_se': [rw[3] for rw in rows],
              't': TD, 'dens_start': dens_start, 'dens_oo': dens_oo, 'hist_edges': [float(e) for e in edges],
              'hist': [float(v) for v in hb]}
out['seconds'] = round(time.time() - t0, 1)
json.dump(out, open('biane_clock.json', 'w'), indent=1)
print('done', out['seconds'])
