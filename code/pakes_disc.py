"""pakes_disc.py: the disc of a pair of zeros and the parabolic sliver of Section ch:boschB (ch/boschB.tex), that is
Lemma bs:lem:disc, Theorem bs:thm:local, the paragraph and Remark bs:rem:rightend after it and Figure bs:fig:disc, and the
bound on the off line part R_off in Section ag:sec:exactly (ch/almost.tex).

Added for the book in October 2026, for the merge of the first author's chapter on Pakes, the Lambert law and the Thorin
condition (misc/reviews/nick-pakes/, triage-B.md, items B16 to B19 and B24 to B26). Adapted from the independent checks
B_disc.py and B_disc_tip.py of that triage (misc/reviews/nick-pakes/check/); no script of the chapter is used.

Notation of Section ch:boschB. psi(s) = G'/G(s) = sum_gamma (s + gamma^2)^(-1) over the zeros gamma of Xi(z) = xi(1/2 + iz)
with Re gamma > 0, counted with multiplicity. A nonreal zero is a + ib or a - ib with a > T0 and 0 < b < 1/2 (the zeros
1/2 -+ b + ia of xi); the two terms of the pair are 1/(s - c) + 1/(s - conj c) with c = b^2 - a^2 + 2iab, and D_gamma is the
open disc {|s - Re c| < Im c}, whose diameter is the segment from conj c to c. The sliver is
E_T0 = {x + iy : x < -T0^2 + T0 + 1/4, y^2 < 1/2 - x}.

The output pakes_disc.json has these parts.
  sympy      the identities in the proofs of Lemma bs:lem:disc and Theorem bs:thm:local, checked symbolically: the imaginary
             parts, the numerator 2y(d + y^2 - c2^2), the poles of a pair, the rightmost point of D_gamma, its derivatives and
             its supremum, the identity behind y^2 < 1/2 - x, the position of the poles, the sign on the real axis, the
             envelope that makes the bound on y sharp, and the height of the sliver at its tip.
  lemma      200000 random tests of the equivalence of Lemma bs:lem:disc in double precision, and 20000 at 50 digits at
             relative distances 1e-30 to 1e-3 from the circle.
  discs      for T0 = 2, 100, 1e6, 3e12, 40000 boundary points of discs D_gamma at 60 digits, half of them with b -> 1/2 and
             a -> T0, angles concentrated at the rightmost point: the largest values of x - (-T0^2 + T0 + 1/4) and of
             y^2 - (1/2 - x) (both negative: every disc lies in E_T0).
  sharp      the gap between the right end and the rightmost point of D_gamma at b = 1/2 - 1e-20, a = T0 + 1e-20, and the
             relative gap (1/2 - x - sup y^2)/(1/2 - x) at x = -(T0 + 10)^2 for the disc with b = 1/2 - 1e-12, a^2 = 3b^2 - x.
  angle      the angle of E_T0 seen from the origin at its tip and at x = -1e30, for T0 = 3e12, against 1/T0, the opening of the
             sector of Corollary bs:cor:uncond.
  right_end  the points of Remark bs:rem:rightend: for T0 = 100 and 3e12, b = 1/2 - 1e-20 and a = T0 + 1e-20, the point of
             D_gamma at height 0.001 and 0.1 to the left of its rightmost point, which lies to the right of -(T0 - 1/2)^2, and
             the imaginary part of the pair there.
  toy        a toy reciprocal Xi function at 30 digits: the first 200 ordinates of zeros_1700.txt on the line and the pair
             60 -+ 0.3i (the zeros 1/2 -+ 0.3 + 60i of xi), so that D_gamma has centre -3599.91 and radius 36: Im psi on a
             grid of 7260 points around D_gamma, on the half circle of radius 1.0001 x 36, just outside the sliver of the toy
             (T0 = 59, below which the toy has only real zeros), and psi on (0, infinity).
  r_off      the constants of the bound |R_off(u)| < 2 Sigma_0 < 3e-12 on |u| <= T0^2/2 in Section ag:sec:exactly, with the
             main term of Sigma_0 of Theorem ag:thm:explicit.
The figure pakes_disc.pdf (Figure bs:fig:disc) has two panels: left, the sliver E_2 with the discs D_gamma of six
hypothetical zeros a + ib with a > 2; right, the set where Im psi > 0 for the toy, on a grid of 481 x 300 points in double
precision, with the circle bounding D_gamma.
Run on Hopper, from a folder that holds zeros_1700.txt (about a minute on one core): python3 pakes_disc.py [zeros_1700.txt];
writes pakes_disc.json and pakes_disc.pdf there (the book's copies are data/pakes_disc.json and fig/pakes_disc.pdf).
"""
import json
import random
import sys
import time

import matplotlib
import mpmath as mp
import numpy as np
import sympy as sp

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

T_START = time.time()
ZFILE = sys.argv[1] if len(sys.argv) > 1 else 'zeros_1700.txt'
out = {}


def S(v, n=6):
    return mp.nstr(v, n)


# 1. symbolic identities
x, y, c1, c2, d, e = sp.symbols('x y c1 c2 d e', real=True)
a, b, T0 = sp.symbols('a b T0', positive=True)
HALF, QUARTER = sp.Rational(1, 2), sp.Rational(1, 4)
s = x + sp.I * y
c = c1 + sp.I * c2
sym = {}
inv1 = sp.expand_complex(1 / (s - c))
inv2 = sp.expand_complex(1 / (s - sp.conjugate(c)))
sym['Im 1/(s-c) + (y-c2)/|s-c|^2'] = str(sp.simplify(sp.im(inv1) + (y - c2) / ((x - c1) ** 2 + (y - c2) ** 2)))
sym['Im 1/(s-conj c) + (y+c2)/|s-conj c|^2'] = str(sp.simplify(sp.im(inv2) + (y + c2) / ((x - c1) ** 2 + (y + c2) ** 2)))
num = (y - c2) * (d + (y + c2) ** 2) + (y + c2) * (d + (y - c2) ** 2)
sym['(y-c2)(d+(y+c2)^2) + (y+c2)(d+(y-c2)^2) - 2y(d+y^2-c2^2)'] = str(sp.expand(num - 2 * y * (d + y ** 2 - c2 ** 2)))
target = -2 * y * ((x - c1) ** 2 + y ** 2 - c2 ** 2) / (((x - c1) ** 2 + (y - c2) ** 2) * ((x - c1) ** 2 + (y + c2) ** 2))
sym['Im pair - (-2y((x-c1)^2+y^2-c2^2)/(|s-c|^2 |s-conj c|^2))'] = str(sp.simplify(sp.together(sp.im(inv1 + inv2) - target)))
cc = b ** 2 - a ** 2 + 2 * sp.I * a * b  # -conj(gamma)^2 for gamma = a + ib
sym['-(a - ib)^2 - c'] = str(sp.expand(-(a - sp.I * b) ** 2 - cc))
sym['-(a + ib)^2 - conj c'] = str(sp.expand(-(a + sp.I * b) ** 2 - sp.conjugate(cc)))
R = b ** 2 - a ** 2 + 2 * a * b  # rightmost point of D_gamma
sym['R - (2b^2 - (a-b)^2)'] = str(sp.expand(R - (2 * b ** 2 - (a - b) ** 2)))
sym['R + (a-b)^2'] = str(sp.expand(R + (a - b) ** 2))
sym['dR/db'] = str(sp.diff(R, b))
sym['dR/da'] = str(sp.diff(R, a))
sym['R(a=T0, b=1/2)'] = str(sp.expand(R.subs({a: T0, b: HALF})))
sym['R(a=T0, b=1/2) + (T0-1/2)^2'] = str(sp.expand(R.subs({a: T0, b: HALF}) + (T0 - HALF) ** 2))
lhs = 4 * b ** 2 * (b ** 2 - x + e) - e ** 2  # (Im c)^2 - (x - Re c)^2 with a^2 = b^2 - x + e
sym['4b^2(b^2-x+e) - e^2 - (8b^4 - 4b^2 x - (e-2b^2)^2)'] = str(sp.expand(lhs - (8 * b ** 4 - 4 * b ** 2 * x - (e - 2 * b ** 2) ** 2)))
sym['8b^4 - 4b^2 x - (1/2 - x - (1-4b^2)(1/2+2b^2-x))'] = str(
    sp.expand(8 * b ** 4 - 4 * b ** 2 * x - (HALF - x - (1 - 4 * b ** 2) * (HALF + 2 * b ** 2 - x))))
sym['1/2 - Re c - (Im c)^2 - (a^2(1-4b^2) + 1/2 - b^2)'] = str(
    sp.expand(HALF - (b ** 2 - a ** 2) - (2 * a * b) ** 2 - (a ** 2 * (1 - 4 * b ** 2) + HALF - b ** 2)))
xr = sp.symbols('xr', real=True)
pair_real = 1 / (xr - cc) + 1 / (xr - sp.conjugate(cc))
sym['pair at real x - 2(x-Re c)/((x-Re c)^2+(Im c)^2)'] = str(
    sp.simplify(sp.together(pair_real - 2 * (xr - (b ** 2 - a ** 2)) / ((xr - (b ** 2 - a ** 2)) ** 2 + (2 * a * b) ** 2))))
sym['envelope: e = 2b^2 gives a^2 = b^2 - x + e ='] = str(sp.expand(b ** 2 - x + 2 * b ** 2))
sym['limit b->1/2 of 8b^4 - 4b^2 x'] = str(sp.limit(8 * b ** 4 - 4 * b ** 2 * x, b, HALF))
sym['1/2 - (-T0^2 + T0 + 1/4), factored'] = str(sp.factor(HALF - (-T0 ** 2 + T0 + QUARTER)))
out['sympy'] = sym
print('sympy', json.dumps(sym, indent=1), flush=True)

# 2. random tests of Lemma bs:lem:disc
rng = random.Random(1)
bad = skipped = 0
N1 = 200000
for _ in range(N1):
    c1v, c2v = rng.uniform(-10, 10), rng.uniform(0.01, 5)
    sv = complex(rng.uniform(-15, 15), rng.uniform(1e-3, 8))
    cv = complex(c1v, c2v)
    if abs(abs(sv - c1v) - c2v) < 1e-9:
        skipped += 1
        continue
    imv = (1 / (sv - cv) + 1 / (sv - cv.conjugate())).imag
    if (imv > 0) != (abs(sv - c1v) < c2v):
        bad += 1
lem = {'double': {'tests': N1, 'skipped within 1e-9 of the circle': skipped, 'disagreements': bad}}
mp.mp.dps = 50
bad = 0
N2 = 20000
for k in range(N2):
    c1v, c2v = mp.mpf(rng.uniform(-10, 10)), mp.mpf(rng.uniform(0.01, 5))
    phi = mp.pi * mp.mpf(rng.uniform(1e-6, 1 - 1e-6))
    rel = mp.mpf(10) ** (-rng.uniform(3, 30))
    inside = bool(k % 2)
    rr = c2v * (1 - rel) if inside else c2v * (1 + rel)
    sv = mp.mpc(c1v + rr * mp.cos(phi), rr * mp.sin(phi))
    cv = mp.mpc(c1v, c2v)
    imv = mp.im(1 / (sv - cv) + 1 / (sv - mp.conj(cv)))
    if (imv > 0) != inside:
        bad += 1
lem['50 digits, relative distance 1e-30 to 1e-3 from the circle'] = {'tests': N2, 'disagreements': bad}
out['lemma'] = lem
print('lemma', lem, flush=True)

# 3. every disc lies in E_T0
mp.mp.dps = 60
rng = random.Random(11)
discs = []
for T0s in ['2', '100', '1e6', '3e12']:
    T = mp.mpf(T0s)
    right_end = -T ** 2 + T + mp.mpf('0.25')
    worst_x = worst_y = -mp.inf
    for k in range(40000):
        if k % 2:
            bv = mp.mpf('0.5') - mp.mpf(10) ** (-rng.uniform(1, 30))
            av = T + mp.mpf(10) ** (-rng.uniform(1, 30))
        else:
            bv = mp.mpf(rng.uniform(1e-6, 0.5 - 1e-12))
            av = T * (1 + mp.mpf(rng.random()) ** 4 * 3) + mp.mpf(rng.random())
        cen, rad = bv ** 2 - av ** 2, 2 * av * bv
        ang = mp.pi * mp.mpf(rng.random()) ** 3
        rr = rad * (1 - mp.mpf('1e-40'))
        px, py = cen + rr * mp.cos(ang), rr * mp.sin(ang)
        worst_x = max(worst_x, px - right_end)
        worst_y = max(worst_y, py ** 2 - (mp.mpf('0.5') - px))
    discs.append({'T0': S(T, 3), 'points': 40000, 'max of x - (-T0^2 + T0 + 1/4)': S(worst_x),
                  'max of y^2 - (1/2 - x)': S(worst_y)})
out['discs'] = discs
print('discs', discs, flush=True)

# 4. sharpness of both bounds
sharp = []
for T0s in ['2', '100', '1e6', '3e12']:
    T = mp.mpf(T0s)
    bv, av = mp.mpf('0.5') - mp.mpf('1e-20'), T + mp.mpf('1e-20')
    gap_right = (-T ** 2 + T + mp.mpf('0.25')) - (bv ** 2 - av ** 2 + 2 * av * bv)
    xx = -(T + 10) ** 2
    bv2 = mp.mpf('0.5') - mp.mpf('1e-12')
    av2 = mp.sqrt(3 * bv2 ** 2 - xx)
    cen2, rad2 = bv2 ** 2 - av2 ** 2, 2 * av2 * bv2
    ytop2 = rad2 ** 2 - (xx - cen2) ** 2  # largest y^2 over this disc at abscissa xx
    sharp.append({'T0': S(T, 3), 'right end minus rightmost point, b = 1/2 - 1e-20, a = T0 + 1e-20': S(gap_right),
                  'x': S(xx, 15), 'a - T0 for a^2 = 3b^2 - x': S(av2 - T),
                  'sup y^2 minus (8b^4 - 4b^2 x)': S(ytop2 - (8 * bv2 ** 4 - 4 * bv2 ** 2 * xx)),
                  'relative gap (1/2 - x - sup y^2)/(1/2 - x)': S(((mp.mpf('0.5') - xx) - ytop2) / (mp.mpf('0.5') - xx))})
out['sharp'] = sharp
print('sharp', sharp, flush=True)

# 5. the angle of E_T0 seen from the origin
T = mp.mpf('3e12')
x0 = -T ** 2 + T + mp.mpf('0.25')
y0 = mp.sqrt(mp.mpf('0.5') - x0)
tip = mp.atan(y0 / (-x0))
along = [{'x': S(xv, 3), 'angle': S(mp.atan(mp.sqrt(mp.mpf('0.5') - xv) / (-xv)), 6),
          'angle * |x|^(1/2)': S(mp.atan(mp.sqrt(mp.mpf('0.5') - xv) / (-xv)) * mp.sqrt(-xv), 10)}
         for xv in [x0, 10 * x0, mp.mpf('-1e26'), mp.mpf('-1e28'), mp.mpf('-1e30')]]
out['angle'] = {'T0': '3e12', 'tip x0 = -T0^2 + T0 + 1/4': S(x0, 20), 'height at the tip': S(y0, 20),
                'height minus (T0 - 1/2)': S(y0 - (T - mp.mpf('0.5'))), 'angle at the tip': S(tip, 15),
                '1/T0, the opening of the sector of Corollary bs:cor:uncond': S(1 / T, 15), 'T0 * angle at the tip - 1': S(T * tip - 1),
                'along the sliver': along}
print('angle', out['angle'], flush=True)

# 6. Remark bs:rem:rightend
ex = []
for T0s in ['100', '3e12']:
    T = mp.mpf(T0s)
    bv, av = mp.mpf('0.5') - mp.mpf('1e-20'), T + mp.mpf('1e-20')
    cen, rad = bv ** 2 - av ** 2, 2 * av * bv
    sv = mp.mpc(cen + rad - mp.mpf('0.1'), mp.mpf('1e-3'))
    cv = mp.mpc(cen, rad)
    pair_im = mp.im(1 / (sv - cv) + 1 / (sv - mp.conj(cv)))
    ex.append({'T0': S(T, 3), 's': S(sv.real, 30) + ' + 0.001i', 'rightmost point of D_gamma': S(cen + rad, 30),
               '-(T0 - 1/2)^2': S(-(T - mp.mpf('0.5')) ** 2, 30), 'x + (T0 - 1/2)^2': S(sv.real + (T - mp.mpf('0.5')) ** 2, 10),
               's in D_gamma': bool(abs(sv - cen) < rad),
               's in {x < -(T0-1/2)^2, y^2 < 1/2 - x}': bool(sv.real < -(T - mp.mpf('0.5')) ** 2 and sv.imag ** 2 < mp.mpf('0.5') - sv.real),
               's in E_T0': bool(sv.real < -T ** 2 + T + mp.mpf('0.25') and sv.imag ** 2 < mp.mpf('0.5') - sv.real),
               'Im of the pair at s': S(pair_im, 6)})
out['right_end'] = ex
print('right_end', ex, flush=True)

# 7. the toy reciprocal Xi function
mp.mp.dps = 30
gam = [mp.mpf(line.split()[0]) for line in open(ZFILE) if line.strip() and not line.startswith('#')][:200]
A0, B0 = mp.mpf(60), mp.mpf('0.3')
cpole = mp.mpc(B0 ** 2 - A0 ** 2, 2 * A0 * B0)
cen, rad = cpole.real, cpole.imag


def psi_toy(z):
    return sum(1 / (z + g * g) for g in gam) + 1 / (z - cpole) + 1 / (z - mp.conj(cpole))


pos, grid, in_half_disc = [], 0, 0
for i in range(-60, 61):
    for j in range(1, 61):
        z = mp.mpc(cen + rad * (i + mp.mpf('0.23')) / 40, rad * (j + mp.mpf('0.37')) / 40)
        grid += 1
        inside = abs(z - cen) < rad
        in_half_disc += inside
        if mp.im(psi_toy(z)) > 0:
            pos.append((z, inside))
circ = [mp.im(psi_toy(cen + rad * mp.mpf('1.0001') * mp.expj(mp.pi * k / 200))) for k in range(1, 200)]
T_toy = mp.mpf(59)  # the toy has only real zeros with 0 < Re gamma <= 59
x0t = -T_toy ** 2 + T_toy + mp.mpf('0.25')
outside = [mp.im(psi_toy(mp.mpc(x0t + mp.mpf('0.01'), (T_toy - mp.mpf('0.5')) * k / 50))) for k in range(1, 51)]
outside += [mp.im(psi_toy(mp.mpc(xv, mp.sqrt(mp.mpf('0.5') - xv) + mp.mpf('0.01'))))
            for xv in [x0t - mp.mpf(10) ** (mp.mpf(k) / 20) for k in range(0, 101)]]
out['toy'] = {'zeros on the line': len(gam), 'largest ordinate': S(gam[-1], 10), 'pair': '60 +- 0.3i',
              'pole c in the upper half plane': S(cpole.real, 10) + ' + ' + S(cpole.imag, 10) + 'i',
              'disc centre': S(cen, 10), 'disc radius': S(rad, 10),
              'grid': 'x = centre + r(i + 0.23)/40, i = -60..60; y = r(j + 0.37)/40, j = 1..60',
              'grid points': grid, 'grid points in D_gamma': int(in_half_disc), 'points with Im psi > 0': len(pos),
              'all of them in D_gamma': all(p[1] for p in pos),
              'largest |z - centre|/r among them': S(max(abs(p[0] - cen) for p in pos) / rad, 6),
              'max Im psi on the half circle of radius 1.0001 r (199 points)': S(max(circ), 6),
              'toy T0, right end of its sliver': [S(T_toy, 3), S(x0t, 10)],
              'max Im psi just outside the toy sliver (151 points)': S(max(outside), 6),
              'psi on (0, infinity)': {S(mp.mpf(v), 3): S(mp.re(psi_toy(mp.mpf(v))), 8) for v in ['1e-6', '1', '100', '1e4', '1e8']}}
print('toy', out['toy'], flush=True)

# 8. the bound on R_off in Section ag:sec:exactly
mp.mp.dps = 30
T = mp.mpf('3e12')
out['r_off'] = {'T0': '3e12', 'main term (log(T0/2pi) + 1)/(2pi T0) of Sigma_0': S((mp.log(T / (2 * mp.pi)) + 1) / (2 * mp.pi * T), 6),
                'Sigma_0 < 1.48e-12 (Theorem ag:thm:explicit), so 2 Sigma_0 <': S(2 * mp.mpf('1.48e-12'), 4),
                'T0^2/2': S(T ** 2 / 2, 6), '|right end of E_T0| = T0^2 - T0 - 1/4': S(T ** 2 - T - mp.mpf('0.25'), 20),
                'the disc |u| <= T0^2/2 misses E_T0': bool(T ** 2 - T - mp.mpf('0.25') > T ** 2 / 2)}
print('r_off', out['r_off'], flush=True)

# 9. Figure bs:fig:disc
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
fig, (axl, axr) = plt.subplots(1, 2, figsize=(6.6, 2.3), gridspec_kw={'width_ratios': [2.6, 2.0]})
T0f = 2.0
x0f = -T0f ** 2 + T0f + 0.25
xs = np.linspace(-33.5, x0f, 600)
yp = np.sqrt(0.5 - xs)
for sgn in (1, -1):
    axl.plot(xs, sgn * yp, 'k--', lw=0.9)
axl.plot([x0f, x0f], [-np.sqrt(0.5 - x0f), np.sqrt(0.5 - x0f)], 'k--', lw=0.9)
th = np.linspace(0, 2 * np.pi, 500)
fig_discs = [(2.02, 0.49), (2.6, 0.49), (3.3, 0.49), (4.2, 0.49), (5.2, 0.49), (3.0, 0.2)]
for av, bv in fig_discs:
    cf, rf = bv ** 2 - av ** 2, 2 * av * bv
    axl.plot(cf + rf * np.cos(th), rf * np.sin(th), color='0.4', lw=0.7)
    axl.plot([cf, cf], [rf, -rf], 'k.', ms=3)
axl.axhline(0, color='k', lw=0.5)
axl.plot([0], [0], 'ko', ms=2.5)
axl.set_xlim(-33.5, 2)
axl.set_ylim(-6.5, 6.5)
axl.set_aspect('equal')
axl.set_xlabel(r'$\mathrm{Re}\,s$')
axl.set_ylabel(r'$\mathrm{Im}\,s$')
g2 = np.array([float(g) ** 2 for g in gam])
cp = complex(cpole)
X = np.linspace(-3648, -3552, 481)
Y = np.linspace(0.05, 45, 300)
XX, YY = np.meshgrid(X, Y)
Z = XX + 1j * YY
P = np.zeros_like(Z)
for g in g2:
    P += 1 / (Z + g)
P += 1 / (Z - cp) + 1 / (Z - np.conj(cp))
IM = P.imag
axr.contourf(XX, YY, IM, levels=[0, float(IM.max()) + 1], colors=['0.75'])
half_th = np.linspace(0, np.pi, 300)
axr.plot(cp.real + cp.imag * np.cos(half_th), cp.imag * np.sin(half_th), 'k-', lw=0.9)
axr.plot([cp.real], [cp.imag], 'kx', ms=5)
axr.set_xlim(-3648, -3552)
axr.set_ylim(0, 45)
axr.set_aspect('equal')
axr.set_xticks([-3640, -3600, -3560])
axr.set_xlabel(r'$\mathrm{Re}\,s$')
fig.tight_layout()
fig.savefig('pakes_disc.pdf')
in_disc_fine = np.abs(Z - cp.real) < cp.imag
out['figure'] = {'left': {'T0': 2, 'right end': x0f, 'discs (a, b)': fig_discs},
                 'right': {'grid': '481 x 300, x in [-3648, -3552], y in [0.05, 45], double precision',
                           'points with Im psi > 0': int((IM > 0).sum()), 'points in D_gamma': int(in_disc_fine.sum()),
                           'points with Im psi > 0 outside D_gamma': int(((IM > 0) & ~in_disc_fine).sum())}}
out['seconds'] = round(time.time() - T_START, 1)
json.dump(out, open('pakes_disc.json', 'w'), indent=1)
print('figure', out['figure'], flush=True)
print('done in', out['seconds'], 's', flush=True)
