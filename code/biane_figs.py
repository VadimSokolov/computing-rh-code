# biane_figs.py: Figures fig:biane:kernels, fig:biane:wick, fig:biane:clock, fig:biane:flow of Section ch:biane
# (ch/biane.tex) and fig:biane:cm, fig:biane:traces of Subsection sec:biane:cm (ch/strip.tex). Added for the book in October 2026;
# they redraw, with corrections, the seven figures of the first author's note on Biane's paper on Polya
# (misc/reviews/nick-biane/), whose counting and kernel figures are merged into one.
#  biane_kernels.pdf: left, the kernels Phi of eq:ch2:Phi, Polya's symmetric asymptote 8 pi^2 cosh(9u/2) e^{-2 pi cosh 2u}
#    and e^{-2 pi cosh 2u}, each divided by its value at u = 0; right, the zero counting functions of Xi, Xi* and P on
#    (0, 100] (from biane_polya.json) against N_0(T) + 7/8 and N_0(T) - 1/4.
#  biane_wick.pdf: left, 1/E exp(-u J^(kappa)_1) for kappa = 1, 0, -1 on both sides of u = 0, with its zeros at minus
#    the Thorin atoms pi^2 n^2/2 - kappa/2; right, 1/E exp(-u H) for the first passage H of Brownian motion with drift 1 at
#    the level 1: real part and modulus, constant beyond the branch point u = -1/2.
#  biane_clock.pdf: the clock from R_0 = 2 pi e^{2.5}: Laplace transform, exact against Monte Carlo, and the histogram of the
#    simulated hitting times against the exact densities of the finite start clock and of tau_oo (from biane_clock.json).
#  biane_flow.pdf: the real zeros in (0, 90) of P_t and Xi_t for -60 <= t <= 0 (from biane_flow.json, where t is lambda).
#  biane_cm.pdf: the weight T^n e^{-Ty} of a condition (n, y), normalised, and the sign factor cos(n phi - delta y) of an off
#    line zero (from biane_traces.json, where y is v).
#  biane_traces.pdf: the Poisson trace K_{1/2}(y) and Polya's trace K^P(y) = sum_n e^{-2 nu_n y}, scaled by 2 pi y, and their
#    difference D(y) against 9/8 and 9/8 - c y log(1/y) (from biane_traces.json, where y is v).
# The script also writes biane_figs.json with the smallest normalised gaps between consecutive real zeros of P_lambda and
# Xi_lambda, with the local density log(t/2 pi)/(2 pi) taken at the midpoint of the gap and at its lower end, and the mean
# speeds (t_k(-0.5) - t_k(0))/0.5 of the real zeros in (0, 90) over -0.5 <= lambda <= 0, to compare with pi/8.
# Run from a folder holding biane_polya.json, biane_flow.json, biane_clock.json and biane_traces.json: python3 biane_figs.py;
# writes the six figures to fig/ (created if absent) and biane_figs.json to the folder itself.
import json
import math
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': 0.25, 'axes.spines.top': False, 'axes.spines.right': False})
C1, C2, C3, C4 = '#1b4f72', '#c0392b', '#27ae60', '#7d3c98'
os.makedirs('fig', exist_ok=True)
TWO_PI = 2 * math.pi
EG = 0.5772156649015329
out = {}

# ---------------------------------------------------------------- 1. kernels and counting functions
def Phi(u):
    u = np.asarray(u, dtype=float)
    tot = np.zeros_like(u)
    for n in range(1, 12):
        tot += (2 * math.pi ** 2 * n ** 4 * np.exp(4.5 * u) - 3 * math.pi * n ** 2 * np.exp(2.5 * u)) * np.exp(-math.pi * n * n * np.exp(2 * u))
    return 2 * tot


uu = np.linspace(0, 1.2, 400)
k1 = Phi(uu) / Phi(np.array([0.0]))[0]
k2 = np.cosh(4.5 * uu) * np.exp(-TWO_PI * (np.cosh(2 * uu) - 1))
k3 = np.exp(-TWO_PI * (np.cosh(2 * uu) - 1))
P = json.load(open('biane_polya.json'))
fc = P['fig_counting']
fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.8))
ax[0].plot(uu, k1, color=C1, lw=1.2, label=r'$\Phi(u)/\Phi(0)$')
ax[0].plot(uu, k2, '--', color=C2, lw=1.2, label=r'$\cosh\frac{9u}{2}\,e^{-2\pi(\cosh2u-1)}$')
ax[0].plot(uu, k3, ':', color=C3, lw=1.4, label=r'$e^{-2\pi(\cosh2u-1)}$')
ax[0].set_xlabel(r'$u$')
ax[0].set_title('kernels, normalised at $u=0$', fontsize=10)
ax[0].legend(fontsize=8, frameon=False)


def steps(z, T):
    xs, ys, n = [0.0], [0], 0
    for t in z:
        xs += [t, t]
        ys += [n, n + 1]
        n += 1
    xs.append(T)
    ys.append(n)
    return xs, ys


T = 100.0
for key, col, ls, lab in [('xi', C1, '-', r'$\Xi$'), ('xistar', C2, '--', r'$\Xi^*$'), ('P', C3, '-', r'$P$')]:
    xs, ys = steps(fc[key], T)
    ax[1].plot(xs, ys, ls, color=col, lw=1.0, label=lab)
tt = np.linspace(1, T, 400)
N0 = tt / TWO_PI * np.log(tt / (TWO_PI * math.e))
ax[1].plot(tt, N0 + 7 / 8, color='0.35', lw=0.8, ls='-.', label=r'$N_0(T)+\frac{7}{8}$')
ax[1].plot(tt, N0 - 1 / 4, color='0.6', lw=0.8, ls=':', label=r'$N_0(T)-\frac{1}{4}$')
ax[1].set_xlim(0, T)
ax[1].set_ylim(-1, 31)
ax[1].set_xlabel(r'$T$')
ax[1].set_title('zeros in $(0,T]$', fontsize=10)
ax[1].legend(fontsize=8, frameon=False, loc='upper left')
fig.tight_layout()
fig.savefig('fig/biane_kernels.pdf')
plt.close(fig)

# ---------------------------------------------------------------- 2. the Wick rotation of the Laplace variable
def recip_S(theta, k):
    q = 2 * theta - k
    sk = {1: math.sin(1.0), 0: 1.0, -1: math.sinh(1.0)}[k]
    if q > 0:
        r = math.sqrt(q)
        return math.sinh(r) / (r * sk)
    if q < 0:
        r = math.sqrt(-q)
        return math.sin(r) / (r * sk)
    return 1 / sk


th = np.linspace(-50, 4, 2000)
fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.8))
for k, col, ls in [(1, C2, '--'), (0, C1, '-'), (-1, C3, ':')]:
    y = [recip_S(t, k) for t in th]
    ax[0].plot(th, y, ls, color=col, lw=1.1, label=r'$\kappa=%d$' % k)
    atoms = [math.pi ** 2 * n * n / 2 - k / 2 for n in range(1, 4)]
    ax[0].plot([-a for a in atoms], [0, 0, 0], '|', color=col, ms=9, mew=1.2)
ax[0].axhline(0, color='0.5', lw=0.6)
ax[0].set_ylim(-0.35, 2.2)
ax[0].set_xlabel(r'$u$')
ax[0].set_title(r'$1/\mathbb{E}\,e^{-uJ^{(\kappa)}_1}$', fontsize=10)
ax[0].legend(fontsize=8, frameon=False, loc='upper left')
th2 = np.linspace(-30, 2, 2000)
z = np.exp(np.sqrt((1 + 2 * th2).astype(complex)) - 1)
ax[1].plot(th2, z.real, color=C1, lw=1.1, label='real part')
ax[1].plot(th2, np.abs(z), '--', color=C2, lw=1.1, label='modulus')
ax[1].axvline(-0.5, color='0.5', lw=0.6, ls=':')
ax[1].axhline(0, color='0.5', lw=0.6)
ax[1].set_ylim(-0.5, 2.2)
ax[1].set_xlabel(r'$u$')
ax[1].set_title(r'$1/\mathbb{E}\,e^{-uH^{(1)}_1}$', fontsize=10)
ax[1].legend(fontsize=8, frameon=False, loc='upper left')
fig.tight_layout()
fig.savefig('fig/biane_wick.pdf')
plt.close(fig)

# ---------------------------------------------------------------- 3. the clock
C = json.load(open('biane_clock.json'))['fig']
fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.8))
th = np.array(C['theta'])
ax[0].semilogy(th, C['LT_exact'], color=C1, lw=1.2, label='exact')
mc, se = np.array(C['LT_mc']), np.array(C['LT_se'])
ax[0].errorbar(th[::4], mc[::4], yerr=2 * se[::4], fmt='o', ms=3, color=C2, lw=0.8, label=r'Monte Carlo, $\pm2$ s.e.')
for nu2 in [95.428869, 154.96485, 215.64804]:
    ax[0].axvline(nu2, color='0.6', lw=0.6, ls=':')
ax[0].set_xlabel(r'$u$')
ax[0].set_title(r'$\mathbb{E}\,e^{-u\tau_r}$, $r=2\pi e^{2.5}$', fontsize=10)
ax[0].legend(fontsize=8, frameon=False)
edges = np.array(C['hist_edges'])
ax[1].bar(edges[:-1], C['hist'], width=np.diff(edges), align='edge', color=C1, alpha=0.35, label='Monte Carlo')
ax[1].plot(C['t'], C['dens_start'], color=C2, lw=1.2, label=r'exact, $\tau_r$')
ax[1].plot(C['t'], C['dens_oo'], '--', color=C3, lw=1.2, label=r'exact, $\tau_\infty$')
ax[1].set_xlabel(r'$x$')
ax[1].set_title('density of the clock', fontsize=10)
ax[1].legend(fontsize=8, frameon=False)
fig.tight_layout()
fig.savefig('fig/biane_clock.pdf')
plt.close(fig)

# ---------------------------------------------------------------- 4. the flow
F = json.load(open('biane_flow.json'))
fz = F['flow_zeros']
fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.9), sharey=True)
for a, key, title in [(ax[0], 'P', r'$P_t$'), (ax[1], 'X', r'$\Xi_t$')]:
    for lam, zs in fz[key].items():
        a.plot([float(lam)] * len(zs), zs, '.', color=C1, ms=2.2)
    a.set_title(title, fontsize=10)
    a.set_xlabel(r'$t$')
    a.set_ylim(0, 90)
    a.set_xlim(-61, 1)
ax[0].set_ylabel(r'real zeros $\theta$ in $(0,90)$')
fig.tight_layout()
fig.savefig('fig/biane_flow.pdf')
plt.close(fig)
gaps = []
for key in ['P', 'X']:
    for lam in ['0.0', '-10.0', '-30.0', '-60.0']:
        zs = fz[key][lam]
        if len(zs) < 2:
            gaps.append([key, lam, len(zs), None, None])
            continue
        mid = min((b - a) * math.log((a + b) / 2 / TWO_PI) / TWO_PI for a, b in zip(zs[:-1], zs[1:]))
        low = min((b - a) * math.log(a / TWO_PI) / TWO_PI for a, b in zip(zs[:-1], zs[1:]))
        gaps.append([key, lam, len(zs), round(mid, 4), round(low, 4)])
out['smallest normalised gap: function, lambda, real zeros in (0,90), density at the midpoint, density at the lower end'] = gaps
print(gaps)
speeds = {}
for key in ['P', 'X']:
    a, b = fz[key]['0.0'], fz[key]['-0.5']
    assert len(a) == len(b)
    sp = [(y - x) / 0.5 for x, y in zip(a, b)]
    i, j = sp.index(min(sp)), sp.index(max(sp))
    speeds[key] = {'zeros': len(sp), 'min': round(min(sp), 4), 'zero with min': round(a[i], 3), 'max': round(max(sp), 4),
                   'zero with max': round(a[j], 3), 'mean': round(sum(sp) / len(sp), 4), 'first zero': round(sp[0], 4)}
speeds['pi/8'] = round(math.pi / 8, 4)
out['mean speed of the real zeros in (0,90) over -0.5 <= lambda <= 0: (zero at -0.5 minus zero at 0)/0.5'] = speeds
print(speeds)

# ---------------------------------------------------------------- 5. complete monotonicity
Tr = json.load(open('biane_traces.json'))
cm = Tr['fig_cm']
fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.8))
TT = np.array(cm['T'])
for key, col, ls in [('100,1', C1, '-'), ('300,10', C2, '--')]:
    n_, v_ = key.split(',')
    ax[0].plot(TT, cm['profiles'][key], ls, color=col, lw=1.2, label=r'$n=%s,\ y=%s$' % (n_, v_))
ax[0].plot(cm['zeros'], np.zeros(len(cm['zeros'])), '|', color='0.4', ms=6)
ax[0].set_xlabel(r'ordinate $\gamma$')
ax[0].set_title(r'weight $\gamma^ne^{-\gamma y}$, normalised', fontsize=10)
ax[0].legend(fontsize=8, frameon=False)
for key, col, ls in [('50,0.25,1', C1, '-'), ('14.13,0.5,0.5', C2, '--')]:
    T_, d_, v_ = key.split(',')
    y = cm['sign'][key]
    ax[1].plot(np.arange(len(y)), y, ls, color=col, lw=1.1, label=r'$\gamma=%s,\ \delta=%s,\ y=%s$' % (T_, d_, v_))
ax[1].axhline(0, color='0.5', lw=0.6)
ax[1].set_xlabel(r'order $n$')
ax[1].set_title(r'sign factor $\cos(n\varphi-\delta y)$', fontsize=10)
ax[1].set_ylim(-1.1, 1.7)
ax[1].legend(fontsize=8, frameon=False, loc='upper center')
fig.tight_layout()
fig.savefig('fig/biane_cm.pdf')
plt.close(fig)

# ---------------------------------------------------------------- 6. the two traces
ft = Tr['fig_traces']
v = np.array(ft['v'])
K, KP = np.array(ft['Theta']), np.array(ft['ThetaP2'])
x = np.log(1 / (TWO_PI * v)) - EG
fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.8))
ax[0].plot(x, TWO_PI * v * K, color=C1, lw=1.2, label=r'$2\pi y\,K_{1/2}(y)$')
ax[0].plot(x, TWO_PI * v * KP, '--', color=C2, lw=1.2, label=r'$2\pi y\,K^P(y)$')
ax[0].plot(x, x, ':', color='0.4', lw=1.0, label=r'$\log\frac{1}{2\pi y}-\gamma_E$')
ax[0].set_xlabel(r'$\log\frac{1}{2\pi y}-\gamma_E$')
ax[0].legend(fontsize=8, frameon=False)
ax[0].set_title('the two traces', fontsize=10)
c = 2 * math.pi - 1 / (6 * math.pi) - 1 / (48 * math.pi)
ax[1].semilogx(v, K - KP, color=C1, lw=1.2, label=r'$D(y)=K_{1/2}(y)-K^P(y)$')
ax[1].semilogx(v, 9 / 8 - c * v * np.log(1 / v), '--', color=C2, lw=1.1, label=r'$\frac{9}{8}-c\,y\log\frac{1}{y}$')
ax[1].axhline(9 / 8, color='0.5', lw=0.8, ls=':')
ax[1].set_xlabel(r'$y$')
ax[1].set_title('their difference', fontsize=10)
ax[1].legend(fontsize=8, frameon=False, loc='lower left')
fig.tight_layout()
fig.savefig('fig/biane_traces.pdf')
plt.close(fig)
json.dump(out, open('biane_figs.json', 'w'), indent=1)
print('done')
