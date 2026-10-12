# maxexp_figs.py: Figures fig:maxexp:laws, fig:maxexp:dominance, fig:maxexp:kernels and fig:maxexp:flow of Sections
# sec:ch16:max and sec:ch16:dbn (ch/ch16.tex, Sections 8.5.4 and 8.5.5). Added for the book in October 2026; they redraw, with
# corrections, Figures 1 to 4 of the first author's note on maxima of exponentials (misc/reviews/nick-maxexp/).
#  maxexp_laws.pdf: left, the Bondesson function N_1(z) of M_1 (Thorin atoms +1 at b_{6k-1}, b_{6k+1} and -1 at b_{6k}, with
#    b_n = pi (n^2 - 1)/12) and its asymptote sqrt(3z/pi)/3; right, the poles and zeros of psi*(theta) = E exp(-theta M*) on
#    Re theta >= -2 pi x 61, from maxexp_psizeros.json (the note's right panel, a counting function of real zeros only, is replaced).
#  maxexp_dominance.pdf: where |A(1-s)| > |A(s)| for Hecke's halves A(s) = int_1^oo Delta(iy) y^{s+9/2} dy, on 0 <= theta <= 25,
#    1/2 < sigma <= 2.5, computed here as in maxexp_hb.py (Gauss Legendre in v = log y, double precision), with the zeros on
#    the critical line.
#  maxexp_kernels.pdf: left, K_Delta(v) = e^{6v} Delta(i e^v) and Polya's kernel e^{-2 pi cosh v}, each divided by its value at
#    v = 0; right, log H(v) = (log F(e^v) + log F(e^{-v}))/2, the multiplier of K_Delta = e^{-2 pi cosh v} H.
#  maxexp_flow.pdf: the real zeros of H^Delta_t below 52 for -1.1 <= t <= 0.5 and the collisions, from maxexp_flow.json.
# Run from a folder holding maxexp_psizeros.json and maxexp_flow.json: python3 maxexp_figs.py; writes the four PDFs to fig/
# there (created if needed), as dbnplot.py does.
import json, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': 0.25, 'axes.spines.top': False, 'axes.spines.right': False})
C1, C2, C3 = '#1b4f72', '#c0392b', '#27ae60'
os.makedirs('fig', exist_ok=True)

# 1. the two laws
fig, ax = plt.subplots(1, 2, figsize=(9.8, 3.7))
b = lambda n: np.pi * (n * n - 1) / 12
atoms = sorted([(b(n), 1) for n in range(5, 200) if n % 6 in (1, 5)] + [(b(6 * k), -1) for k in range(1, 34)])
zmax = 1500.0
xs, ys, N = [0.0], [0], 0
for z, w in atoms:
    if z > zmax: break
    xs += [z, z]; ys += [N, N + w]; N += w
xs.append(zmax); ys.append(N)
ax[0].plot(xs, ys, color=C1, lw=1.1, label=r'$N_1(z)$')
zz = np.linspace(0, zmax, 400)
ax[0].plot(zz, np.sqrt(3 * zz / np.pi) / 3, '--', color=C2, lw=1.1, label=r'$\frac{1}{3}\sqrt{3z/\pi}$')
ax[0].set_xlabel(r'$z$'); ax[0].set_title(r'The Bondesson function of $M_1$', fontsize=10); ax[0].legend(fontsize=8, loc='upper left')
P = json.load(open('maxexp_psizeros.json'))
poles = np.arange(1, 61)
ax[1].plot(-poles, np.zeros(60), '|', color='0.35', ms=9, mew=1.0, label='poles')
rz = [float(x) / (2 * np.pi) for x in P['real_zeros']]
ax[1].plot(rz, np.zeros(len(rz)), 'o', color=C1, ms=4, label='real zeros')
nr = [(z['re'] / (2 * np.pi), z['im']) for z in P['nonreal_zeros_upper_half']]
ax[1].plot([x for x, y in nr] * 2, [y for x, y in nr] + [-y for x, y in nr], 'D', color=C2, ms=5, mfc='none', mew=1.2, label='nonreal zeros')
ax[1].set_xlim(-62, 1); ax[1].set_ylim(-4.5, 4.5)
ax[1].set_xlabel(r'$\mathrm{Re}\,\theta/2\pi$'); ax[1].set_ylabel(r'$\mathrm{Im}\,\theta$')
ax[1].set_title(r'Poles and zeros of $\psi^*(\theta)=\mathrm{E}\,e^{-\theta M^*}$', fontsize=10); ax[1].legend(fontsize=8, loc='upper right')
fig.tight_layout(); fig.savefig('fig/maxexp_laws.pdf'); plt.close(fig)

# kernel helpers (double precision)
def logF(y):
    y = np.asarray(y, dtype=float)
    small = y < 1
    yy = np.where(small, 1 / y, y)
    s = np.zeros_like(yy)
    for n in range(1, 60):
        s += 24 * np.log1p(-np.exp(-2 * np.pi * n * yy))
    return np.where(small, -12 * np.log(y) + 2 * np.pi * y - 2 * np.pi / y + s, s)
def K(v):
    v = np.abs(np.asarray(v, dtype=float))
    return np.exp(6 * v - 2 * np.pi * np.exp(v) + logF(np.exp(v)))

# 2. dominance map
X, WX = np.polynomial.legendre.leggauss(40)
def nodes(a, bb=4.5, h=0.25):
    m = int(np.ceil((bb - a) / h)); e = np.linspace(a, bb, m + 1)
    v = np.concatenate([(e[i] + e[i + 1]) / 2 + (e[i + 1] - e[i]) / 2 * X for i in range(m)])
    w = np.concatenate([(e[i + 1] - e[i]) / 2 * WX for i in range(m)])
    return v, w * K(v)
v1, k1 = nodes(0.0)
TH = np.linspace(0, 25, 501); SG = np.linspace(0.5, 2.5, 161)[1:]
E1 = np.exp(1j * TH[:, None] * v1[None, :]); E2 = np.exp(-1j * TH[:, None] * v1[None, :])
R = np.array([np.log(np.abs(E2 @ (k1 * np.exp(-(s - 0.5) * v1)))) - np.log(np.abs(E1 @ (k1 * np.exp((s - 0.5) * v1)))) for s in SG])
f0 = E1 @ k1
zl = [TH[i] - f0.real[i] * (TH[i + 1] - TH[i]) / (f0.real[i + 1] - f0.real[i]) for i in range(len(TH) - 1) if f0.real[i] * f0.real[i + 1] < 0]
fig, ax = plt.subplots(figsize=(7.5, 3.2))
ax.contourf(TH, SG, (R > 0).astype(float), levels=[0.5, 1.5], colors=['0.78'])
ax.contour(TH, SG, R, levels=[0], colors=[C1], linewidths=1.0)
ax.plot(zl, [0.5] * len(zl), 'o', color=C2, ms=6, clip_on=False, zorder=5, label=r'zeros of $\xi(\frac{1}{2}+i\theta,\Delta)$')
ax.set_xlim(0, 25); ax.set_ylim(0.5, 2.5)
ax.set_xlabel(r'$\theta$'); ax.set_ylabel(r'$\sigma$')
ax.set_title(r'Shaded: $|A(1-s)|>|A(s)|$, $s=\sigma+i\theta$', fontsize=10); ax.legend(fontsize=8, loc='upper left')
fig.tight_layout(); fig.savefig('fig/maxexp_dominance.pdf'); plt.close(fig)
print('dominance figure: zeros on the line', [round(z, 3) for z in zl])

# 3. kernels
v = np.linspace(-2.5, 2.5, 1001)
fig, ax = plt.subplots(1, 2, figsize=(9.8, 3.7))
ax[0].plot(v, K(v) / K(0.0), color=C1, lw=1.3, label=r'$K_\Delta(v)/K_\Delta(0)$')
ax[0].plot(v, np.exp(-2 * np.pi * (np.cosh(v) - 1)), '--', color=C2, lw=1.3, label=r'$e^{-2\pi(\cosh v-1)}$')
ax[0].set_xlabel(r'$v$'); ax[0].set_title('The two kernels', fontsize=10); ax[0].legend(fontsize=8)
logH = (logF(np.exp(v)) + logF(np.exp(-v))) / 2
ax[1].plot(v, logH, color=C3, lw=1.3)
ax[1].set_xlabel(r'$v$'); ax[1].set_title(r'$\log H(v)$, $H=K_\Delta(v)\,e^{2\pi\cosh v}$', fontsize=10)
fig.tight_layout(); fig.savefig('fig/maxexp_kernels.pdf'); plt.close(fig)
print('log H(0) =', float(logH[500]), ' check K(0) = Delta(i) =', float(K(0.0)))

# 4. the flow
Fl = json.load(open('maxexp_flow.json'))
T = sorted(Fl['trajectories'], key=float)
fig, ax = plt.subplots(figsize=(7.5, 3.9))
for ta, tb in zip(T, T[1:]):
    za, zb = np.array(Fl['trajectories'][ta]), np.array(Fl['trajectories'][tb])
    if len(za) == 0 or len(zb) == 0: continue
    for z in za:
        j = int(np.argmin(np.abs(zb - z)))
        if abs(zb[j] - z) < 0.6 and int(np.argmin(np.abs(za - zb[j]))) == int(np.argmin(np.abs(za - z))):
            ax.plot([float(ta), float(tb)], [z, zb[j]], color=C1, lw=0.9)
for e in Fl['events']:
    if e['event'] != 'collision': continue
    tc, xc = float(e['t_*']), float(e['theta_*'])
    above = min((t for t in T if float(t) >= tc), key=float)
    za = np.array(Fl['trajectories'][above])
    for z in sorted(za, key=lambda z: abs(z - xc))[:2]:
        ax.plot([float(above), tc], [z, xc], color=C1, lw=0.9)
    ax.plot(tc, xc, 'x', color=C2, ms=7, mew=1.6)
ax.axvline(0, color='0.4', lw=0.8, ls=':')
ax.set_xlim(-1.1, 0.5); ax.set_ylim(0, 54)
ax.set_xlabel(r'$t$'); ax.set_ylabel(r'zeros $\theta$ of $H^\Delta_t$')
ax.set_title(r'Real zeros below $52$ and their collisions ($\times$)', fontsize=10)
fig.tight_layout(); fig.savefig('fig/maxexp_flow.pdf'); plt.close(fig)
print('done')
