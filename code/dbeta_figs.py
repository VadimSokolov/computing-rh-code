# dbeta_figs.py: Figures fig:dbeta:thorin and fig:dbeta:arg of Section sec:dbeta:beta (ch/ch16.tex), the Dirichlet beta
# function L(s, chi_4) in Thorin form. Added for the book in October 2026, from the first author's note on the Dirichlet beta
# function in Thorin form (misc/reviews/nick-dbeta/, Figures 1 and 2 there), redrawn from an independent computation with
# Arb (python-flint) at 53 bits. xi(s, chi_4) = (4/pi)^{(s+1)/2} Gamma((s+1)/2) L(s, chi_4), v_alpha(t) = Re (xi'/xi)(alpha+it).
# dbeta_thorin.pdf. Left: v_alpha(t) = pi rho_{alpha-1/2}(t) for alpha = 1, 0.75, 0.55 on [0, 30], logarithmic scale, with the
# ordinates of the zeros (from dbeta_zeros.json). Right: at alpha = 3/2, the drift and archimedean piece
# b_alpha + G_alpha(t) = A^inf_4(alpha, t) + (L'/L)(alpha), the prime piece varpi_alpha(t) = Re (L'/L)(alpha+it) - (L'/L)(alpha), and v_alpha.
# dbeta_argument.pdf. Left: the cumulative Thorin mass (1/pi) arg xi(alpha+it, chi_4) for alpha = 1, 0.75, 0.55, continued
# from the real axis through the principal arguments of quotients of values at steps of 0.005 (the largest increment is
# printed, and asserted to be below 0.3 in modulus), against the zero count N(t). Right: v_{0.75}(t), its value b_{0.75} at
# t = 0, and the comparison bound (0.25/2.65)(A^inf_4(3.15, t) - S_4(3.15)) of Corollary cor:dbeta:lower with alpha' = 3.15.
# Run on Hopper in the folder that holds dbeta_zeros.json (written by dbeta_zeros.py): python3 dbeta_figs.py; writes
# fig/dbeta_thorin.pdf and fig/dbeta_argument.pdf there, and the extremes of the plotted curves to dbeta_figs.json.
import json
import math
import os

import matplotlib
import numpy as np
from flint import acb, acb_series, arb, ctx, dirichlet_char

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': 0.25, 'axes.spines.top': False,
                     'axes.spines.right': False})
ctx.prec = 53
CHI4 = dirichlet_char(4, 3)
PI = arb.pi()
LOG4PI = (arb(4) / PI).log()
Z = json.load(open('dbeta_zeros.json'))
ZEROS = [float(x) for x in Z['chi4']['ordinates']]
out = {}
os.makedirs('fig', exist_ok=True)


def dlogL(s):
    c = acb_series.dirichlet_l(acb_series([s, 1], prec=2), CHI4).coeffs()
    return c[1] / c[0]


def v(a, t):
    s = acb(a, t)
    return float((LOG4PI / 2 + ((s + 1) / 2).digamma().real / 2 + dlogL(s).real).mid())


def A(a, t):
    return float((LOG4PI / 2 + ((acb(a, t) + 1) / 2).digamma().real / 2).mid())


def Ssum(a):
    c = acb_series.zeta(acb_series([acb(a), 1], prec=2)).coeffs()
    return float((-(c[1] / c[0]).real - arb(2).log() / (arb(2) ** a - 1)).mid())


def xi4(s):
    return ((s + 1) / 2 * LOG4PI).exp() * ((s + 1) / 2).gamma() * acb.dirichlet_l(s, CHI4)


col = {1: '#1b4f72', 0.75: '#c0392b', 0.55: '#27ae60'}
ts = np.arange(0, 30.0001, 0.005)

# Figure 1
fig, ax = plt.subplots(1, 2, figsize=(9.8, 3.7))
for a in [1, 0.75, 0.55]:
    vv = [v(a, t) for t in ts]
    out['v_%s' % a] = {'min': min(vv), 'at': float(ts[int(np.argmin(vv))]), 'max': max(vv)}
    ax[0].semilogy(ts, vv, color=col[a], lw=1.0, label=r'$\alpha=%s$' % a)
for z in ZEROS:
    if z <= 30:
        ax[0].axvline(z, color='0.5', lw=0.5, ls='--')
ax[0].set_xlim(0, 30)
ax[0].set_xlabel(r'$\theta$')
ax[0].set_title(r'$v_\alpha(\theta)=\mathrm{Re}\,(\xi^\prime/\xi)(\alpha+i\theta,\chi_4)$', fontsize=10)
ax[0].legend(fontsize=8, loc='lower right')
a = 1.5
bb = float(dlogL(acb(a)).real.mid())
t2 = np.arange(0, 30.0001, 0.01)
P = [float(dlogL(acb(a, t)).real.mid()) - bb for t in t2]
CG = [A(a, t) + bb for t in t2]
VV = [p + c for p, c in zip(P, CG)]
out['decomposition_3/2'] = {'min P': min(P), 'at': float(t2[int(np.argmin(P))]), 'max P': max(P), 'min v': min(VV)}
ax[1].plot(t2, CG, color='#1b4f72', lw=1.2, label=r'drift and archimedean, $b_\alpha+G_\alpha$')
ax[1].plot(t2, P, color='#c0392b', lw=1.2, label=r'primes, $\varpi_\alpha$')
ax[1].plot(t2, VV, color='k', lw=1.2, label=r'sum, $v_\alpha$')
ax[1].axhline(0, color='k', lw=0.6)
ax[1].set_xlim(0, 30)
ax[1].set_xlabel(r'$\theta$')
ax[1].set_title(r'Sine squared pieces at $\alpha=3/2$', fontsize=10)
ax[1].legend(fontsize=8)
fig.tight_layout()
fig.savefig('fig/dbeta_thorin.pdf')

# Figure 2
fig, ax = plt.subplots(1, 2, figsize=(9.8, 3.7))
maxinc = 0.0
for a in [1, 0.75, 0.55]:
    prev = xi4(acb(a, 0))
    arg = 0.0
    M = [0.0]
    for t in ts[1:]:
        cur = xi4(acb(a, float(t)))
        inc = float((cur / prev).arg().mid())
        maxinc = max(maxinc, abs(inc))
        arg += inc
        M.append(arg / math.pi)
        prev = cur
    out['M_%s(30)' % a] = M[-1]
    ax[0].plot(ts, M, color=col[a], lw=1.1, label=r'$\alpha=%s$' % a)
out['largest arg increment'] = maxinc
assert maxinc < 0.3
N = [sum(1 for z in ZEROS if z <= t) for t in ts]
ax[0].step(ts, N, where='post', color='k', lw=0.9, label=r'$N(\theta)$')
ax[0].set_xlim(0, 30)
ax[0].set_xlabel(r'$\theta$')
ax[0].set_title(r'Cumulative Thorin mass $\pi^{-1}\arg\xi(\alpha+i\theta,\chi_4)$', fontsize=10)
ax[0].legend(fontsize=8, loc='upper left')
v75 = [v(0.75, t) for t in t2]
b75 = v75[0]
S315 = Ssum(3.15)
bound = [0.25 / 2.65 * (A(3.15, t) - S315) for t in t2]
cross = next(t for t, bd in zip(t2, bound) if bd > b75)
out['bound_0.75'] = {'b': b75, 'bound exceeds b from t': float(cross), 'bound at 30': bound[-1]}
ax[1].plot(t2, v75, color='k', lw=1.1, label=r'$v_{0.75}(\theta)$')
ax[1].axhline(b75, color='#1b4f72', lw=1.0, ls=':', label=r'$b_{0.75}=v_{0.75}(0)$')
ax[1].plot(t2, bound, color='#c0392b', lw=1.2, ls='--', label=r"comparison bound, $\alpha'=3.15$")
ax[1].set_xlim(0, 30)
ax[1].set_ylim(0, 0.6)
ax[1].set_xlabel(r'$\theta$')
ax[1].set_title(r'The minimum at the origin, $\alpha=0.75$', fontsize=10)
ax[1].legend(fontsize=8, loc='upper left')
fig.tight_layout()
fig.savefig('fig/dbeta_argument.pdf')
print(json.dumps(out, indent=1))
json.dump(out, open('dbeta_figs.json', 'w'), indent=1)
