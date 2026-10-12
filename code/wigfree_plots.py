# wigfree_plots.py: the three figures of the free probability section (misc/reviews/nick-wigner/merge/wigfree.tex), drawn from
# the outputs of wigfree_jacobi.py, wigfree_upsilon.py and wigfree_heat.py. Added for the book in October 2026, from the first
# author's note on the semicircle law and the Thorin framework (misc/reviews/nick-wigner/source-v7/, version 7).
# fig/wigfree_laws.pdf (fig:wigfree:laws): left, the Thorin measure of the Upsilon image of the centre clock in the chart
#   theta = sqrt(u), whose distribution function is W(1/theta^2), against the zero count N(theta), the distribution function of
#   the Thorin measure of the clock itself; right, the densities of the clock T, of its Upsilon image and of its free
#   counterpart Lambda_BP(T) under the Bercovici Pata bijection.
# fig/wigfree_jacobi.pdf (fig:wigfree:jacobi): the Jacobi parameters alpha_n and 2 beta_n of the measure
#   nu = sum_j a_j^{-1} delta_{1/a_j}, n = 1..49, against half the right end L_n of the local arcsine interval from the Abel
#   inversion (with n(L_n) = n and with n(L) = n - 1/2) and against the cruder rule (1/2) gamma_{2n/pi}^{-2}.
# fig/wigfree_heat.pdf (fig:wigfree:heat): left, the real parts of the six roots of e^{-(t/2) D^2} p for the example of the
#   section, 0 <= t <= 0.4; right, the imaginary part of its complex pair, against sqrt(0.64 - t) for the pair on its own.
# Run from a folder that holds wigfree_jacobi.json, wigfree_upsilon.json and wigfree_heat.json: python3 wigfree_plots.py;
# writes the three PDF files there (a few seconds).
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': 0.25, 'axes.spines.top': False, 'axes.spines.right': False})

J = json.load(open('wigfree_jacobi.json'))
U = json.load(open('wigfree_upsilon.json'))
H = json.load(open('wigfree_heat.json'))

# ---------------- the three laws ----------------
fig, ax = plt.subplots(1, 2, figsize=(7.8, 3.3))
th = U['thorin of the Upsilon image, theta grid']
theta = np.array(th['theta'])
ax[0].plot(theta, th['N(theta)'], color='0.45', lw=1.0, drawstyle='steps-post', label=r'$N(\theta)$, Thorin measure of $T$')
ax[0].plot(theta, th['W(1/theta^2)'], color='C0', lw=1.4, label=r'$W(1/\theta^2)$, Thorin measure of $\Upsilon(T)$')
ax[0].set_xlabel(r'$\theta=\sqrt{u}$')
ax[0].set_title(r'Thorin distribution functions, mass in $(0,\theta^2]$', fontsize=10)
ax[0].legend(fontsize=8, frameon=False, loc='upper left')
ax[0].set_xlim(0, theta.max())
dens = U['densities']
for key, lab, col, ls in (('T', r'$T$', 'C0', '-'), ('Upsilon(T)', r'$\Upsilon(T)$', 'C1', '--'), ('Lambda(T)', r'$\Lambda_{\rm BP}(T)$', 'C2', '-.')):
    x = np.array(dens[key]['x'])
    f = np.array(dens[key]['f'])
    keep = x <= 0.1
    xx, ff = x[keep], f[keep]
    if key == 'Lambda(T)':
        # the free law vanishes to the left of its support [c, oo)
        xx = np.concatenate(([0.0, xx[0]], xx))
        ff = np.concatenate(([0.0, 0.0], ff))
    ax[1].plot(xx, ff, color=col, ls=ls, lw=1.3, label=lab)
ax[1].set_xlabel(r'$x$')
ax[1].set_xlim(0, 0.1)
ax[1].set_ylim(bottom=0)
ax[1].set_title('densities', fontsize=10)
ax[1].legend(fontsize=8, frameon=False, loc='upper right')
fig.tight_layout()
fig.savefig('wigfree_laws.pdf')
plt.close(fig)

# ---------------- the Jacobi parameters ----------------
S = J['series for the figure']
n = np.array(S['n'])
fig, ax = plt.subplots(figsize=(6, 3.4))
ax.loglog(n, np.array(S['alpha'], dtype=float), 'o', ms=3.2, color='C0', label=r'$\alpha_n$')
ax.loglog(n, np.array(S['2beta'], dtype=float), 'x', ms=4, color='C1', label=r'$2\beta_n$')
ax.loglog(n, np.array(S['abel'], dtype=float), '-', color='k', lw=1.1, label=r'$\frac{1}{2}L_n$, $n(L_n)=n$')
ax.loglog(n, np.array(S['abel_mid'], dtype=float), '--', color='0.4', lw=1.1, label=r"$\frac{1}{2}L_n'$, $n(L_n')=n-\frac{1}{2}$")
ax.loglog(n, np.array(S['crude'], dtype=float), ':', color='C3', lw=1.3, label=r'$\frac{1}{2}\gamma_{2n/\pi}^{-2}$')
ax.set_xlabel(r'$n$')
ax.set_title(r'Jacobi parameters of $\nu=\sum_ja_j^{-1}\delta_{1/a_j}$', fontsize=10)
ax.legend(fontsize=8, frameon=False, loc='lower left')
fig.tight_layout()
fig.savefig('wigfree_jacobi.pdf')
plt.close(fig)

# ---------------- the heat flow ----------------
fig, ax = plt.subplots(1, 2, figsize=(7.8, 3.3))
P = H['paths of the roots of e^{-(t/2) D^2} p']
t = np.array(P['t'])
re = np.array(P['re'])
im = np.array(P['im'])
for k in range(re.shape[1]):
    # plot each root as a point cloud in the (t, Re z) plane, filled where the root is real
    real = np.abs(im[:, k]) < 1e-7
    ax[0].plot(t[real], re[real, k], '.', ms=1.6, color='C0')
    ax[0].plot(t[~real], re[~real, k], '.', ms=1.6, color='C3')
tstar = float(H['least positive root t* of the discriminant (50 digits)'])
ax[0].axvline(tstar, color='0.5', lw=0.8, ls='--')
ax[0].set_xlabel(r'$t$')
ax[0].set_ylabel(r'$\mathrm{Re}\,z$')
ax[0].set_title(r'roots of $e^{-(t/2)D^2}p$ (red: complex pair)', fontsize=10)
# right: the imaginary part of the complex pair, against sqrt(0.64 - t) for the pair (x - 1)^2 + 0.64 on its own
pair = np.abs(im).max(axis=1)
ax[1].plot(t, pair, '-', color='C3', lw=1.4, label='the pair in $p$')
tt = np.linspace(0, 0.64, 400)
ax[1].plot(tt, np.sqrt(0.64 - tt), '--', color='0.4', lw=1.1, label=r'the pair alone, $\sqrt{0.64-t}$')
ax[1].axvline(tstar, color='0.5', lw=0.8, ls='--')
ax[1].set_xlabel(r'$t$')
ax[1].set_ylabel(r'$|\mathrm{Im}\,z|$')
ax[1].set_xlim(0, 0.66)
ax[1].set_ylim(bottom=0)
ax[1].set_title('imaginary part of the complex pair', fontsize=10)
ax[1].legend(fontsize=8, frameon=False, loc='upper right')
fig.tight_layout()
fig.savefig('wigfree_heat.pdf')
plt.close(fig)
print('ok')
