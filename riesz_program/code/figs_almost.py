# Added for the book (not part of the authors' package): draws Figures ag:fig:almost and ag:fig:cancel of Chapter
# ch:almost (fig/ag_almost.pdf and fig/ag_cancel.pdf) from almost.json, more.json and high.json, with the density written
# v_alpha/pi as in the book; the package has the notes' figures notes/fig_almost.pdf and notes/fig_cancel.pdf, which
# label it rho_alpha, but not the script that drew them, and the two new figures have the same content.
#   ag_almost, left: the zero density exponents of almost.py (checked against the table 'exp' of almost.json).
#   ag_almost, middle and right: v_alpha/pi at alpha = 0.6 and 0.7 on [284, 300], around the hypothetical zero at height
#       tj[4] = 292.40, by the formula of almost.py: the 1700 ordinates of g_*.npy on the critical line and the fifteen
#       quadruples at real parts 0.75 and 0.25 and heights tj of almost.json. almost.json keeps these curves only at
#       step 0.5, too coarse for the dips, so they are evaluated again at step 0.005, the grid of the notes' figure, and
#       compared with the archived values at the archived points.
#   ag_cancel, left: the negative masses of more.json (more2.py) against the basepoint.
#   ag_cancel, middle and right: the backgrounds of high.json (high.py) plus the kernels of one hypothetical pair, added
#       as in high.py, checked against the minima and negative masses stored in high.json.
# Run from code/ after setup.sh (it reads the JSON and npy files from the working directory): python3 figs_almost.py
# Writes ag_almost.pdf, ag_cancel.pdf and figs_almost.json (the checks and the depth of the dips at height 292.40).
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': 0.25, 'axes.spines.top': False, 'axes.spines.right': False})
BLUE, RED, ORANGE = '#1b4f72', '#c0392b', '#d68910'

A = json.load(open('almost.json'))
M = json.load(open('more.json'))
H = json.load(open('high.json'))
g = np.concatenate([np.load(f) for f in ['g_1_400.npy', 'g_401_1000.npy', 'g_1001_1700.npy']])
tj = A['synthetic']['tj']
off = [(0.75, t) for t in tj]
P = lambda a, x: a / (np.pi * (a * a + x * x))
out = {'n_ordinates': len(g), 'n_hypothetical_heights': len(tj)}


def rho(alpha, off, th):
    """v_alpha/pi on the grid th, exactly as rho() of almost.py (same terms, same order of summation)."""
    r = np.zeros_like(th)
    for gg in g:
        r += P(alpha - 0.5, th - gg) + P(alpha - 0.5, th + gg)
    for (b, t) in off:   # quadruple: b+it, 1-b+it (and the conjugates through +theta)
        r += P(alpha - b, th - t) + P(alpha - (1 - b), th - t) + P(alpha - b, th + t) + P(alpha - (1 - b), th + t)
    return r


# zero density exponents, as in almost.py
ingham = lambda s: 3 * (1 - s) / (2 - s)
gm = lambda s: 15 * (1 - s) / (3 + 5 * s)
dh = lambda s: 2 * (1 - s)
d = 0.0
for row in A['exp']:
    a = row['alpha']
    d = max(d, abs(ingham(a) - row['ingham']), abs(dh(a) - row['dh']))
    if row['gm'] is not None:
        d = max(d, abs(gm(a) - row['gm']))
out['exp_table_max_abs_diff'] = d

# the densities near the hypothetical zero at height tj[4], and the check against the archived curves (step 0.5)
t0 = tj[4]
th = np.linspace(284.0, 300.0, 3201)          # step 0.005
curves, out['almost_check'], out['dip'] = {}, {}, {}
for alpha in (0.6, 0.7):
    C = A['curve_%s' % alpha]
    tha = np.array(C['th'])
    out['almost_check'][str(alpha)] = dict(points=len(tha),
                                           max_abs_diff_with_zeros=float(np.max(np.abs(rho(alpha, off, tha) - np.array(C['r'])))),
                                           max_abs_diff_line_only=float(np.max(np.abs(rho(alpha, [], tha) - np.array(C['r0'])))))
    r0, r = rho(alpha, [], th), rho(alpha, off, th)
    curves[alpha] = (r0, r)
    i = int(np.argmin(r))
    thf = np.linspace(t0 - 0.05, t0 + 0.05, 100001)   # step 1e-6 around the hypothetical zero
    rf = rho(alpha, off, thf)
    k = int(np.argmin(rf))
    out['dip'][str(alpha)] = dict(height=t0, min_on_plot_grid=float(r[i]), at=float(th[i]), min_fine=float(rf[k]), at_fine=float(thf[k]),
                                  line_only_at_height=float(rho(alpha, [], np.array([t0]))[0]),
                                  negative_kernel_width=0.75 - alpha, negative_kernel_height=1 / (np.pi * (0.75 - alpha)))
    print('alpha %.1f: archived points agree to %.1e (with zeros) and %.1e (line only); min on the plot grid %.4f at %.3f, fine min %.4f at %.5f'
          % (alpha, out['almost_check'][str(alpha)]['max_abs_diff_with_zeros'], out['almost_check'][str(alpha)]['max_abs_diff_line_only'],
             r[i], th[i], rf[k], thf[k]))

fig, ax = plt.subplots(1, 3, figsize=(12.5, 3.5))
s = np.linspace(0.5, 1, 300)
s7 = np.linspace(0.7, 1, 181)
ax[0].plot(s, ingham(s), color=ORANGE, lw=1.5, label=r'Ingham $3(1-\alpha)/(2-\alpha)$')
ax[0].plot(s7, gm(s7), color=RED, lw=1.5, label=r'Guth and Maynard $15(1-\alpha)/(3+5\alpha)$')
ax[0].plot(s, dh(s), ':', color='gray', lw=1.5, label=r'density hypothesis $2(1-\alpha)$')
ax[0].axhline(1, color='k', lw=0.7)
ax[0].set_xlabel(r'$\alpha$'); ax[0].set_ylabel(r'exponent $e(\alpha)$'); ax[0].legend(fontsize=6.5, loc='lower left')
ax[0].set_title(r'$N(\alpha,T)\leq T^{e(\alpha)+o(1)}$: power saving $1-e(\alpha)$', fontsize=9)
for axk, alpha, col in [(ax[1], 0.6, BLUE), (ax[2], 0.7, RED)]:
    r0, r = curves[alpha]
    axk.plot(th, r0, color='gray', lw=0.8, label='zeros on the line only')
    axk.plot(th, r, color=col, lw=1.0, label=r'with hypothetical zeros at $\beta=0.75$')
    axk.axvline(t0, color='k', lw=0.5, ls=':')
    axk.axhline(0, color='k', lw=0.6)
    axk.set_xlabel(r'$\theta$'); axk.legend(fontsize=6.5, loc='upper left')
    axk.set_title(r'$v_\alpha/\pi$ at $\alpha=%s$ near a hypothetical zero' % alpha, fontsize=9)
fig.tight_layout(); fig.savefig('ag_almost.pdf'); plt.close(fig)

fig, ax = plt.subplots(1, 3, figsize=(12.5, 3.5))
for beta, col in [(0.6, RED), (0.75, ORANGE), (0.9, BLUE)]:
    rows = [c for c in M['curve'] if c['beta'] == beta]
    ax[0].plot([c['alpha'] for c in rows], [c['neg'] for c in rows], 'o-', color=col, lw=1.5, label=r'$\beta=%s$' % beta)
    ax[0].axvline(beta, color=col, lw=0.6, ls=':')
ax[0].axhline(len(tj), color='gray', lw=0.7, ls='--', label='number of hypothetical zeros')
ax[0].set_xlabel(r'basepoint $\alpha$'); ax[0].set_ylabel('negative mass on [0, 2000]'); ax[0].legend(fontsize=6.5, loc='upper left')
ax[0].set_title(r'Negative mass jumps to zero at $\alpha=\beta$', fontsize=9)
al = 0.6
out['high_check'] = {}
for axk, key, lab in [(ax[1], '300.0', '300'), (ax[2], '100000.0', '$10^5$')]:
    h, gam0 = H[key], float(key)
    tt, bg = np.array(h['th']), np.array(h['bg'])
    dth = tt[1] - tt[0]
    chk = dict(background_min=float(bg.min()), background_mean=float(bg.mean()), rows=[])
    for row in h['rows']:     # every hypothetical pair of high.py, to check the stored minima and negative masses
        w = row['beta'] - al
        r = bg + P(-w, tt - gam0) + P(al - (1 - row['beta']), tt - gam0)
        chk['rows'].append(dict(beta=row['beta'], min=float(r.min()), archived_min=row['min'],
                                neg=float(np.sum(np.clip(-r, 0, None)) * dth), archived_neg=row['neg']))
    chk['max_abs_diff'] = max(max(abs(c['min'] - c['archived_min']), abs(c['neg'] - c['archived_neg'])) for c in chk['rows'])
    out['high_check'][key] = chk
    print('height %s: background min %.4f mean %.4f; stored minima and negative masses agree to %.1e'
          % (key, chk['background_min'], chk['background_mean'], chk['max_abs_diff']))
    axk.plot(tt - gam0, bg, color='gray', lw=1.0, label=r'background, $\alpha=0.6$')
    for beta, col in [(0.7, RED), (0.9, BLUE)]:
        w = beta - al
        axk.plot(tt - gam0, bg + P(-w, tt - gam0) + P(al - (1 - beta), tt - gam0), color=col, lw=1.0, label=r'with a zero at $\beta=%s$' % beta)
    axk.axhline(0, color='k', lw=0.6); axk.set_ylim(-2, 4)
    axk.set_xlabel(r'$\theta-\gamma_0$'); axk.legend(fontsize=6.5, loc='lower left')
    axk.set_title(r'$v_{0.6}/\pi$ near a hypothetical zero, height %s' % lab, fontsize=9)
fig.tight_layout(); fig.savefig('ag_cancel.pdf'); plt.close(fig)

out['more_curve'] = M['curve']
json.dump(out, open('figs_almost.json', 'w'), indent=1)
print('ok')
