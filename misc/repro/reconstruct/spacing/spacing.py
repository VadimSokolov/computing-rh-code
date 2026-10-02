# RECONSTRUCTED SCRIPT, NOT THE AUTHORS' SCRIPT. No script in the archive writes data/spacing.json; this one
# is meant to reproduce it. Chapter ch:spacing (ch/ch15.tex) quotes the file, and code/plots3.py (part 9) draws
# Figure fig:ch15:cells from its leaks. Reconstructed in misc/repro/reconstruct/spacing/ from ch/ch15.tex,
# code/plots3.py, code/analyse.py and the round 1 review U10 (misc/reviews/round1/U10-review.md, U10-verify.md).
#
# spacing.json holds, in the archived order:
#   coef    kappa_cell of Section sec:ch15:first with the zero sum truncated to the first 1700 ordinates,
#           the sum of |leaks| (5.415 in ch15);
#   mean, var, min   of the 1699 unfolded spacings x_{k+1}-x_k, x_k = vartheta(gamma_k)/pi with vartheta the
#           Riemann Siegel theta function (0.99998, 0.146 and 0.089 in Section sec:ch15:spacings); var is the
#           population variance (numpy's default, ddof=0);
#   leaks   the 79 first order cell errors (mu_eps(C_k)-1)/eps, k=1..79, with Theta=200, K=N(200)=79, b_0=0,
#           b_k=(gamma_k+gamma_{k+1})/2 for 1<=k<79 and b_79=Theta (the cells of code/analyse.py):
#           leaks_k = -(1/pi) sum_{j<=1700} [1/(b_k-g_j)-1/(b_{k-1}-g_j)+1/(b_k+g_j)-1/(b_{k-1}+g_j)],
#           evaluated as [A(b_{k-1})-A(b_k)]/pi + [C(b_{k-1})-C(b_k)]/pi with A(c) = sum_j 1/(c-g_j) (kernels at
#           the ordinates) and C(c) = sum_j 1/(c+g_j) (kernels at their reflections).
#
# Conventions of the archive: the ordinates are float64 in g_1_400.npy, g_401_1000.npy and g_1001_1700.npy
# (as in code/analyse.py and code/plots3.py); if these files are absent the script regenerates them with
# mpmath.zetazero at 30 digits and saves them (as misc/repro/zeros.py does). vartheta is mpmath.siegeltheta at
# 30 digits (as in code/analyse.py and code/kenttail.py), rounded to float64 and divided by np.pi; the sums are
# float64 numpy sums. With the zero lists of misc/repro/zeros.py these float64 recipes reproduce all 83 numbers
# of data/spacing.json bit for bit (variants.py, variants2.py and variants3.py show which recipes do not).
# Also written, and not part of the archived file: spacing_extra.json (checks of the ordinates against Arb, the
# same quantities at 30 digits, and the Xi'/Xi form of kappa_cell asked for at ch15:24), and spacing.pdf, a
# candidate for fig/spacing.pdf (Figure fig:ch15:spacing), which has no drawing script in the archive either.
# Run: python3 spacing.py   (in a folder holding the three .npy files, or none of them)
import json, os, time
import numpy as np, mpmath as mp
from multiprocessing import Pool

FILES = ['g_1_400.npy', 'g_401_1000.npy', 'g_1001_1700.npy']
NZ, THETA = 1700, 200.0


def zetazero30(n):
    mp.mp.dps = 30
    return float(mp.zetazero(n).imag)


def ordinates():
    if all(os.path.exists(f) for f in FILES):
        return np.concatenate([np.load(f) for f in FILES]), 'loaded from ' + ', '.join(FILES)
    with Pool(int(os.environ.get('SLURM_CPUS_PER_TASK', os.cpu_count()))) as p:
        g = np.array(p.map(zetazero30, range(1, NZ + 1), chunksize=10))
    np.save(FILES[0], g[:400]); np.save(FILES[1], g[400:1000]); np.save(FILES[2], g[1000:1700])
    return g, 'computed with mpmath.zetazero at 30 digits and saved'


if __name__ == '__main__':
    t0 = time.time()
    g, src = ordinates()
    assert len(g) == NZ and np.all(np.diff(g) > 0)
    print('ordinates:', src, '; gamma_1 =', g[0], ', gamma_1700 =', g[-1], flush=True)

    # unfolded spacings (Section sec:ch15:spacings)
    mp.mp.dps = 30
    th = np.array([float(mp.siegeltheta(t)) for t in g])
    x = th / np.pi
    s = np.diff(x)

    # first order cell errors from the truncated zero sum (Section sec:ch15:first)
    gw = g[g < THETA]
    assert len(gw) == 79
    b = np.concatenate([[0.0], (gw[:-1] + gw[1:]) / 2, [THETA]])
    A = np.array([np.sum(1 / (c - g)) for c in b])    # A+C = Xi'/Xi(b) under RH, truncated at gamma_1700
    C = np.array([np.sum(1 / (c + g)) for c in b])
    leaks = (A[:-1] - A[1:]) / np.pi + (C[:-1] - C[1:]) / np.pi
    out = dict(coef=float(np.sum(np.abs(leaks))), mean=float(np.mean(s)), var=float(np.var(s)), min=float(np.min(s)),
               leaks=[float(v) for v in leaks])
    json.dump(out, open('spacing.json', 'w'))
    print('spacing.json: coef %.12g mean %.12g var %.12g min %.12g, %d leaks' % (out['coef'], out['mean'], out['var'],
          out['min'], len(leaks)), flush=True)

    # checks and additions, not in the archived file
    from flint import acb, arb, ctx
    ctx.prec = 128
    za = np.array([float(z.imag.mid()) for z in acb.zeta_zeros(1, NZ)])
    tha = np.array([float((acb(0.25, float(t) / 2).lgamma().imag - arb(float(t)) / 2 * arb.pi().log()).mid())
                    for t in g])                                       # vartheta = Im log Gamma(1/4+it/2)-(t/2)log pi
    G = [mp.mpf(float(t)) for t in g]                                  # the float64 ordinates, exactly
    X = [mp.siegeltheta(t) / mp.pi for t in G]
    S = [X[k + 1] - X[k] for k in range(NZ - 1)]
    m30 = mp.fsum(S) / len(S)
    B = [mp.mpf(0)] + [(G[k] + G[k + 1]) / 2 for k in range(78)] + [mp.mpf(THETA)]
    Lt = [mp.fsum(1 / (c - t) + 1 / (c + t) for t in G) for c in B]
    lk30 = [-(Lt[k + 1] - Lt[k]) / mp.pi for k in range(79)]

    def xi_log_der(c):   # Xi'/Xi(c) = -Im xi'/xi(1/2+ic), xi'/xi(s) = 1/s+1/(s-1)-log(pi)/2+psi(s/2)/2+zeta'/zeta(s)
        z = mp.mpc(0.5, c)
        return -mp.im(1 / z + 1 / (z - 1) - mp.log(mp.pi) / 2 + mp.digamma(z / 2) / 2 + mp.zeta(z, 1, 1) / mp.zeta(z))
    R = [xi_log_der(c) for c in B]
    lkxi = [-(R[k + 1] - R[k]) / mp.pi for k in range(79)]
    k = int(np.argmin(s))
    extra = dict(
        note='reconstruction checks, not part of data/spacing.json; indices are 1 based',
        ordinates=src, gamma_1=float(g[0]), gamma_79=float(g[78]), gamma_80=float(g[79]), gamma_1700=float(g[-1]),
        arb_zeros_identical_float64=int(np.sum(za == g)), arb_zeros_max_abs_diff=float(np.max(np.abs(za - g))),
        arb_theta_identical_float64=int(np.sum(tha == th)), arb_theta_max_abs_diff=float(np.max(np.abs(tha - th))),
        n_spacings=int(len(s)), var_ddof1=float(np.var(s, ddof=1)),
        min_pair_index=[k + 1, k + 2], min_pair_ordinates=[float(g[k]), float(g[k + 1])],
        N_to_minus_third=float(len(s) ** (-1 / 3)),
        mp30=dict(coef=float(mp.fsum(abs(v) for v in lk30)), mean=float(m30),
                  var=float(mp.fsum((v - m30) ** 2 for v in S) / len(S)), min=float(min(S)),
                  leaks=[float(v) for v in lk30]),
        sum_leaks=float(np.sum(leaks)), leaks_min=float(leaks.min()), leaks_argmin=int(np.argmin(leaks)) + 1,
        leaks_max=float(leaks.max()), leaks_argmax=int(np.argmax(leaks)) + 1,
        kappa_cell_xi=float(mp.fsum(abs(v) for v in lkxi)), sum_leaks_xi=float(mp.fsum(lkxi)),
        minus_xi_log_der_200_over_pi=float(-R[-1] / mp.pi), leaks_xi=[float(v) for v in lkxi],
        leaks_xi_min=float(min(lkxi)), leaks_xi_max=float(max(lkxi)))
    json.dump(extra, open('spacing_extra.json', 'w'), indent=1)
    print('Arb zeros: %d of %d identical in float64, max abs diff %.3g; Arb theta: %d identical, max abs diff %.3g'
          % (extra['arb_zeros_identical_float64'], NZ, extra['arb_zeros_max_abs_diff'],
             extra['arb_theta_identical_float64'], extra['arb_theta_max_abs_diff']))
    print('smallest spacing %.10f between gamma_%d = %.4f and gamma_%d = %.4f' % (s[k], k + 1, g[k], k + 2, g[k + 1]))
    print('sum of leaks %.7f; kappa_cell from Xi\'/Xi %.7f; -Xi\'(200)/(pi Xi(200)) %.7f'
          % (extra['sum_leaks'], extra['kappa_cell_xi'], extra['minus_xi_log_der_200_over_pi']), flush=True)

    # candidate for fig/spacing.pdf (Figure fig:ch15:spacing): 40 bins over the range of the spacings; the style
    # (size, colours, widths, grid, fonts) was read from the operators of the book's PDF (matplotlib 3.10.8)
    import matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 3.4))
    cnt, edges, _ = ax.hist(s, bins=40, density=True, color='lightsteelblue', edgecolor='white',
                            label='unfolded spacings, 1700 zeros')
    u = np.linspace(0, 3, 400)
    ax.plot(u, 32 / np.pi ** 2 * u ** 2 * np.exp(-4 * u ** 2 / np.pi), color='crimson', lw=1.5, label='GUE Wigner surmise')
    ax.plot(u, np.exp(-u), '--', color='gray', lw=1.5, label='Poisson')
    ax.set_xlabel('normalised spacing'); ax.set_title('Spacings of the first 1700 zeros'); ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig('spacing.pdf'); plt.close()
    json.dump(dict(edges=[float(e) for e in edges], density=[float(c) for c in cnt],
                   counts=[int(c) for c in np.histogram(s, bins=edges)[0]]), open('spacing_hist.json', 'w'))
    print('done in %.1f s' % (time.time() - t0))
