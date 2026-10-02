# NEW SCRIPT, not in the archive (round 2 item "kappa_cell", AUTHORS comment at ch/ch15.tex:24, review U10 B15).
# It computes the first order cell constant of Section sec:ch7:cumulative (ch/ch07.tex) and Section
# sec:ch15:first (ch/ch15.tex) from its definition,
#     kappa_cell = (1/pi) sum_{k=1}^{79} |l(b_k) - l(b_{k-1})|,   l = Xi'/Xi,
#     b_0 = 0, b_k = (gamma_k + gamma_{k+1})/2 for 1 <= k <= 78, b_79 = Theta = 200,
# which the book prints as 5.44 (ch01, ch03, ch07, ch15) and the AUTHORS comment as 5.43665, and which no script
# or data file of the archive holds. The other numbers of the paragraph ch15:22 and of the caption of
# fig:ch15:cells are recomputed beside it. Nothing here replaces an archived file; the closest archived code is
# the "mass derivative check" of code/plots.py, which evaluates l(200) by the Polya route and writes
# data/extra.json (XipXi200, pred_coeff).
#
# A  certified: l(b) = -Im xi'/xi(1/2+ib) and l''(b) = Im (log xi)'''(1/2+ib) in Arb ball arithmetic (400 bits),
#    xi'/xi(s) = 1/s + 1/(s-1) - log(pi)/2 + psi(s/2)/2 + zeta'/zeta(s), with the zeros gamma_1..gamma_80 from
#    acb.zeta_zeros (certified), at the exact cells and at the float64 cells of code/analyse.py.
# B  the book's own route: code/polya_flint.py with setup(420, 200, 3) and xi_pair, exactly as code/plots.py
#    computes XipXi200, at the float64 cells; this is the block proposed for code/plots.py (see REPORT.md).
# C  the truncated zero sum of Section sec:ch15:first: the float64 recipe that reproduces data/spacing.json
#    (coef 5.415, leaks), the same sum at high precision with Arb zeros, and the sum completed by the smooth
#    density theta'(t)/pi beyond the Gram point T* = g_1699 (the tail model of code/analyse.py and ch07).
# D  the exact cell distance d_cell(eps) = sum_k |mu_eps(C_k) - 1| from mu_eps(C_k) - 1 = (phi(b_k) - phi(b_{k-1}))/pi,
#    phi(b) = arg xi(1/2+eps+ib) - arg xi(1/2+ib) continued along the horizontal segment (proof of prop:mass),
#    at the eight resolutions of data/analysis.json and at 0.002 and 0.001, against tv_cells, cell_min, cell_max.
# Inputs: polya_flint.py (a copy of code/polya_flint.py) beside the script; the archived data files are read from
# /scratch/vsokolov/rh_book_repro/book_data and the float64 zero lists g_*.npy from the earlier rerun in
# /scratch/vsokolov/rh_book_repro/run (both read only). Output: kappa_cell.json. Run: python3 kappa_cell.py
import json, os, sys, time
import numpy as np
import mpmath as mp
from multiprocessing import Pool
from flint import arb, acb, acb_series, ctx

BOOK = '/scratch/vsokolov/rh_book_repro/book_data/'
RUN = '/scratch/vsokolov/rh_book_repro/run/'
THETA, K, PREC = 200, 79, 400
EPS_M = [1000, 500, 250, 100, 50, 20, 10, 5, 2, 1]   # eps = m/1000
T0 = time.time()


def log(*a):
    print('[%6.1fs]' % (time.time() - T0), *a, flush=True)


def s_of(b):
    return acb(arb(1) / 2, b)


def ell_pair(b):
    """(l(b), l''(b)) for l = Xi'/Xi, from the Taylor series of log xi at s = 1/2 + ib (constant terms dropped)."""
    S = acb_series([s_of(b), 1], prec=4)
    L = S.log() + (S - 1).log() - S * (arb.pi().log() / 2) + (S / 2).lgamma() + S.zeta().log()
    c = L.coeffs() + [acb(0)] * 4
    return -c[1].imag, 6 * c[3].imag     # l = i (log xi)', l'' = -i (log xi)''' (both real)


def s3(x, n=40):
    return x.str(n, radius=True)


def constants(ell, ell2):
    """cell errors e_k (coefficient of eps), cubic coefficients f_k (of eps^3), kappa and the eps^2 slope term."""
    pi = arb.pi()
    e = [-(ell[k] - ell[k - 1]) / pi for k in range(1, K + 1)]
    f = [(ell2[k] - ell2[k - 1]) / (6 * pi) for k in range(1, K + 1)]
    assert all(not v.contains(0) for v in e), 'a cell error is not separated from 0'
    kappa = sum((abs(v) for v in e), arb(0))
    c2 = sum((f[k] if e[k] > 0 else -f[k] for k in range(K)), arb(0))
    return e, f, kappa, c2


def phi_path(bf):
    """phi(b, eps) = arg xi(1/2+eps+ib) - arg xi(1/2+ib), continuous in eps, for eps = m/1000, m = 1..1000."""
    ctx.prec = 128
    b = arb(bf)
    s0 = s_of(b)
    g0 = s0.log().imag + (s0 - 1).log().imag + (s0 / 2).lgamma().imag
    zprev, acc, out, dmax = s0.zeta(), arb(0), {}, 0.0
    for m in range(1, 1001):
        s = acb(arb(1) / 2 + arb(m) / 1000, b)
        z = s.zeta()
        d = (z / zprev).arg()                          # principal argument of one small step
        dmax = max(dmax, abs(float(d.mid())))
        acc += d
        zprev = z
        if m in EPS_M:   # s(s-1) and Gamma(s/2) are continuous in the upper half plane; pi^{-s/2} adds nothing
            out[m] = (s.log().imag + (s - 1).log().imag + (s / 2).lgamma().imag - g0 + acc)
    return {m: v.str(30, radius=True) for m, v in out.items()}, dmax


if __name__ == '__main__':
    res = dict(note='kappa_cell of ch07/ch15 from its definition; see the header of kappa_cell.py')
    P1 = json.load(open(BOOK + 'part1.json'))
    SP = json.load(open(BOOK + 'spacing.json'))
    EX = json.load(open(BOOK + 'extra.json'))
    AN = json.load(open(BOOK + 'analysis.json'))

    # ---- A: certified evaluation at the exact cells ------------------------------------------------------
    ctx.prec = PREC
    zz = acb.zeta_zeros(1, 80)
    gam = [z.imag for z in zz]
    assert all(abs(float(z.real.mid()) - 0.5) < 1e-100 for z in zz)
    assert gam[78] < THETA < gam[79], 'N(200) is not 79'
    mp.mp.dps = 60
    ref = [mp.mpf(x) for x in P1['ref']]
    dref = max(abs(mp.mpf(g.mid().str(60, radius=False)) - r) for g, r in zip(gam, ref))
    log('Arb zeros gamma_1..gamma_80: gamma_79 =', gam[78].str(25), ', gamma_80 =', gam[79].str(25),
        '; max |Arb - part1.json ref| =', mp.nstr(dref, 3))
    b_ex = [arb(0)] + [(gam[k] + gam[k + 1]) / 2 for k in range(78)] + [arb(THETA)]
    L_ex = [ell_pair(b) for b in b_ex]
    ell = [p[0] for p in L_ex]; ell2 = [p[1] for p in L_ex]
    assert ell[0].contains(0) and ell2[0].contains(0)
    ell[0] = arb(0); ell2[0] = arb(0)                 # Xi is even: l(0) = l''(0) = 0
    e, f, kappa, c2 = constants(ell, ell2)
    tele = sum(e, arb(0))
    log('A exact cells: kappa_cell =', s3(kappa, 45))
    log('A sum of cell errors =', s3(tele, 30), '; -l(200)/pi =', s3(-ell[K] / arb.pi(), 30),
        '; l(200) =', s3(ell[K], 30))
    log('A eps^2 coefficient of d_cell/eps: c2 =', s3(c2, 20))
    ef = [float(v.mid()) for v in e]
    order = sorted(range(K), key=lambda k: -ef[k])
    res['A_exact_cells'] = dict(
        prec_bits=PREC, max_abs_diff_zeros_vs_part1_ref=float(dref),
        kappa_cell=s3(kappa, 45), kappa_cell_float=float(kappa.mid()),
        sum_cell_errors=s3(tele, 30), minus_l200_over_pi=s3(-ell[K] / arb.pi(), 30), l200=s3(ell[K], 30),
        c2_eps2_coefficient=s3(c2, 20), c2_float=float(c2.mid()),
        cell_errors=ef, cell_errors_cubic=[float(v.mid()) for v in f],
        l_at_boundaries=[float(v.mid()) for v in ell], l2_at_boundaries=[float(v.mid()) for v in ell2],
        boundaries=[b.str(30, radius=False) for b in b_ex],
        min=min(ef), argmin=int(np.argmin(ef)) + 1, max=max(ef), argmax=int(np.argmax(ef)) + 1,
        largest_five=[[k + 1, ef[k]] for k in order[:5]], n_positive=int(sum(v > 0 for v in ef)),
        max_radius_cell_error=max(float(v.rad()) for v in e))

    # ---- A2: the float64 cells of code/analyse.py --------------------------------------------------------
    near_hp = np.array([np.longdouble(s) for s in P1['ref']])[:80]
    cells = np.concatenate([[0.0], [float((near_hp[k] + near_hp[k + 1]) / 2) for k in range(78)], [THETA]])
    L64 = [ell_pair(arb(float(c))) for c in cells[1:]]
    ell64 = [arb(0)] + [p[0] for p in L64]; ell64_2 = [arb(0)] + [p[1] for p in L64]
    e64, f64, kappa64, c2_64 = constants(ell64, ell64_2)
    dcell = max(abs(float((arb(float(cells[k])) - b_ex[k]).mid())) for k in range(K + 1))
    log('A2 float64 cells of analyse.py: kappa_cell =', s3(kappa64, 30), '; max |cell - exact| =', '%.2e' % dcell,
        '; kappa difference', '%.2e' % float((kappa64 - kappa).mid()))
    res['A2_float64_cells'] = dict(kappa_cell=s3(kappa64, 30), max_abs_cell_shift=dcell,
                                   kappa_minus_exact=float((kappa64 - kappa).mid()), cells=[float(c) for c in cells])

    # ---- B: the Polya route of the book (code/polya_flint.py, as code/plots.py) -----------------------------
    sys.path.insert(0, os.getcwd())
    from polya_flint import setup, xi_pair
    t1 = time.time()
    S = setup(420, 200, 3)                             # as code/plots.py
    Z = np.load(RUN + 'plotcurves.npz')['zeros']       # float64 of the 80 zeros, as code/plots.py reads them
    # ---------- block proposed for code/plots.py, after out=dict(XipXi200=...) ----------
    bc = [0.0] + [float((Z[k] + Z[k + 1]) / 2) for k in range(78)] + [200.0]
    Lp = []
    for t in bc:
        xa, xb = xi_pair(S, acb(0, t)); Lp.append(float((-(xb.imag)) / xa.real))
    leaks = [-(Lp[k + 1] - Lp[k]) / np.pi for k in range(79)]
    kappa_plots = float(np.sum(np.abs(leaks)))
    # ------------------------------------------------------------------------------------
    LpA = []
    for t in bc:                                      # same route kept in balls, to compare digits
        xa, xb = xi_pair(S, acb(0, t)); LpA.append(-(xb.imag) / xa.real)
    ctx.prec = PREC
    ep = [-(LpA[k] - LpA[k - 1]) / arb.pi() for k in range(1, K + 1)]
    kappa_polya = sum((abs(v) for v in ep), arb(0))
    LA = [arb(0)] + [ell_pair(arb(float(t)))[0] for t in bc[1:]]
    dmax_polya = max(abs(float((LpA[k] - LA[k]).mid())) for k in range(1, K + 1))
    log('B Polya route (setup(420,200,3)) at the float64 cells of plots.py: kappa_cell =', s3(kappa_polya, 30),
        '; float64 block gives', repr(kappa_plots), '; l(200) =', repr(Lp[-1]), 'vs extra.json', EX['XipXi200'],
        '; max |l_Polya - l_Arb| =', '%.2e' % dmax_polya, '(%.1fs)' % (time.time() - t1))
    res['B_polya_route'] = dict(kappa_cell_balls=s3(kappa_polya, 30), kappa_cell_plots_block=kappa_plots,
                                l200_plots_block=Lp[-1], extra_json_XipXi200=EX['XipXi200'],
                                max_abs_diff_l_polya_vs_arb=dmax_polya,
                                cell_leaks_plots_block=leaks, sum_leaks_plots_block=float(np.sum(leaks)))

    # ---- C: truncated zero sums --------------------------------------------------------------------------
    g = np.concatenate([np.load(RUN + f) for f in ['g_1_400.npy', 'g_401_1000.npy', 'g_1001_1700.npy']])
    gw = g[g < THETA]; assert len(gw) == K
    bt = np.concatenate([[0.0], (gw[:-1] + gw[1:]) / 2, [THETA]])
    A = np.array([np.sum(1 / (c - g)) for c in bt]); C = np.array([np.sum(1 / (c + g)) for c in bt])
    lk = (A[:-1] - A[1:]) / np.pi + (C[:-1] - C[1:]) / np.pi
    same = int(sum(float(x) == float(y) for x, y in zip(lk, SP['leaks'])))
    coef = float(np.sum(np.abs(lk)))
    log('C float64 recipe: coef', repr(coef), 'vs spacing.json', repr(SP['coef']), '; leaks identical', same, 'of 79')
    srt = sorted(range(K), key=lambda k: -lk[k])
    res['C_truncated_float64'] = dict(coef=coef, coef_identical_to_spacing_json=coef == SP['coef'],
                                      leaks_identical_to_spacing_json=same, sum_leaks=float(np.sum(lk)),
                                      min=float(lk.min()), argmin=int(np.argmin(lk)) + 1, max=float(lk.max()),
                                      argmax=int(np.argmax(lk)) + 1,
                                      largest_five=[[k + 1, float(lk[k])] for k in srt[:5]],
                                      n_positive=int(np.sum(lk > 0)))
    # the same sum at high precision with the Arb zeros (exact cells), and the smooth tail beyond T* = g_1699
    ctx.prec = 128
    z17 = acb.zeta_zeros(1, 1700)
    mp.mp.dps = 40
    G = [mp.mpf(z.imag.mid().str(40, radius=False)) for z in z17]
    assert max(abs(G[k] - mp.mpf(float(g[k]))) for k in range(1700)) < 1e-12
    B = [mp.mpf(b.mid().str(40, radius=False)) for b in b_ex]
    H = [mp.fsum(2 * c / (c * c - t * t) for t in G) for c in B]          # 1/(c-t) + 1/(c+t)
    H2 = [mp.fsum(2 / (c - t) ** 3 + 2 / (c + t) ** 3 for t in G) for c in B]   # l'' under RH, truncated
    ell2m = [mp.mpf(v.mid().str(40, radius=False)) for v in ell2]
    dH2 = max(abs(H2[k] - ell2m[k]) for k in range(K + 1))
    log('C check of the series for l\'\': max |sum_{j<=1700} 2/(b-g)^3+2/(b+g)^3 - l\'\'(b)| =', mp.nstr(dH2, 3),
        '; max |l\'\'| =', mp.nstr(max(abs(v) for v in ell2m), 5))
    kap_tr = mp.fsum(abs(H[k] - H[k - 1]) for k in range(1, K + 1)) / mp.pi
    mp.mp.dps = 30
    Tstar = mp.findroot(lambda T: mp.siegeltheta(T) / mp.pi + 1 - 1700, 2197.5)
    assert G[1699] < Tstar < mp.mpf(acb.zeta_zero(1701).imag.mid().str(30, radius=False))
    mp.mp.dps = 40
    tail = [mp.quad(lambda t: 2 * c / (c * c - t * t) * mp.siegeltheta(t, 1) / mp.pi,
                    [Tstar, 2 * Tstar, 4 * Tstar, 16 * Tstar, mp.inf]) for c in B]
    Ht = [H[k] + tail[k] for k in range(K + 1)]
    kap_tail = mp.fsum(abs(Ht[k] - Ht[k - 1]) for k in range(1, K + 1)) / mp.pi
    ellm = [mp.mpf(v.mid().str(40, radius=False)) for v in ell]
    dH = max(abs(H[k] + tail[k] - ellm[k]) for k in range(1, K + 1))
    kap_m = mp.mpf(kappa.mid().str(40, radius=False))
    log('C high precision truncated sum (1700 Arb zeros, exact cells): kappa_1700 =', mp.nstr(kap_tr, 15),
        '; with the smooth tail beyond T* =', mp.nstr(Tstar, 12), ':', mp.nstr(kap_tail, 15),
        '; kappa - kappa_1700 =', mp.nstr(kap_m - kap_tr, 8), '(%.4f percent of kappa)' % float(100 * (kap_m - kap_tr) / kap_m),
        '; max |truncated + tail - l| =', mp.nstr(dH, 3))
    res['C_truncated_hp'] = dict(kappa_1700=mp.nstr(kap_tr, 20), kappa_1700_plus_tail=mp.nstr(kap_tail, 20),
                                 Tstar=mp.nstr(Tstar, 15), kappa_minus_kappa_1700=mp.nstr(kap_m - kap_tr, 12),
                                 percent_of_kappa=float(100 * (kap_m - kap_tr) / kap_m),
                                 percent_of_kappa_1700=float(100 * (kap_m - kap_tr) / kap_tr),
                                 kappa_minus_kappa_1700_plus_tail=mp.nstr(kap_m - kap_tail, 5),
                                 max_abs_l_minus_truncated_plus_tail=mp.nstr(dH, 5),
                                 max_abs_l2_minus_truncated_l2=mp.nstr(dH2, 5),
                                 sum_cell_errors_1700=mp.nstr(-(H[K] - H[0]) / mp.pi, 15),
                                 tail_at_boundaries=[float(v) for v in tail])

    # ---- D: exact cell distance at finite resolution -----------------------------------------------------
    t1 = time.time()
    with Pool(int(os.environ.get('SLURM_CPUS_PER_TASK', '8'))) as pool:
        paths = pool.map(phi_path, [float(c) for c in cells[1:]])
    ctx.prec = 128
    dstep = max(p[1] for p in paths)
    rows = []
    kf = float(kappa.mid()); c2f = float(c2.mid())
    for m in EPS_M:
        eps = m / 1000
        phi = [arb(0)] + [arb(p[0][m]) for p in paths]
        mu1 = [(phi[k] - phi[k - 1]) / arb.pi() for k in range(1, K + 1)]
        d = sum((abs(v) for v in mu1), arb(0))
        dm = [float(v.mid()) for v in mu1]
        mass = K + float((phi[K] / arb.pi()).mid())    # M_eps(200) = N(200) + phi(200)/pi
        row = dict(eps=eps, d_cell=float(d.mid()), d_cell_rad=float(d.rad()), slope=float(d.mid()) / eps,
                   slope_minus_kappa=float(d.mid()) / eps - kf, first_two_terms=kf + c2f * eps * eps,
                   cell_min=1 + min(dm), cell_max=1 + max(dm), mass=mass,
                   cell_errors=[v / eps for v in dm] if m in (5, 1) else None)
        a = [r for r in AN['rows'] if abs(r['eps'] - eps) < 1e-12]
        if a:
            row.update(archived_tv_cells=a[0]['tv_cells'], archived_slope=a[0]['tv_cells'] / eps,
                       archived_minus_exact=a[0]['tv_cells'] - row['d_cell'],
                       archived_cell_min=a[0]['cell_min'], archived_cell_max=a[0]['cell_max'],
                       archived_mass=a[0]['mass'], archived_mass_minus_exact=a[0]['mass'] - mass)
        rows.append(row)
        log('D eps %-6g d_cell %.12f slope %.9f slope-kappa %+.3e (c2 eps^2: %+.3e) cells [%.7f, %.7f] mass %.12f%s' % (
            eps, row['d_cell'], row['slope'], row['slope_minus_kappa'], c2f * eps * eps, row['cell_min'],
            row['cell_max'], mass, '' if not a else
            '; archived %.12f (diff %+.2e), cells [%.7f, %.7f], mass diff %+.1e' % (
                a[0]['tv_cells'], row['archived_minus_exact'], a[0]['cell_min'], a[0]['cell_max'],
                a[0]['mass'] - mass)))
    r05 = [r for r in AN['rows'] if r['eps'] == 0.005][0]['tv_cells'] / 0.005
    r10 = [r for r in AN['rows'] if r['eps'] == 0.01][0]['tv_cells'] / 0.01
    rich = (4 * r05 - r10) / 3
    log('D max arg step %.3f rad; Richardson (eps 0.01, 0.005) of the archived slopes: %.9f (kappa %.9f) (%.1fs)' % (
        dstep, rich, kf, time.time() - t1))
    res['D_finite_eps'] = dict(rows=rows, max_arg_step=dstep, richardson_archived=rich,
                               note='mu_eps(C_k)-1 = (phi(b_k)-phi(b_{k-1}))/pi, float64 cells of code/analyse.py')
    json.dump(res, open('kappa_cell.json', 'w'), indent=1)
    log('done')
