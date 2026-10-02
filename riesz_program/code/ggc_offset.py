# Added for the book (not part of the authors' package): splits the relative error of every row of the Thorin form in Table ri:tab:ggccoef (ch/riesz.tex, computed by ggc_coef.py, output ggc_coef.json) into the error of the Weyl tail beyond theta = 200 and the quadrature error of ggc_coef.py on [0, 200].
#
# With w = (2k-1)^2, rho = v_1/pi and v_1(theta) = Re (xi'/xi)(1 + i theta), the Thorin integral
#   I(w) = int_0^inf log(1 + w/theta^2) rho(theta) dtheta
# equals -log c_k exactly, c_k = pi^k/(2k(2k-1)(k-1)! zeta(2k)). The table prints kappa_k e^{-I} (Thorin form, exact value
# 1/zeta(2k)) and e^{-I}/e^{-b_1(2k-1)} (GGC factor, exact value c_k e^{b_1(2k-1)}), kappa_k = 2k(2k-1)(k-1)!/pi^k.
# ggc_coef.py truncates I(w): rho on [0, 200] by the trapezoid rule on its grid (step 0.025 from 0.025 to 60, one step
# 0.05 to 60.05, then step 0.1000357 to 200; the constant b_1/pi on [0, 0.025]) and the Weyl density (1/2pi) log(theta/2pi)
# beyond 200. This script computes, for k = 1..6,
#   exact      c_k, 1/zeta(2k), the stable factor e^{-b_1(2k-1)} and the exact GGC factor c_k e^{b_1(2k-1)}, in Arb;
#   truncated  the truncated integral without quadrature error: v_1 in Arb ball arithmetic, Gauss Legendre panels of width
#              0.25 (graded towards 0) on [0, 200] with 24 nodes, checked against 12 nodes, the constant part b_1/pi
#              integrated in closed form, and the Weyl tail in closed form, checked against mpmath quad;
#   authors    the scheme of ggc_coef.py rerun line by line on its own grid, checked against ggc_coef.json, and split into
#              its pieces [0, 0.025], [0.025, 60], [60, 200] and the tail;
#   extension  rho integrated on to 500, 1000 and 2000 with the Weyl tail beyond: the error left by the truncated integral
#              shrinks, which confirms that it is the error of the Weyl tail.
# Relative errors are value/reference - 1. Weyl tail error = truncated/exact - 1; quadrature error = ggc_coef.json/truncated
# - 1; their product (1 + .)(1 + .) is 1 + the error of the table. Output: ggc_offset.json (reads ggc_coef.json).
import json, math, os, time
import numpy as np, mpmath as mp
from multiprocessing import Pool
from mpmath.calculus.quadrature import GaussLegendre
from flint import arb, acb, acb_series, ctx

ctx.prec = 128
mp.mp.dps = 30
KS = list(range(1, 7))
WS = [(2 * k - 1) ** 2 for k in KS]
NODES = {d: GaussLegendre(mp.mp).calc_nodes(d, mp.mp.prec) for d in (3, 4)}   # 12 and 24 nodes on [-1, 1]


def v1(th):
    """Re (xi'/xi)(1 + i th) for th > 0, as an Arb ball (as in v1min.py)."""
    s = acb(1, th)
    c = acb_series([s, 1], prec=2).zeta().coeffs()
    return 1 / (1 + th * th) - arb.pi().log() / 2 + acb(arb(1) / 2, th / 2).digamma().real / 2 + (c[1] / c[0]).real


def to_mp(x):
    return mp.mpf(x.mid().str(40, radius=False))


B1 = 1 + arb.const_euler() / 2 - (4 * arb.pi()).log() / 2          # v_1(0) = b_1
B1M = to_mp(B1)


def panel(task):
    """(1/pi) int_a^b log(1 + w/theta^2) (v_1 - b_1 if subtract else v_1) dtheta for every w, by Gauss Legendre."""
    a, b, degree, subtract = task
    mp.mp.dps = 30
    a, b = mp.mpf(a), mp.mpf(b)
    h, m = (b - a) / 2, (a + b) / 2
    acc, rad = [mp.mpf(0)] * len(WS), 0.0
    for x, wt in NODES[degree]:
        t = m + h * x
        v = v1(arb(mp.nstr(t, 45)))
        rad = max(rad, float(v.rad()))
        vm = to_mp(v) - (B1M if subtract else 0)
        for i, w in enumerate(WS):
            acc[i] += wt * mp.log(1 + w / t ** 2) * vm
    return [mp.nstr(h * s / mp.pi, 32) for s in acc], rad


def authors_r(t):
    """The value r of ggc_coef.py at a grid point t > 0, with its own settings (dps 25, mp.diff of log xi)."""
    mp.mp.dps = 25
    xi = lambda s: mp.mpf(0.5) if s == 1 else s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)
    dl = lambda s: mp.diff(lambda z: mp.log(xi(z)), s)
    return float(mp.re(dl(mp.mpc(1, t))) / mp.pi)


def arb_r(t):
    """v_1(t)/pi from Arb at a float grid point, for the authors' grid."""
    return float((v1(arb(float(t))) / arb.pi()).mid())


def F(L, w):
    """int_0^L log(1 + w/theta^2) dtheta."""
    L, w = mp.mpf(L), mp.mpf(w)
    return L * mp.log(1 + w / L ** 2) + 2 * mp.sqrt(w) * mp.atan(L / mp.sqrt(w))


def weyl_tail(L, w, terms=80):
    """int_L^inf log(1 + w/theta^2) (1/2pi) log(theta/2pi) dtheta in closed form (series in w/L^2, valid for w < L^2)."""
    L, w = mp.mpf(L), mp.mpf(w)
    lg, s = mp.log(L / (2 * mp.pi)), mp.mpf(0)
    for m in range(1, terms + 1):
        s += (-1) ** (m + 1) * w ** m / m * L ** (1 - 2 * m) * (lg / (2 * m - 1) + mp.mpf(1) / (2 * m - 1) ** 2)
    return s / (2 * mp.pi)


def arg_zeta_1(T, steps=1800):
    """arg zeta(1 + iT) by continuous variation along the horizontal line from 10 + iT, where zeta is within 1e-3 of 1."""
    vals = [float(acb(arb(10 - 9 * j / steps), T).zeta().arg().mid()) for j in range(steps + 1)]
    return float(np.unwrap(vals)[-1])


def panels(bounds, width, degree, subtract):
    """Panels of the given width between consecutive bounds (each bound a multiple of width), as tasks."""
    out = []
    for a, b in zip(bounds[:-1], bounds[1:]):
        n = int(round((b - a) / width))
        out += [(str(a + j * width), str(a + (j + 1) * width), degree, subtract) for j in range(n)]
    return out


def run(pool, tasks):
    res = pool.map(panel, tasks, chunksize=4)
    tot = [mp.fsum(mp.mpf(r[0][i]) for r in res) for i in range(len(WS))]
    return tot, max(r[1] for r in res)


if __name__ == '__main__':
    t_start = time.time()
    ncpu = int(os.environ.get('SLURM_CPUS_PER_TASK', '8'))
    pool = Pool(ncpu)
    J = json.load(open('ggc_coef.json'))
    out = {'b1': float(B1.mid()), 'cpus': ncpu}

    # exact values in Arb
    ex = []
    for k in KS:
        z = acb(2 * k).zeta().real
        kap = 2 * k * (2 * k - 1) * arb(math.factorial(k - 1)) / arb.pi() ** k
        c = 1 / (kap * z)
        ex.append(dict(k=k, w=(2 * k - 1) ** 2, inv_zeta=to_mp(1 / z), kappa=to_mp(kap), c=to_mp(c),
                       stable=to_mp((-B1 * (2 * k - 1)).exp()), ggc_exact=to_mp(c * (B1 * (2 * k - 1)).exp()),
                       I_exact=to_mp(-c.log())))

    # the truncated integral on [0, 200] without quadrature error: graded panels near 0, then width 0.25;
    # pieces [0, 0.025], [0.025, 60], [60, 200] match the grid of ggc_coef.py
    graded0 = ['0', '1e-6', '1e-5', '1e-4', '1e-3', '5e-3', '0.025']
    graded1 = ['0.025', '0.1', '0.25']
    pieces = {}
    for deg in (4, 3):
        t0 = time.time()
        p0 = [(a, b, deg, True) for a, b in zip(graded0[:-1], graded0[1:])]
        p1 = [(a, b, deg, True) for a, b in zip(graded1[:-1], graded1[1:])] + panels([0.25, 60.0], 0.25, deg, True)
        p2 = panels([60.0, 200.0], 0.25, deg, True)
        S0, r0 = run(pool, p0)
        S1, r1 = run(pool, p1)
        S2, r2 = run(pool, p2)
        pieces[deg] = (S0, S1, S2, max(r0, r1, r2))
        print(f'[0, 200] with {3 * 2 ** (deg - 1)} nodes per panel: {len(p0) + len(p1) + len(p2)} panels, '
              f'{time.time() - t0:.1f} s, largest Arb radius of v_1 {max(r0, r1, r2):.1e}', flush=True)
    S0, S1, S2, rad = pieces[4]
    gl_diff = max(abs(pieces[4][j][i] - pieces[3][j][i]) for j in range(3) for i in range(len(WS)))
    out['check_gauss_legendre_24_vs_12_max_abs_diff'] = float(gl_diff)
    out['check_arb_radius_v1_max'] = rad

    # extension: rho on [200, 2000] (12 nodes per panel of width 0.25), checked with 24 nodes on [1900, 2000]
    t0 = time.time()
    cuts = [200.0, 500.0, 1000.0, 2000.0]
    ext = {}
    acc = [mp.mpf(0)] * len(WS)
    for a, b in zip(cuts[:-1], cuts[1:]):
        E, _ = run(pool, panels([a, b], 0.25, 3, False))
        acc = [acc[i] + E[i] for i in range(len(WS))]
        ext[int(b)] = list(acc)
    Ea, _ = run(pool, panels([1900.0, 2000.0], 0.25, 3, False))
    Eb, _ = run(pool, panels([1900.0, 2000.0], 0.25, 4, False))
    out['check_extension_1900_2000_12_vs_24_max_abs_diff'] = float(max(abs(Ea[i] - Eb[i]) for i in range(len(WS))))
    print(f'extension to 2000: {time.time() - t0:.1f} s', flush=True)

    # the scheme of ggc_coef.py, rerun on its own grid
    t0 = time.time()
    mp.mp.dps = 25
    b1f = float(1 + mp.euler / 2 - mp.log(4 * mp.pi) / 2)
    th = np.concatenate([np.linspace(0, 60, 2401), np.linspace(60.05, 200, 1400)])
    r = np.array([b1f / np.pi] + pool.map(authors_r, list(th[1:]), chunksize=20))
    ra = np.array([b1f / np.pi] + pool.map(arb_r, list(th[1:]), chunksize=20))
    print(f'authors grid: {time.time() - t0:.1f} s; max |r(mp.diff) - r(Arb)| = {np.max(np.abs(r - ra)):.2e}', flush=True)
    out['authors_grid'] = dict(nodes=len(th), first_steps=[float(th[1] - th[0]), float(th[2400] - th[2399]),
                                                          float(th[2401] - th[2400]), float(th[2402] - th[2401])],
                               max_abs_diff_r_mpdiff_vs_arb=float(np.max(np.abs(r - ra))))
    weyl = lambda t: mp.log(t / (2 * mp.pi)) / (2 * mp.pi)
    h0 = mp.mpf(0.025)
    origin_pred = float(2 * h0 * (1 - mp.log(2 * mp.pi) / 2) * B1M / mp.pi)   # trapezoid error at the log singularity
    out['origin_prediction_in_I'] = origin_pred

    rows = []
    for i, k in enumerate(KS):
        w = WS[i]
        E = ex[i]
        # authors' scheme, as in ggc_coef.py
        mp.mp.dps = 25
        f = np.log1p(w / th[1:] ** 2) * r[1:]
        T = np.trapezoid(f, th[1:])
        Ta = np.trapezoid(f[:2400], th[1:2401])
        Tb = np.trapezoid(f[2399:], th[2400:])
        P0 = r[0] * float(mp.quad(lambda x: mp.log(1 + w / x ** 2), [0, th[1]]))
        Wa = float(mp.quad(lambda x: mp.log(1 + w / x ** 2) * weyl(x), [200, mp.inf]))
        I_rerun = T + P0 + Wa
        Ck = 2 * k * (2 * k - 1) * float(mp.factorial(k - 1)) / float(mp.pi) ** k
        stable_f = np.exp(-b1f * (2 * k - 1))
        thor_rerun, rest_rerun = Ck * np.exp(-I_rerun), np.exp(-I_rerun) / stable_f
        fa = np.log1p(w / th[1:] ** 2) * ra[1:]
        T_arb = np.trapezoid(fa, th[1:])
        # the truncated integral without quadrature error
        mp.mp.dps = 30
        c = mp.mpf(B1M) / mp.pi
        A0 = S0[i] + c * F('0.025', w)
        A1 = S1[i] + c * (F(60, w) - F('0.025', w))
        A2 = S2[i] + c * (F(200, w) - F(60, w))
        Wt = weyl_tail(200, w)
        Wq = mp.quad(lambda x: mp.log(1 + w / x ** 2) * mp.log(x / (2 * mp.pi)) / (2 * mp.pi), [200, 400, 1000, mp.inf])
        I_tr = A0 + A1 + A2 + Wt
        thor_tr, ggc_tr = E['kappa'] * mp.exp(-I_tr), mp.exp(-I_tr) / E['stable']
        # the table (ggc_coef.json)
        Jr = J['coef'][i]
        I_json = mp.log(mp.mpf(Jr['C']) / mp.mpf(Jr['thorin']))
        err_table = mp.mpf(Jr['thorin']) / E['inv_zeta'] - 1
        err_table_ggc = mp.mpf(Jr['rest']) / E['ggc_exact'] - 1
        err_weyl = mp.exp(E['I_exact'] - I_tr) - 1
        err_quad = mp.exp(I_tr - I_json) - 1
        dI = {'origin_0_0.025': float(P0 - A0), 'trapezoid_0.025_60': float(Ta - A1),
              'trapezoid_60_200': float(Tb - A2), 'weyl_tail_quad': float(Wa - Wt)}
        extension = {}
        for L in (500, 1000, 2000):
            I_L = A0 + A1 + A2 + ext[L][i] + weyl_tail(L, w)
            e = mp.exp(E['I_exact'] - I_L) - 1
            extension[str(L)] = dict(rel_err=float(e), rel_err_over_w=float(e / w))
        row = dict(k=k, w=w, inv_zeta_2k=float(E['inv_zeta']), c_k=float(E['c']), kappa_k=float(E['kappa']),
                   stable=float(E['stable']), ggc_exact=float(E['ggc_exact']), I_exact=float(E['I_exact']),
                   I_trunc_0_200=float(A0 + A1 + A2), weyl_tail_200=float(Wt), weyl_tail_series_minus_quad=float(Wt - Wq),
                   I_trunc=float(I_tr), thorin_trunc=float(thor_tr), ggc_trunc=float(ggc_tr),
                   thorin_json=Jr['thorin'], rest_json=Jr['rest'], I_json=float(I_json),
                   rel_err_table=float(err_table), rel_err_table_ggc_column=float(err_table_ggc),
                   rel_err_weyl_tail=float(err_weyl), rel_err_weyl_tail_over_w=float(err_weyl / w),
                   rel_err_quadrature=float(err_quad),
                   check_product=float((1 + err_weyl) * (1 + err_quad) - 1 - err_table),
                   quadrature_error_in_I_by_piece=dI, quadrature_error_in_I_total=float(I_json - I_tr),
                   trapezoid_0_025_60_minus_origin_prediction=dI['trapezoid_0.025_60'] - origin_pred,
                   rel_err_quadrature_without_origin=float(mp.exp(I_tr - I_json + origin_pred) - 1),
                   rerun=dict(I=I_rerun, thorin=thor_rerun, rest=rest_rerun,
                              rel_diff_thorin_vs_json=thor_rerun / Jr['thorin'] - 1,
                              rel_diff_rest_vs_json=rest_rerun / Jr['rest'] - 1,
                              I_trapezoid_arb_minus_mpdiff=float(T_arb - T)),
                   extension=extension)
        rows.append(row)
        print(f"k={k}  1/zeta(2k)={row['inv_zeta_2k']:.7f}  table={Jr['thorin']:.7f}  truncated={row['thorin_trunc']:.7f}"
              f"  err table {row['rel_err_table']:+.4e} = weyl {row['rel_err_weyl_tail']:+.4e}"
              f" (/w {row['rel_err_weyl_tail_over_w']:.4e}) + quad {row['rel_err_quadrature']:+.4e}"
              f"  | GGC trunc {row['ggc_trunc']:.7f} exact {row['ggc_exact']:.7f}"
              f"  | rerun vs json {row['rerun']['rel_diff_thorin_vs_json']:+.1e}", flush=True)
        print('      quadrature error in I by piece: ' + ', '.join(f'{a} {b:+.4e}' for a, b in dI.items()) +
              f"; extension rel err/w at 500, 1000, 2000: " +
              ', '.join(f"{extension[L]['rel_err_over_w']:+.3e}" for L in ('500', '1000', '2000')), flush=True)
    out['rows'] = rows

    # first order in w, the error of the Weyl tail from a cutoff L is w int_L^inf (rho - Weyl) theta^{-2} dtheta, and
    # rho - Weyl = (1/pi) d/dtheta arg zeta(1 + i theta) + (11/12)/(pi theta^2) + O(theta^{-4}); integrating by parts,
    # it is -w arg zeta(1 + iL)/(pi L^2) + (11/12) w/(3 pi L^3) + (2w/pi) int_L^inf arg zeta(1 + i theta) theta^{-3} dtheta
    first = {}
    for L in (200, 500, 1000, 2000):
        a = arg_zeta_1(L)
        pred = -a / (math.pi * L ** 2) + (11 / 12) / (3 * math.pi * L ** 3)
        obs = rows[0]['rel_err_weyl_tail_over_w'] if L == 200 else rows[0]['extension'][str(L)]['rel_err_over_w']
        first[str(L)] = dict(arg_zeta_1_plus_iL=a, first_order_prediction_over_w=pred, observed_over_w_k1=obs,
                             observed_minus_prediction=obs - pred)
        print(f'cutoff {L}: arg zeta(1 + i{L}) = {a:+.6f}, predicted rel err/w {pred:+.4e}, observed (k=1) {obs:+.4e}',
              flush=True)
    out['weyl_tail_first_order'] = first
    out['seconds'] = time.time() - t_start
    json.dump(out, open('ggc_offset.json', 'w'), indent=1)
    print(f"origin prediction 2h(1 - log(2pi)/2) b_1/pi = {origin_pred:.5e}; GL 24 vs 12 max diff {float(gl_diff):.1e};"
          f" extension 12 vs 24 on [1900, 2000] {out['check_extension_1900_2000_12_vs_24_max_abs_diff']:.1e};"
          f" total {out['seconds']:.0f} s")
