# heat_ladder.py: the ladder functions along the de Bruijn Newman flow, Table pk:tab:heatladder in Section pk:sec:ladder
# (ch/pakes2.tex). Added for the book in October 2026, for the merge of the first author's chapter on Pakes, the Lambert
# law and the Thorin condition (misc/reviews/nick-pakes/, triage-C.md; the chapter's table is recomputed here on a window
# that does not depend on t). With H_t(z) = int_0^oo e^{t u^2} Phi_dBN(u) cos(z u) du (Section sec:ch21:flow), so that
# H_0(z) = Xi(z/2)/8, and for t >= Lambda the Thorin atoms a_k(t) = (z_k(t)/2)^2 of R_t(u) = H_t(0)/H_t(2i sqrt u)
# (Theorem rz:thm:dbn), the ladder functions of R_t are W^[j](x; t) = sum_k beta_j(a_k(t) x) e^{-a_k(t) x}, with
# beta_0 = 1, beta_{j+1}(u) = u beta_j(u) - u beta_j'(u). For t = 0, 0.2, 1, 4, 10, 20, 30, 40, 50 the script
#  1. computes H_t by Gauss Legendre quadrature on [0, 5/4] with 1536 nodes in 480 bit arithmetic (python-flint, Arb; the
#     integrand is below 1e-160 beyond u = 5/4 for t <= 50), and checks it against 1024 nodes at z = 100, 300, 590 and, at
#     t = 0, against Xi(z/2)/8 computed from the zeta function;
#  2. finds the zeros of H_t on (0, 600] by a sign scan of step 0.2 and 52 bisections, and at t = 0 compares them with
#     2 gamma_k (data/zeros_1700.txt) and their number with N(300) = 138 (Arb);
#  3. evaluates W^[j](x; t), j <= 5, from these atoms in double precision on y = -log x in [0, 8] with step 0.001, and
#     records the atoms a_1 and a_2/a_1, the number of zeros below z = 10 (a scan that starts at z = 10 misses the first
#     zero for t = 40 and t = 50), the ratios W^[2]/W^[0] and W^[3]/W^[0] at y = 6, 7, 8, whether they decrease on the grid and where
#     they are least, the sets where W^[j] < 0, the largest n with W^[1], ..., W^[n+1] >= 0 on the grid (the class L_n of
#     Urbanik that R_t reaches on the window), and an estimate of the part of W^[j] at y = 8 that the zeros above z = 600
#     would add (from the density (1/4pi) log(z/4pi) of the zeros of H_0, relative to W^[j]);
#  4. locates where the class drops from L_2 to L_1 on the window: W^[3] >= 0 on the grid at t = 20 and not at t = 30, and
#     the least value of W^[3] on the grid is computed for t = 21, ..., 29 and then for the tenths of the first unit
#     interval on which it changes sign (pools of nine processes).
# For y <= 0 every term is positive, since a_1(t) x >= a_1(t) exceeds the largest root of beta_j for j <= 5 (computed and
# checked below). The computation is numerical: the quadrature error is checked, not bounded, and the truncation at
# z = 600 is estimated, not bounded.
# Run on Hopper (about three minutes on 9 cores), from a folder containing zeros_1700.txt: python3 heat_ladder.py; writes
# heat_ladder.json there.
import json
import math
import time
from multiprocessing import Pool

from flint import arb, acb, ctx

PREC = 480
TS = ['0', '0.2', '1', '4', '10', '20', '30', '40', '50']
U = arb(5) / 4
ZHI, ZSTEP = 600.0, 0.2
YS = [round(0.001 * i, 3) for i in range(8001)]
BETA = [[1], [0, 1], [0, -1, 1], [0, 1, -3, 1], [0, -1, 7, -6, 1], [0, 1, -15, 25, -10, 1]]


def nodes(N):
    ctx.prec = PREC
    xs, ws = [], []
    for k in range(N):
        x, w = arb.legendre_p_root(N, k, weight=True)
        xs.append(x)
        ws.append(w)
    return xs, ws


def Phi(u):  # Phi_dBN(u) = sum_n (2 pi^2 n^4 e^{9u} - 3 pi n^2 e^{5u}) exp(-pi n^2 e^{4u}); n <= 12 suffices
    p = arb.pi()
    e4, e5, e9 = (4 * u).exp(), (5 * u).exp(), (9 * u).exp()
    return sum((2 * p ** 2 * n ** 4 * e9 - 3 * p * n ** 2 * e5) * (-p * n * n * e4).exp() for n in range(1, 13))


def make_H(t, N):
    ctx.prec = PREC
    xs, ws = nodes(N)
    us = [U * (x + 1) / 2 for x in xs]
    fs = [U * w / 2 * (arb(t) * u * u).exp() * Phi(u) for u, w in zip(us, ws)]
    return lambda z: sum(f * (arb(z) * u).cos() for f, u in zip(fs, us))


def sgn(v):
    return 1 if v > 0 else (-1 if v < 0 else 0)


def zeros_of(t):
    ctx.prec = PREC
    H = make_H(t, 1536)
    grid = [ZSTEP * i for i in range(int(round(ZHI / ZSTEP)) + 1)]
    sg = [sgn(H(z)) for z in grid]
    zs = []
    for i in range(len(grid) - 1):
        if sg[i] * sg[i + 1] < 0:
            lo, hi, slo = arb(repr(grid[i])), arb(repr(grid[i + 1])), sg[i]
            for _ in range(52):
                mid = (lo + hi) / 2
                sm = sgn(H(mid))
                if sm == 0:
                    break
                if sm == slo:
                    lo = mid
                else:
                    hi = mid
            zs.append(float(((lo + hi) / 2).mid()))
    H2 = make_H(t, 1024)
    quad = {str(z): float(abs((H(z) - H2(z)) / H(z)).mid()) for z in (100, 300, 590)}
    res = dict(n_zeros=len(zs), unresolved_signs=sum(1 for v in sg if v == 0), H0_positive=sg[0] > 0,
               quad_rel_diff_1024_vs_1536=quad, zeros=zs)
    if t == '0':
        xi = {}
        for z in (100, 300, 590):
            sv = acb(arb(1) / 2, arb(z) / 2)
            x8 = (sv * (sv - 1) / 2 * (-(sv / 2) * acb.pi().log()).exp() * (sv / 2).gamma() * sv.zeta()).real / 8
            xi[str(z)] = dict(H=H(z).str(15), Xi_over_8=x8.str(15), rel_diff=float(abs((H(z) - x8) / x8).mid()))
        res['xi_check'] = xi
        res['N_300'] = float(arb(300).zeta_nzeros().mid())
    return res


def ladder(atoms):
    W = [[0.0] * len(YS) for _ in range(6)]
    for i, y in enumerate(YS):
        x = math.exp(-y)
        for a in atoms:
            u = a * x
            if u > 745:
                break
            e = math.exp(-u)
            for j in range(6):
                p = 0.0
                for c in reversed(BETA[j]):
                    p = p * u + c
                W[j][i] += p * e
    return W


def neg_sets(v):
    res, i = [], 0
    while i < len(v):
        if v[i] < 0:
            k = i
            while k + 1 < len(v) and v[k + 1] < 0:
                k += 1
            res.append((YS[i], YS[k + 1] if k + 1 < len(v) else YS[k]))
            i = k + 1
        else:
            i += 1
    return res


def tail_estimate(j, y, W):  # the zeros above z = 600, with density (1/4pi) log(z/4pi), at y; relative to |W^[j](y)|
    x = math.exp(-y)
    tot, z, h = 0.0, 600.0, 0.05
    while True:
        u = z * z / 4 * x
        p = sum(abs(c) * u ** k for k, c in enumerate(BETA[j]))
        term = p * math.exp(-u) * math.log(z / (4 * math.pi)) / (4 * math.pi) * h
        tot += term
        if u > 200:
            break
        z += h
    return tot / abs(W)


def worker(t):
    t1 = time.time()
    r = zeros_of(t)
    atoms = [z * z / 4 for z in r['zeros']]
    W = ladder(atoms)
    row = {k: r[k] for k in r if k != 'zeros'}
    row['zeros_1_to_5'] = r['zeros'][:5]
    row['n_zeros_below_10'] = sum(1 for z in r['zeros'] if z < 10)
    row['a1'] = atoms[0]
    row['a2_over_a1'] = atoms[1] / atoms[0]
    if t == '0':
        gam = [float(l) for l in open('zeros_1700.txt') if l.strip()]
        row['max_abs_z_over_2_minus_gamma'] = max(abs(z / 2 - g) for z, g in zip(r['zeros'], gam))
    row['ratios_at_y_6_7_8'] = {str(y): [W[2][k] / W[0][k], W[3][k] / W[0][k]] for y, k in ((6, 6000), (7, 7000), (8, 8000))}
    row['ratios_decrease_on_0_8'] = [all(W[j][k + 1] / W[0][k + 1] < W[j][k] / W[0][k] for k in range(8000)) for j in (2, 3)]
    row['negative_sets'] = {j: neg_sets(W[j]) for j in range(6)}
    row['ratio_least_on_0_8'] = {j: min((W[j][k] / W[0][k], YS[k]) for k in range(8001)) for j in (2, 3)}
    row['urbanik_index_on_window'] = max([n for n in range(5) if all(neg_sets(W[j]) == [] for j in range(1, n + 2))],
                                         default=-1)
    row['W_at_y8'] = [W[j][8000] for j in range(6)]
    row['tail_estimate_rel_at_y8'] = [tail_estimate(j, 8.0, W[j][8000]) for j in range(6)]
    row['seconds'] = round(time.time() - t1, 1)
    return t, row


def w3_min(t):  # the least value of W^[3](x; t) on the grid y in [0, 8], where it is attained, and the set where W^[3] < 0
    r = zeros_of(t)
    W = ladder([z * z / 4 for z in r['zeros']])
    m = min((W[3][k], YS[k]) for k in range(len(YS)))
    return t, dict(w3_min=m[0], at_y=m[1], w3_negative=neg_sets(W[3]), z1=r['zeros'][0], n_zeros=r['n_zeros'],
                   a2_over_a1=(r['zeros'][1] / r['zeros'][0]) ** 2, unresolved_signs=r['unresolved_signs'])


def largest_root(j):  # the largest positive root of beta_j, by a scan of step 0.001 on (0, 20] and 60 bisections
    f = lambda u: sum(c * u ** k for k, c in enumerate(BETA[j]))
    last = None
    for i in range(1, 20001):
        a, b = 0.001 * (i - 1) + 1e-9, 0.001 * i
        if f(b) == 0:  # a root on the grid, as u = 1 for beta_2
            last = b
        elif f(a) * f(b) < 0:
            for _ in range(60):
                m = (a + b) / 2
                if f(a) * f(m) <= 0:
                    b = m
                else:
                    a = m
            last = (a + b) / 2
    assert last is not None and all(f(last + 0.01 * i) > 0 for i in range(1, 2001))  # positive beyond it, up to 20 above
    return last


def first_negative(keys, res):  # the first key, in order, with a negative least value of W^[3]
    return next((k for k in keys if res[k]['w3_min'] < 0), None)


if __name__ == '__main__':
    T0 = time.time()
    with Pool(len(TS)) as pool:
        res = dict(pool.map(worker, TS))
    out = dict(prec_bits=PREC, nodes=1536, U='5/4', z_range=[0, ZHI], y_window=[0, 8], per_t=res,
               urbanik_index={t: res[t]['urbanik_index_on_window'] for t in TS})
    # for y <= 0 every term of W^[j], j <= 5, is positive when a_1(t) exceeds the largest root of beta_j
    out['beta_largest_root'] = {j: largest_root(j) for j in range(2, 6)}
    out['a1_exceeds_largest_roots'] = all(res[t]['a1'] > max(out['beta_largest_root'].values()) for t in TS)
    # where the index drops from 2 to 1 on the window, between t = 20 and t = 30
    assert res['20']['urbanik_index_on_window'] == 2 and res['30']['urbanik_index_on_window'] == 1
    K1 = [str(k) for k in range(21, 30)]
    with Pool(len(K1)) as pool:
        r1 = dict(pool.map(w3_min, K1))
    k1 = first_negative(K1 + ['30'], {**r1, '30': dict(w3_min=-1)})
    lo = int(k1) - 1
    K2 = ['%d.%d' % (lo, i) for i in range(1, 10)]
    with Pool(len(K2)) as pool:
        r2 = dict(pool.map(w3_min, K2))
    k2 = first_negative(K2 + [str(lo + 1)], {**r2, str(lo + 1): dict(w3_min=-1)})
    i2 = (K2 + [str(lo + 1)]).index(k2)
    below = str(lo) if i2 == 0 else K2[i2 - 1]
    out['index_drop'] = dict(w3_min_integers=r1, w3_min_tenths=r2, last_nonnegative=below, first_negative=k2,
                             monotone=all((r1[a]['w3_min'] < 0) <= (r1[b]['w3_min'] < 0) for a, b in zip(K1, K1[1:]))
                             and all((r2[a]['w3_min'] < 0) <= (r2[b]['w3_min'] < 0) for a, b in zip(K2, K2[1:])))
    out['seconds'] = round(time.time() - T0, 1)
    json.dump(out, open('heat_ladder.json', 'w'), indent=1)
    for t in TS:
        r = res[t]
        print(t, r['n_zeros'], 'a1 = %.8f, a2/a1 = %.7f' % (r['a1'], r['a2_over_a1']),
              {y: [round(v, 4) for v in r['ratios_at_y_6_7_8'][y]] for y in r['ratios_at_y_6_7_8']},
              r['ratios_decrease_on_0_8'], 'W4 < 0:', r['negative_sets'][4], 'W5 < 0:', r['negative_sets'][5],
              'quad:', r['quad_rel_diff_1024_vs_1536'], 'tail:', ['%.1e' % v for v in r['tail_estimate_rel_at_y8']], flush=True)
    print('t = 0:', res['0'].get('max_abs_z_over_2_minus_gamma'), res['0'].get('N_300'), json.dumps(res['0'].get('xi_check')))
    print('index (L_n reached on the window):', out['urbanik_index'])
    print('largest roots of beta_j:', out['beta_largest_root'], 'a_1(t) above them:', out['a1_exceeds_largest_roots'])
    d = out['index_drop']
    print('least W3 on the window:', {k: round(v['w3_min'], 5) for k, v in {**d['w3_min_integers'], **d['w3_min_tenths']}.items()})
    print('W3 >= 0 on the grid at t =', d['last_nonnegative'], 'and not at t =', d['first_negative'], 'monotone:', d['monotone'],
          'seconds', out['seconds'])
