# wigbiane_flow.py: Laguerre's quantity along the backward de Bruijn Newman flow up to each collision, for the section on Biane's
# theorem and the collisions of zeros (ch/ch21.tex, after Section sec:ch21:hankel; draft
# misc/reviews/nick-wigner/merge/wigbiane.tex): the critical point and Biane bounds along the flow, the corrections of the local
# bound, and two sided brackets for the collision times. Added for the book in October 2026, from the first author's note on
# the semicircle law and the Thorin framework (misc/reviews/nick-wigner/, version 7, Section 8), recomputed from H_t itself
# (the note follows the zero dynamics of 400 to 1200 zeros, a truncation).
# H_t(z) = int_0^oo e^{t u^2} Phi_dBN(u) cos(z u) du (Section sec:ch21:flow) by the trapezoid rule with step h = 1/360 on
# [0, 1.35] in Arb at 800 bits. The sum equals H_t(z) plus the aliasing terms (H_t(2 pi k/h - z) + H_t(2 pi k/h + z))/2,
# k != 0, below 10^{-100} |H_t(z)| for |Re z| <= 830, and the integrand beyond 1.35 is below 10^{-290}. With
# C_n(w) = sum_j w_j u_j^n w^j, w_j = h e^{t u_j^2} Phi_dBN(u_j) (halved at j = 0), the n-th z derivative is
# (i^n C_n(e^{izh}) + (-i)^n C_n(e^{-izh}))/2, that is Re(i^n C_n(e^{izh})) for real z, and d/dt H_t = -H_t'' holds exactly.
# The collisions are read from dbn_collisions.json (written by code/dbn_collisions.py). For each pair (k, k + 1) listed there
# the double zero H_t = H_t' = 0 is refined by Newton's method in (z, t). On a grid of about 140 times s = -t in [0, s_c)
# (uniform to 0.95 s_c, geometric towards s_c, and the multiples of 0.1) the script follows the two zeros by Newton's method,
# the minimiser u of Lambda = -(log H_t)'' over their gap (12 samples, then safeguarded Newton on Lambda' = 0) and the critical
# point c (H_t'(c) = 0). From the logarithmic derivatives at u it records m = Lambda(u) = S_2, S_4 = Lambda''(u)/6,
# eta_dir = 6 S_4/m^2 - 3 (so that dm/ds = (1 + eta_dir) m^2), P = (u - z_k)^{-2} + (u - z_{k+1})^{-2}, R = m - P,
# Q = 6 S_4 - 3 P^2, and L = Lambda(c). It also follows every pair a +- ib that has left the axis by time t (the collisions of
# dbn_collisions.json, by complex Newton's method from z_c + i (2(t_c - t))^{1/2} and then from the previous position moved by
# the velocity -H''/H'; the mirror pairs -a +- ib included) and adds, for the zeros 2 gamma_j, 207 <= j <= 2000, beyond the
# census, a model tail in which every such zero has left the axis with b^2 = 4|t|. From these it computes B_p =
# sum b^2/|u - a|^p, W = sum |u - z|^{-4} over the nonreal zeros, and the correction eta of the local bound. Per pair it
# reports m(0), L(0), R(0), rho = sup R, the largest eta_dir, eta and W, the largest -Q, the bracket
# 1/((1 + eta_bar) m(0)) <= s_c <= s*, s* = (r_1 - r_2)^{-1} log((m(0) - r_2)/(m(0) - r_1)), r_{1,2} = 3 rho +-
# (9 rho^2 + 6 W_bar)^{1/2}, in the zero free form (eta_dir, and 6 W_bar replaced by the largest -Q) and in the form with the
# nonreal zeros (eta, W), the largest ratios of L(s) and m(s) to L(0)/(1 - s L(0)) and m(0)/(1 - s m(0)), the largest
# b^2/(2(t_c - t)) over the tracked pairs, central difference checks of dm/ds and dL/ds at s = 0, and L(s), m(s) at
# s = 0.1, 0.2, ... Reference ordinates from Arb (acb.zeta_zeros). Run from the folder holding dbn_collisions.json:
# python3 wigbiane_flow.py; writes wigbiane_flow.json there. Uses SLURM_CPUS_PER_TASK processes (default 1).
import json
import math
import multiprocessing as mpr
import os
import time

from flint import acb, arb, arb_poly, ctx

t0 = time.time()
out = {}
PREC = 800
ctx.prec = PREC
HU = arb(1) / 360
NU = 486
I = acb(0, 1)
IP = [acb(1), acb(0, 1), acb(-1), acb(0, -1)]


def rec(k, v):
    out[k] = v
    print(k, str(v)[:2000], round(time.time() - t0, 1), flush=True)
    json.dump(out, open('wigbiane_flow.json', 'w'), indent=1)


def S(x, n=12):
    return x.mid().str(n, radius=False)


def F(x):
    return float(x.mid())


def phi_dbn(u):
    # Phi_dBN(u) = sum_n (2 pi^2 n^4 e^{9u} - 3 pi n^2 e^{5u}) exp(-pi n^2 e^{4u})
    pi = arb.pi()
    e4, e5, e9 = (4 * u).exp(), (5 * u).exp(), (9 * u).exp()
    tot = arb(0)
    eps = arb(2) ** (-(PREC + 40))
    n = 1
    while True:
        term = (2 * pi * pi * n ** 4 * e9 - 3 * pi * n * n * e5) * (-pi * n * n * e4).exp()
        tot += term
        if n >= 2 and abs(term) < eps * abs(tot):
            break
        n += 1
    return tot.mid()


U = [HU * j for j in range(NU + 1)]
PHI = [phi_dbn(u) for u in U]
U2 = [u * u for u in U]


def polys(t, nmax=4):
    t = arb(t)
    w = [HU * (t * U2[j]).exp() * PHI[j] for j in range(NU + 1)]
    w[0] = w[0] / 2
    res = []
    cur = w
    for n in range(nmax + 1):
        res.append(arb_poly(cur))
        if n < nmax:
            cur = [cur[j] * U[j] for j in range(NU + 1)]
    return res


def ders_real(z, C, nmax=4):
    om = (I * (z * HU)).exp()
    vals = []
    for n in range(nmax + 1):
        c = C[n](om)
        vals.append((c.real, -c.imag, -c.real, c.imag)[n % 4])
    return vals


def ders_cplx(z, C, nmax=2):
    om = (I * z * HU).exp()
    omi = 1 / om
    return [(IP[n % 4] * C[n](om) + IP[(-n) % 4] * C[n](omi)) / 2 for n in range(nmax + 1)]


def logders(z, C):
    # Lambda, Lambda', Lambda'' and q = H'/H at the real point z
    H = ders_real(z, C, 4)
    a1, a2, a3, a4 = H[1] / H[0], H[2] / H[0], H[3] / H[0], H[4] / H[0]
    f2 = a2 - a1 * a1
    f3 = a3 - 3 * a1 * a2 + 2 * a1 ** 3
    f4 = a4 - 4 * a1 * a3 - 3 * a2 * a2 + 12 * a1 * a1 * a2 - 6 * a1 ** 4
    return -f2, -f3, -f4, a1


def newton_real(z, C):
    z = arb(z).mid()
    for _ in range(80):
        H = ders_real(z, C, 1)
        dz = (H[0] / H[1]).mid()
        z = (z - dz).mid()
        if abs(F(dz)) < 1e-60 * (1 + abs(F(z))):
            return z, True
    return z, False


def cnewton(z, C):
    for _ in range(80):
        H = ders_cplx(z, C, 1)
        dz = H[0] / H[1]
        z = acb((z.real - dz.real).mid(), (z.imag - dz.imag).mid())
        if abs(complex(F(dz.real), F(dz.imag))) < 1e-50 * (1 + abs(F(z.real))):
            H = ders_cplx(z, C, 2)
            v = -H[2] / H[1]
            return z, acb(v.real.mid(), v.imag.mid()), True
    return z, None, False


def double_zero(z, t):
    z, t = arb(z), arb(t)
    for _ in range(40):
        C = polys(t, 3)
        H = ders_real(z, C, 3)
        a, b, c, d = H[1], -H[2], H[2], -H[3]
        det = a * d - b * c
        dz = ((d * H[0] - b * H[1]) / det).mid()
        dt = ((-c * H[0] + a * H[1]) / det).mid()
        z, t = (z - dz).mid(), (t - dt).mid()
        if abs(F(dz)) < 1e-60 and abs(F(dt)) < 1e-60:
            return z, t, True
    return z, t, False


def minimiser(lo, hi, C, u0=None):
    xs = [(lo + (hi - lo) * arb(j) / 13).mid() for j in range(1, 13)]
    if u0 is not None and F(lo) < F(u0) < F(hi):
        xs.append(u0)
    xs.sort(key=F)
    vals = [F(logders(x, C)[0]) for x in xs]
    i = min(range(len(xs)), key=lambda j: vals[j])
    a = xs[i - 1] if i > 0 else lo
    b = xs[i + 1] if i + 1 < len(xs) else hi
    u = xs[i]
    for _ in range(100):
        L0, L1, L2, q = logders(u, C)
        if F(L1) > 0:
            b = u
        else:
            a = u
        un = (u - L1 / L2).mid() if F(L2) > 0 else ((a + b) / 2).mid()
        if not (F(a) < F(un) < F(b)):
            un = ((a + b) / 2).mid()
        if abs(F(un - u)) < 1e-45:
            return un
        u = un
    return u


def critical(lo, hi, C, c0=None):
    a, b = lo, hi
    c = c0 if (c0 is not None and F(lo) < F(c0) < F(hi)) else ((lo + hi) / 2).mid()
    for _ in range(150):
        H = ders_real(c, C, 2)
        if F(H[1] / H[0]) > 0:
            a = c
        else:
            b = c
        cn = (c - H[1] / H[2]).mid()
        if not (F(a) < F(cn) < F(b)):
            cn = ((a + b) / 2).mid()
        if abs(F(cn - c)) < 1e-45:
            return cn
        c = cn
    return c


def sstar(m0, rho, Wb):
    root = math.sqrt(9 * rho * rho + 6 * Wb)
    r1, r2 = 3 * rho + root, 3 * rho - root
    if m0 <= r1:
        return None
    if r1 - r2 < 1e-13 * m0:
        return 1 / m0
    return math.log((m0 - r2) / (m0 - r1)) / (r1 - r2)


def tail_sums(uf, s, gam_tail):
    # model tail: each zero +-2 gamma_j, 207 <= j <= 2000, half of a pair with b^2 = 4|t|
    b2 = 4 * s
    B = {4: 0.0, 5: 0.0, 6: 0.0}
    W = 0.0
    if s <= 0:
        return B, W
    for g in gam_tail:
        for a in (2 * g, -2 * g):
            d = abs(uf - a)
            d4 = d ** 4
            B[4] += 0.5 * b2 / d4
            B[5] += 0.5 * b2 / (d4 * d)
            B[6] += 0.5 * b2 / (d4 * d * d)
            W += 1 / (d * d + b2) ** 2
    return B, W


def run_pair(args):
    col, gk, allc, gam_tail = args
    zc, tc, ok = double_zero(arb(col['z']), arb(col['t']))
    res = {'pair': col['pair'], 'refined z_c': S(zc, 16), 'refined t_c': S(tc, 14), 'Newton converged': ok,
           '|t_c - census|': abs(F(tc) - float(col['t'])), '|z_c - census|': abs(F(zc) - float(col['z']))}
    sc = -F(tc)
    xs = [j / 100 for j in range(96)] + [1 - 0.05 * 2 ** (-j / 2) for j in range(1, 37)]
    grid = sorted(set([x * sc for x in xs] + [0.1 * j for j in range(1, 40) if 0.1 * j < 0.999 * sc]))
    za, zb = arb(repr(2 * gk[0])), arb(repr(2 * gk[1]))
    za0 = zb0 = None
    hist = []
    u = c = None
    cplx = {}
    worst_b2 = 0.0
    prev_s = None
    for s in grid:
        C = polys(arb(-s) if s else arb(0), 4)
        pa, pb = za, zb
        if prev_s is not None and sc - s < 0.02 * sc:
            half = math.sqrt(2 * (sc - s))
            if 2 * half < F(zb - za):
                pa, pb = zc - half, zc + half
        za, okA = newton_real(pa, C)
        zb, okB = newton_real(pb, C)
        if not (okA and okB and F(za) < F(zb)):
            res['failed at s'] = s
            break
        if za0 is None:
            za0, zb0 = za, zb
        u = minimiser(za, zb, C, u)
        L0, L1, L2, q = logders(u, C)
        c = critical(za, zb, C, c)
        Lc = logders(c, C)
        m = F(L0)
        Pa = 1 / (u - za) ** 2 + 1 / (u - zb) ** 2
        # R, Q and eta_dir in ball arithmetic: near the collision m and P are large and nearly equal
        st = {'s': s, 'm': m, 'R': F(L0 - Pa), 'P': F(Pa), 'eta_dir': F(L2 / (L0 * L0) - 3), 'Q': F(L2 - 3 * Pa * Pa),
              'L': F(Lc[0]), 'S3 rel': abs(F(L1)) / 2 / m ** 1.5}
        B = {4: 0.0, 5: 0.0, 6: 0.0}
        W = 0.0
        npairs = 0
        uf = F(u)
        for (ids, zcs, tcs) in allc:
            tcf = float(tcs)
            key = tuple(ids)
            if tcf <= -s or key == tuple(col['pair']):
                continue
            if key not in cplx:
                g = acb(arb(zcs), arb(repr(math.sqrt(2 * (tcf + s)))))
            else:
                zo, vo, so = cplx[key]
                g = zo + vo * (s - so)
            zz, v, okc = cnewton(g, C)
            if not okc or F(zz.imag) <= 0:
                res.setdefault('complex Newton failures', []).append([list(ids), s])
                cplx.pop(key, None)
                continue
            cplx[key] = (zz, v, s)
            npairs += 1
            a, b = F(zz.real), F(zz.imag)
            if tcf + s > 1e-4:
                worst_b2 = max(worst_b2, b * b / (2 * (tcf + s)))
            for aa in (a, -a):
                d = abs(uf - aa)
                for p in B:
                    B[p] += b * b / d ** p
                W += 2 / (d * d + b * b) ** 2
        Bt, Wt = tail_sums(uf, s, gam_tail)
        for p in B:
            B[p] += Bt[p]
        W += Wt
        eta = (36 * m * B[4] + 108 * B[4] ** 2 + 72 * math.sqrt(m + 6 * B[4]) * B[5] + 120 * B[6]) / m ** 2
        st.update({'W': W, 'eta': eta, 'pairs': npairs, 'tail share of B4': Bt[4] / B[4] if B[4] > 0 else 0.0})
        hist.append(st)
        prev_s = s
    if not hist:
        return res
    H0 = hist[0]
    m0, L0v = H0['m'], H0['L']
    rho = max(x['R'] for x in hist)
    etad = max(x['eta_dir'] for x in hist)
    etab = max(x['eta'] for x in hist)
    Wb = max(x['W'] for x in hist)
    qb = max(0.0, max(-x['Q'] for x in hist))
    res.update({'s_c': sc, 'last s': hist[-1]['s'], 'm(0)': m0, 'census m(0)': float(col['m0']), 'L(0) = Lambda(c)': L0v,
                'census Lambda(c)': float(col['Lc']), 'R(0)': H0['R'], 'rho = sup R': rho, 'largest eta_dir': etad,
                'largest eta (B form)': etab, 'largest W': Wb, 'largest -Q (0 if negative)': qb,
                'largest |S3|/S2^(3/2)': max(x['S3 rel'] for x in hist),
                'complex pairs (census) at the end': hist[-1]['pairs'],
                'tail share of B4 where eta is largest': max(hist, key=lambda x: x['eta'])['tail share of B4'],
                'tail share of B4 at the last point': hist[-1]['tail share of B4'],
                'largest b^2/(2(t_c - t)) of the tracked pairs': worst_b2})
    res['lower 1/m(0)'] = 1 / m0
    res['lower 1/((1 + max(eta_dir, 0)) m(0))'] = 1 / ((1 + max(etad, 0)) * m0)
    res['lower 1/((1 + eta) m(0))'] = 1 / ((1 + etab) * m0)
    res['upper s* (rho, -Q)'] = sstar(m0, rho, qb / 6)
    res['upper s* (rho, W)'] = sstar(m0, rho, Wb)
    res['frozen R predictor s*(R(0))'] = sstar(m0, H0['R'], 0.0)
    g0 = F(zb0 - za0)
    res['gap g at s = 0 (z units)'] = g0
    res['first order predictor g^2/8 + g^4 R(0)/32'] = g0 ** 2 / 8 + g0 ** 4 * H0['R'] / 32
    res['largest L(s)(1 - s L(0))/L(0), s < 1/L(0)'] = max([x['L'] * (1 - x['s'] * L0v) / L0v for x in hist if x['s'] * L0v < 1])
    res['largest m(s)(1 - s m(0))/m(0), s < 1/m(0)'] = max([x['m'] * (1 - x['s'] * m0) / m0 for x in hist if x['s'] * m0 < 1])
    res['s, L(s), L(0)/(1 - s L(0)), m(s), m(0)/(1 - s m(0)) at s = 0.1 j'] = [
        [round(x['s'], 6), x['L'], L0v / (1 - x['s'] * L0v) if x['s'] * L0v < 1 else None, x['m'],
         m0 / (1 - x['s'] * m0) if x['s'] * m0 < 1 else None]
        for x in hist if x['s'] > 0 and abs(x['s'] * 10 - round(x['s'] * 10)) < 1e-9]
    res['1/L(0), 1/m(0)'] = [1 / L0v, 1 / m0]
    # central differences at s = 0
    hstep = 1e-8
    C0 = polys(arb(0), 4)
    z0a, _ = newton_real(arb(repr(2 * gk[0])), C0)
    z0b, _ = newton_real(arb(repr(2 * gk[1])), C0)
    u0 = minimiser(z0a, z0b, C0)
    c0 = critical(z0a, z0b, C0)
    Lu, Lc0 = logders(u0, C0), logders(c0, C0)
    vals = {}
    for sg in (1, -1):
        Cd = polys(arb(-sg * hstep), 4)
        da, _ = newton_real(z0a, Cd)
        db, _ = newton_real(z0b, Cd)
        vals[sg] = (logders(minimiser(da, db, Cd, u0), Cd)[0], logders(critical(da, db, Cd, c0), Cd)[0])
    res['dm/ds at 0: difference, Lambda\'\' - 2 Lambda^2 at u'] = [F((vals[1][0] - vals[-1][0]) / (2 * arb(hstep))),
                                                                  F(Lu[2] - 2 * Lu[0] ** 2)]
    res['dL/ds at 0: difference, Lambda\'\' - 2 Lambda^2 - Lambda\'^2/Lambda at c'] = [
        F((vals[1][1] - vals[-1][1]) / (2 * arb(hstep))), F(Lc0[2] - 2 * Lc0[0] ** 2 - Lc0[1] ** 2 / Lc0[0])]
    res['trajectory (s, m, R, eta_dir, eta, W, L), every sixth point'] = [
        [x['s'], x['m'], x['R'], x['eta_dir'], x['eta'], x['W'], x['L']] for x in hist[::6]]
    return res


if __name__ == '__main__':
    nproc = int(os.environ.get('SLURM_CPUS_PER_TASK', '1'))
    d = json.load(open('dbn_collisions.json'))
    rows = d[[k for k in d if k.startswith('collisions among')][0]]

    def pr(x):
        return [int(i) for i in (json.loads(x) if isinstance(x, str) else x)]

    cols = [{'pair': pr(r['pair']), 'z': r['z'], 't': r['t_coll'], 'm0': r['m(0)'], 'Lc': r['Lambda(c)']} for r in rows]
    other = d['summary']['other collisions (pairs of nonadjacent indices, or k > 200)']
    allc = [(c_['pair'], c_['z'], c_['t']) for c_ in cols] + [(pr(o[0]), o[1], o[2]) for o in other]
    ctx.prec = 128
    gam = [float(z.imag.mid()) for z in acb.zeta_zeros(1, 2000)]
    ctx.prec = PREC
    gam_tail = gam[206:]
    rec('pairs', len(cols))
    rec('nonreal pairs followed', len(allc))
    rec('tail zeros (count, first, last ordinate)', [len(gam_tail), gam_tail[0], gam_tail[-1]])
    pool = mpr.get_context('fork').Pool(nproc)
    order = sorted(range(len(cols)), key=lambda i: float(cols[i]['t']))
    jobs = [(cols[i], (gam[cols[i]['pair'][0] - 1], gam[cols[i]['pair'][1] - 1]), allc, gam_tail) for i in order]
    results = pool.map(run_pair, jobs, chunksize=1)
    rec('results', results)
    summ = []
    for r in results:
        if 'm(0)' not in r:
            continue
        low = r['lower 1/((1 + max(eta_dir, 0)) m(0))']
        up = r['upper s* (rho, -Q)']
        summ.append({'pair': r['pair'], 't_c': r['refined t_c'], 'm(0)': r['m(0)'], 'R(0)': r['R(0)'], 'rho': r['rho = sup R'],
                     'eta_dir max': r['largest eta_dir'], 'eta max': r['largest eta (B form)'], 'W max': r['largest W'],
                     'lower': low, 'lower (eta)': r['lower 1/((1 + eta) m(0))'], 'upper (Q)': up,
                     'upper (W)': r['upper s* (rho, W)'], 's_c': r['s_c'],
                     'in bracket': bool(low <= r['s_c'] and (up is None or r['s_c'] <= up))})
    rec('summary', summ)
    ok = [r for r in results if 'm(0)' in r]
    sets = {'k <= 120, t > -1.6': [r for r in ok if r['pair'][0] <= 120 and r['s_c'] < 1.6],
            'k <= 200, t > -1.6': [r for r in ok if r['s_c'] < 1.6], 'all': ok}
    stats = {}
    for name, rs in sets.items():
        if not rs:
            continue
        early = [1 - 1 / (r['m(0)'] * r['s_c']) for r in rs]
        ferr = [abs(r['frozen R predictor s*(R(0))'] / r['s_c'] - 1) for r in rs if r['frozen R predictor s*(R(0))']]
        foerr = [abs(r['first order predictor g^2/8 + g^4 R(0)/32'] / r['s_c'] - 1) for r in rs]
        fr = [(abs(r['frozen R predictor s*(R(0))'] / r['s_c'] - 1), 6 * r['R(0)'] / r['m(0)'], r['pair']) for r in rs
              if r['frozen R predictor s*(R(0))']]
        by_th = {}
        for th in (0.3, 0.4, 0.5, 0.6, 0.7, 0.8):
            sel = [x for x in fr if x[1] < th]
            by_th['6R(0)/m(0) < %g' % th] = [len(sel), max(x[0] for x in sel) if sel else None]
        stats[name] = {'pairs': len(rs), 'mean of 1 - 1/(m(0) s_c) (how early 1/m(0) is)': sum(early) / len(early),
                       'range of s_c m(0)': [min(r['m(0)'] * r['s_c'] for r in rs), max(r['m(0)'] * r['s_c'] for r in rs)],
                       'range of s_c L(0)': [min(r['L(0) = Lambda(c)'] * r['s_c'] for r in rs),
                                             max(r['L(0) = Lambda(c)'] * r['s_c'] for r in rs)],
                       'range of s_c/(g^2/8)': [min(8 * r['s_c'] / r['gap g at s = 0 (z units)'] ** 2 for r in rs),
                                                max(8 * r['s_c'] / r['gap g at s = 0 (z units)'] ** 2 for r in rs)],
                       'frozen R predictor: count and largest error by 6R(0)/m(0)': by_th,
                       'frozen R predictor: largest error, 6R(0)/m(0), pair': max(fr) if fr else None,
                       'frozen R predictor: pairs with 6R(0)/m(0) >= r1 (undefined)': len(rs) - len(fr),
                       'largest frozen R predictor error': max(ferr) if ferr else None,
                       'largest first order predictor error': max(foerr),
                       'all in the zero free bracket': all(x['in bracket'] for x in summ if x['pair'] in [r['pair'] for r in rs])}
    rec('statistics', stats)
    out['seconds'] = time.time() - t0
    json.dump(out, open('wigbiane_flow.json', 'w'), indent=1)
    print('done', round(time.time() - t0, 1))
