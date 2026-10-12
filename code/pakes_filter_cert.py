# pakes_filter_cert.py: ball arithmetic certificates (Arb, through python-flint) for Section 11.5.5 (label pk:sec:filter):
# the starred entries of Tables pk:tab:tilt and pk:tab:chars, the four cumulant inequality of Lemma pk:lem:fourcum for
# chi_{-8} (Proposition pk:prop:chi8), the first order argument of Proposition pk:prop:bpytilt (two values of G') and the
# positivity bounds S_a < 1 of Lemma pk:lem:thetapos. Added for the book in October 2026, from Part D of the triage of the
# first author's chapter on Pakes (misc/reviews/nick-pakes/triage-D.md); the numerical survey is pakes_filter.py.
#
# Hankel certificates. For a GGC the numbers a_k = kappa_{k+1}(lam)/k!, kappa_k(lam) the cumulants of the Esscher transform
# e^{-lam b} P(B in db)/E e^{-lam B}, are moments of a positive measure, so K_n = [a_{i+j}] and K'_n = [a_{i+j+1}], i, j < n,
# are positive semidefinite, and so are the congruent rescalings [a_{i+j+s}/sc^{i+j+s}], sc > 0 rational. The script encloses
# the moments in balls, computes the cumulants and the rescaled matrix in ball arithmetic, takes an exact rational vector v
# close to the eigenvector of the least eigenvalue of the midpoint matrix, and encloses v^T M v; an enclosure inside
# (-oo, 0) proves that M is not positive semidefinite, hence that the variable is not a GGC.
# Moments, two routes. (1) At lam = 0 from completed L values, with no kernel and no quadrature: E B^k = Lambda(2k)/Lambda(0)
# and E B^{-k} = Lambda(1+2k)/Lambda(0), Lambda(s) = (q/pi)^{(s+p)/2} Gamma((s+p)/2) q^{-s} sum_r a_r zeta(s, r/q) (Hurwitz),
# Lambda(0) = sqrt(q) L(0) = sqrt(q) sum_r a_r (1/2 - r/q) for odd sequences and Lambda(0) = Lambda(1) = sqrt(q) L(1),
# L(1) = -(1/q) sum_r a_r psi(r/q), for even ones; for xi, Lambda(s) = s(s-1) pi^{-s/2} Gamma(s/2) zeta(s) = 2 xi(s).
# (2) For lam > 0 and for the Gaussian tilts, m_k = int e^{2ky} e^{-lam e^{2y}} e^{t y^2/4} e^{-+y/2} K(y) dy, K the even
# kernel (Phi of eq:ch2:Phi for xi, 2 e^{(p+1/2)y} theta(e^{2y}) for a theta series), by acb.integral on [-5, 0] and [0, 5]:
# on [0, 5] the series in x = e^{2y}, on [-5, 0] the same series at e^{-2y} (K is even by the functional equation), both
# holomorphic for |Im y| < pi/4, where the integrand returns a non-finite ball outside; the series are summed with a rigorous
# tail bound. Tails: for |y| >= 5, k <= 25, -20 <= t <= 5, lam >= 0 and q <= 24 the integrand is at most
# 40 exp(55|y| + 5y^2/4 - (pi/24) e^{2|y|}) <= e^{-2577}, with logarithmic derivative below -5000, so each tail is below
# 2^-3000, which is added to each half integral.
# Run on Hopper: bash misc/tools/hopper_run.sh -c 16 -g pakes_filter_cert.json code/pakes_filter_cert.py -- 16
# (python3 pakes_filter_cert.py [ncpu]; writes pakes_filter_cert.json in the current folder; a few minutes).
import json
import sys
import time
from fractions import Fraction
from math import comb, factorial
from multiprocessing import Pool

PREC = 256
YCUT = 5

# quadrature certificates: (kernel, side, t, lam, order, shift); shift 1 is K'_n, shift 0 is K_n
QUAD_CASES = [
    ('xi', 'fwd', '-20', '2', 8, 1), ('xi', 'fwd', '-5', '2', 7, 1), ('xi', 'fwd', '-5', '2', 7, 0),
    ('xi', 'fwd', '-1', '2', 7, 1), ('xi', 'fwd', '-1/2', '2', 7, 1), ('xi', 'fwd', '-1/10', '2', 8, 1),
    ('xi', 'fwd', '-1/100', '2', 8, 1), ('xi', 'fwd', '-1/100', '0', 9, 1), ('xi', 'fwd', '1/100', '2', 5, 1),
    ('xi', 'fwd', '1/5', '2', 4, 1), ('xi', 'fwd', '1', '2', 4, 1), ('xi', 'fwd', '5', '2', 3, 1),
    ('xi', 'fwd', '0', '2', 7, 1), ('xi', 'fwd', '0', '2', 12, 1),
    ('chi24', 'fwd', '0', '2', 11, 1), ('chi8', 'rev', '0', '2', 7, 1), ('chi12', 'rev', '0', '10', 9, 1),
    ('chi-7', 'fwd', '0', '0', 3, 1), ('DH', 'rev', '0', '0', 2, 1),
]
# L value certificates at lam = 0: (kernel, side, order, shift)
LV_CASES = [
    ('chi-7', 'fwd', 3, 1), ('chi-8', 'fwd', 2, 1), ('chi13', 'fwd', 3, 1), ('DH', 'fwd', 6, 1),
    ('xi', 'rev', 3, 1), ('chi-7', 'rev', 2, 1), ('chi-8', 'rev', 2, 1), ('chi13', 'rev', 3, 1), ('DH', 'rev', 2, 1),
    ('chi24', 'rev', 6, 1), ('chi5', 'rev', 10, 1),
]
FOURCUM = [('chi-8', 'fwd'), ('chi-8', 'rev'), ('chi-7', 'rev'), ('DH', 'rev'), ('xi', 'fwd'), ('chi-4', 'fwd')]
SNAMES = ['DH', 'chi-3', 'chi-4', 'chi5', 'chi8', 'chi12', 'chi24', 'chi-7', 'chi-8', 'chi13']


def setup():
    from flint import ctx
    ctx.prec = PREC


_NF = None


def nonfinite():
    global _NF
    if _NF is None:
        from flint import arb, acb
        for mk in (lambda: acb(arb('nan'), arb('nan')), lambda: acb(arb(1) / arb(0, 1))):
            try:
                x = mk()
                if not x.is_finite():
                    _NF = x
                    break
            except Exception:
                pass
    return _NF


def ball(r):
    from flint import arb
    return r * arb(0, 1)


def kron(d, n):
    if n == 0:
        return 0
    res = 1
    while n % 2 == 0:
        n //= 2
        if d % 2 == 0:
            return 0
        if d % 8 in (3, 5):
            res = -res
    a, m = d % n, n
    while a != 0:
        while a % 2 == 0:
            a //= 2
            if m % 8 in (3, 5):
                res = -res
        a, m = m, a
        if a % 4 == 3 and m % 4 == 3:
            res = -res
        a %= m
    return res if m == 1 else 0


def case_data(name):
    # (q, {residue: arb coefficient}, parity), or None for xi
    from flint import arb
    if name == 'xi':
        return None
    if name == 'DH':
        five = arb(5)
        kap = ((10 - 2 * five.sqrt()).sqrt() - 2) / (five.sqrt() - 1)
        return 5, {1: arb(1), 2: kap, 3: -kap, 4: arb(-1)}, 1
    d = int(name[3:])
    q = abs(d)
    return q, {r: arb(kron(d, r)) for r in range(1, q + 1) if kron(d, r) != 0}, (1 if d < 0 else 0)


def series_acb(cd, x):
    # sum_{n>=1} n^p a_n e^{-pi n^2 x/q} (theta) or sum (2 pi^2 n^4 x^{9/4} - 3 pi n^2 x^{5/4}) e^{-pi n^2 x} (xi, x = e^{2y}),
    # with a rigorous bound for the tail; None if Re x is not known to be positive
    from flint import arb, acb
    PI = arb.pi()
    rho = x.real.lower()
    if not (rho > 0):
        return None
    if cd is None:
        q, p, coef = 1, None, None
    else:
        q, coef, p = cd
    s = acb(0)
    n = 1
    first = None
    X = abs(x).upper()
    while True:
        if cd is None:
            term = n ** 2 * (2 * PI ** 2 * n ** 2 * x - 3 * PI) * (-PI * n * n * x).exp()
        else:
            a = coef.get(n % q)
            term = (n ** p) * a * (-PI * n * n * x / q).exp() if a is not None else acb(0)
        s += term
        if first is None:
            first = abs(term).upper()
        m = n + 1
        # |term_m| <= T_m and T_{m+1}/T_m <= r for all later m
        if cd is None:
            T = m ** 2 * (2 * PI ** 2 * m ** 2 * X + 3 * PI) * (-PI * m * m * rho).exp()
            r = (arb(m + 1) / m) ** 4 * (-PI * (2 * m + 1) * rho).exp()
        else:
            T = (m ** p) * (-PI * m * m * rho / q).exp()
            r = (arb(m + 1) / m) ** p * (-PI * (2 * m + 1) * rho / q).exp()
        if r < arb(0.5) and T < first * arb(2) ** (-(PREC + 40)):
            break
        n += 1
        if n > 5000:
            return None
    tail = 2 * T
    return s + acb(ball(tail), ball(tail))


def K_acb(cd, y):
    # the even kernel at y with Re y >= 0 (the caller reflects): Phi(y) = 2 sum (...) for xi, 2 e^{(p+1/2) y} theta(e^{2y})
    x = (2 * y).exp()
    S = series_acb(cd, x)
    if S is None:
        return None
    if cd is None:
        return 2 * (5 * y / 2).exp() * S      # 2 sum (2 pi^2 n^4 e^{9y/2} - 3 pi n^2 e^{5y/2}) e^{-pi n^2 e^{2y}}
    q, coef, p = cd
    return 2 * ((p + 0.5) * y).exp() * S


def moment_task(args):
    name, side, tq, lq, k, half = args
    setup()
    from flint import arb, acb
    cd = case_data(name)
    t = arb(Fraction(tq).numerator) / Fraction(tq).denominator
    lam = arb(Fraction(lq).numerator) / Fraction(lq).denominator
    sg = 1 if side == 'fwd' else -1

    def f(z, analytic):
        Kz = K_acb(cd, z if half > 0 else -z)
        if Kz is None:
            return nonfinite()
        return (2 * k * z).exp() * (-lam * (2 * z).exp()).exp() * (t * z * z / 4).exp() * (-sg * z / 2).exp() * Kz
    a, b = (0, YCUT) if half > 0 else (-YCUT, 0)
    t0 = time.time()
    I = acb.integral(f, a, b, rel_tol=arb(2) ** (-(PREC - 26)))
    re = I.real + ball(arb(2) ** (-3000))
    return (args, re.str(PREC // 3 + 10, radius=True), time.time() - t0)


def lv_moments(name, side, K):
    # normalised moments E B^k (fwd) or E B^{-k} (rev), k = 0..K, from completed L values
    setup()
    from flint import arb
    pi = arb.pi()
    if name == 'xi':
        def Lam(s):
            s = arb(s)
            return s * (s - 1) * pi ** (-s / 2) * (s / 2).gamma() * s.zeta()
        lam0 = arb(1)
    else:
        q, coef, p = case_data(name)

        def Lam(s):
            s = arb(s)
            L = arb(q) ** (-s) * sum((coef[r] * s.zeta(arb(r) / q) for r in coef), arb(0))
            return (arb(q) / pi) ** ((s + p) / 2) * ((s + p) / 2).gamma() * L
        if p == 1:
            lam0 = arb(q).sqrt() * sum((coef[r] * (arb(1) / 2 - arb(r) / q) for r in coef), arb(0))
        else:
            lam0 = arb(q).sqrt() * (-sum((coef[r] * (arb(r) / q).digamma() for r in coef), arb(0)) / q)
    if side == 'fwd':
        m = [arb(1)] + [Lam(2 * k) / lam0 for k in range(1, K + 1)]
    else:
        m = [arb(1)] + [Lam(1 + 2 * k) / lam0 for k in range(1, K + 1)]
    return m, lam0


def cumulants(m):
    from flint import arb
    K = len(m) - 1
    mu = [x / m[0] for x in m]
    c = [arb(0)] * (K + 1)
    for j in range(1, K + 1):
        c[j] = mu[j] - sum((comb(j - 1, i - 1) * c[i] * mu[j - i] for i in range(1, j)), arb(0))
    return c


def certify(m, n, shift):
    setup()
    from flint import arb
    import mpmath as mp
    c = cumulants(m)
    K = len(m) - 1
    a = [c[j + 1] / factorial(j) for j in range(K)]
    mp.mp.dps = PREC // 3

    def tomp(x):
        return mp.mpf(x.mid().str(PREC // 3, radius=False))
    sc = Fraction(mp.nstr(tomp(a[1]) / tomp(a[0]), 30))
    scb = arb(sc.numerator) / sc.denominator
    M = [[a[i + j + shift] / scb ** (i + j + shift) for j in range(n)] for i in range(n)]
    Mm = mp.matrix([[tomp(M[i][j]) for j in range(n)] for i in range(n)])
    E, Q = mp.eigsy(Mm)
    idx = min(range(n), key=lambda i: E[i])
    v = [Fraction(mp.nstr(Q[i, idx], 40)) for i in range(n)]
    vb = [arb(x.numerator) / x.denominator for x in v]
    qf = arb(0)
    for i in range(n):
        for j in range(n):
            qf += vb[i] * vb[j] * M[i][j]
    vv = sum((x * x for x in vb), arb(0))
    rmax = max(float(M[i][j].rad()) for i in range(n) for j in range(n))
    return {'least_eig_midpoint': mp.nstr(E[idx], 12), 'vMv': qf.str(12, radius=True), 'v_norm2': vv.str(10, radius=True),
            'certified_negative': bool(qf < 0), 'max_entry_radius': '%.2e' % rmax,
            'cumulants_certified_positive': bool(all(x > 0 for x in c[1:])), 'K': K}


def fourcum(name, side):
    m, lam0 = lv_moments(name, side, 4)
    c = cumulants(m)
    D = 3 * c[3] ** 2 - 2 * c[2] * c[4]
    return {'kappa2': c[2].str(12, radius=True), 'kappa3': c[3].str(12, radius=True), 'kappa4': c[4].str(12, radius=True),
            '3k3^2-2k2k4': D.str(12, radius=True), 'certified_positive_(not_GGC)': bool(D > 0),
            'certified_negative': bool(D < 0), 'Lambda0': lam0.str(15, radius=True)}


def lv_task(args):
    name, side, n, shift = args
    m, lam0 = lv_moments(name, side, 2 * n + 1)
    out = certify(m, n, shift)
    out['Lambda0'] = lam0.str(15, radius=True)
    return (args, out)


def perturb():
    # G' at two points of the first gap (pi, 4 pi), Proposition pk:prop:bpytilt
    setup()
    from flint import arb
    pi = arb.pi()
    g = arb.const_euler()
    res = {}
    for xs in ('1.015', '2.47'):
        fx = Fraction(xs)
        x = pi * arb(fx.numerator) / fx.denominator
        M = 1
        A = arb(0)
        Ad = arb(0)
        for m in range(1, M + 1):
            s = x - pi * m * m
            am = 4 * pi ** 2 * m ** 4
            bm = 6 * pi * m ** 2
            A += am * (s.log() - 1 + g) / s ** 2 + bm * (s.log() + g) / s
            Ad += am * (3 - 2 * g - 2 * s.log()) / s ** 3 + bm * (1 - g - s.log()) / s ** 2
        A *= pi / 8
        Ad *= pi / 8
        r = (pi * x).sqrt()
        phi = r ** 2 / r.sin() ** 2
        phid = pi * (r.sin() - r * r.cos()) / r.sin() ** 3
        Gd = (Ad * phi - A * phid) / phi ** 2
        res[xs + 'pi'] = {'ImPsi': A.str(15, radius=True), 'phi0': phi.str(15, radius=True), 'G': (A / phi).str(15, radius=True),
                          'Gprime': Gd.str(12, radius=True), 'positive': bool(Gd > 0), 'negative': bool(Gd < 0)}
    return res


def sbound(name):
    # S_a = sum_{n>=2} n^p |a_n| e^{-pi (n^2 - 1)/q}, enclosed (partial sum to 400 and a geometric tail bound)
    setup()
    from flint import arb
    pi = arb.pi()
    q, coef, p = case_data(name)
    S = arb(0)
    N = 400
    for n in range(2, N + 1):
        a = coef.get(n % q)
        if a is not None:
            S += (n ** p) * abs(a) * (-pi * (n * n - 1) / q).exp()
    T = arb(N + 1) ** p * (-pi * ((N + 1) ** 2 - 1) / q).exp()
    r = (arb(N + 2) / (N + 1)) ** p * (-pi * (2 * N + 3) / q).exp()
    S += ball(T / (1 - r))
    return {'S': S.str(10, radius=True), 'S_below_1': bool(S < 1)}


if __name__ == '__main__':
    ncpu = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    setup()
    from flint import arb, acb
    t0 = time.time()
    out = {'prec_bits': PREC}
    I = acb.integral(lambda z, an: (-z * z).exp(), 0, 1, rel_tol=arb(2) ** (-230))
    out['self_test_int_exp'] = I.real.str(25, radius=True)
    out['nonfinite_ball'] = str(nonfinite())
    need = {}
    for (nm, side, t, lam, n, shift) in QUAD_CASES:
        key = (nm, side, t, lam)
        need[key] = max(need.get(key, 0), 2 * n + 1)
    tasks = [(key[0], key[1], key[2], key[3], k, half) for key, K in need.items() for k in range(K + 1) for half in (1, -1)]
    tasks.sort(key=lambda x: -x[4])
    with Pool(ncpu) as P:
        res = P.map(moment_task, tasks, chunksize=1)
        lv = P.map(lv_task, LV_CASES, chunksize=1)
    mom = {}
    tmax = 0
    for (args, sball, dt) in res:
        nm, side, t, lam, k, half = args
        mom.setdefault((nm, side, t, lam), {}).setdefault(k, []).append(arb(sball))
        tmax = max(tmax, dt)
    out['quadrature_integrals'] = len(tasks)
    out['longest_integral_s'] = round(tmax, 1)
    qres = []
    for (nm, side, t, lam, n, shift) in QUAD_CASES:
        key = (nm, side, t, lam)
        ms = [mom[key][k][0] + mom[key][k][1] for k in range(need[key] + 1)]
        r = certify(ms, n, shift)
        r.update({'kernel': nm, 'side': side, 't': t, 'lam': lam, 'order': n, 'matrix': "K'" if shift else 'K',
                  'm0': ms[0].str(15, radius=True)})
        qres.append(r)
    out['quadrature_certificates'] = qres
    # the two routes agree: normalised moments of chi_{-7} (B) and of DH (1/B) at lam = 0, by quadrature and by L values
    agree = {}
    for (nm, side) in (('chi-7', 'fwd'), ('DH', 'rev')):
        key = (nm, side, '0', '0')
        mq = [mom[key][k][0] + mom[key][k][1] for k in range(need[key] + 1)]
        ml, _ = lv_moments(nm, side, need[key])
        agree['%s %s' % (nm, side)] = all((mq[k] / mq[0]).overlaps(ml[k]) for k in range(need[key] + 1))
    out['routes_overlap'] = agree
    lres = []
    for (args, r) in lv:
        nm, side, n, shift = args
        r.update({'kernel': nm, 'side': side, 'lam': '0', 'order': n, 'matrix': "K'" if shift else 'K'})
        lres.append(r)
    out['lvalue_certificates'] = lres
    out['four_cumulant'] = {'%s %s' % (nm, side): fourcum(nm, side) for nm, side in FOURCUM}
    out['perturbation'] = perturb()
    out['S_bounds'] = {nm: sbound(nm) for nm in SNAMES}
    out['time_s'] = round(time.time() - t0)
    with open('pakes_filter_cert.json', 'w') as f:
        json.dump(out, f, indent=1)
    print('self test', out['self_test_int_exp'], 'nonfinite', out['nonfinite_ball'])
    print('%d quadrature integrals, longest %.1f s; quadrature and L value moments overlap: %s' % (len(tasks), tmax, out['routes_overlap']))
    for r in qres:
        print('%-6s %s t=%-7s lam=%-3s %s_%d: v^T M v = %s (midpoint least eigenvalue %s), negative: %s; entry radius <= %s; kappa_1..kappa_%d > 0: %s' % (
            r['kernel'], r['side'], r['t'], r['lam'], r['matrix'], r['order'], r['vMv'], r['least_eig_midpoint'],
            r['certified_negative'], r['max_entry_radius'], r['K'], r['cumulants_certified_positive']))
    for r in lres:
        print('L values %-6s %s lam=0 %s_%d: v^T M v = %s, negative: %s; kappa_1..kappa_%d > 0: %s; Lambda(0) = %s' % (
            r['kernel'], r['side'], r['matrix'], r['order'], r['vMv'], r['certified_negative'], r['K'],
            r['cumulants_certified_positive'], r['Lambda0']))
    for k, v in out['four_cumulant'].items():
        print('four cumulants', k, v)
    for k, v in out['perturbation'].items():
        print('perturbation', k, v)
    for k, v in out['S_bounds'].items():
        print('S bound', k, v)
    print('time %d s' % out['time_s'])
