# pakes_filter.py: the infinite divisibility and GGC (Hankel) tests of the criticality filter, Section 11.5.5
# (label pk:sec:filter), for the variable of Biane, Pitman and Yor, its Gaussian tilts and the variables of
# Dirichlet L functions and of the Davenport and Heilbronn function: Tables pk:tab:tilt and pk:tab:chars and the
# survey of real characters in the text. Added for the book in October 2026, from Part D of the triage of the first
# author's chapter on Pakes and the Thorin condition (misc/reviews/nick-pakes/triage-D.md), recomputed independently.
#
# Variables. For a completed function Lambda(s) = Lambda(1 - s) with Lambda(1/2 + s) = int K(y) e^{sy} dy, K even
# and positive, the variable B with E B^{s/2} = Lambda(s)/Lambda(0) has log B = 2Y, where Y has the density
# e^{-y/2} K(y)/Lambda(0). For xi, K is the Polya kernel Phi of eq:ch2:Phi, Lambda = 2 xi and B = pi S_2/2
# (Theorem rz:thm:bpy). For a real sequence a of period q and parity p (a_1 = 1) whose theta series
# theta(x) = sum_{n>=1} n^p a_n e^{-pi n^2 x/q} satisfies theta(1/x) = x^{p+1/2} theta(x), K(y) = 2 e^{(p+1/2)y} theta(e^{2y})
# as in Section sec:ch16:eta, and Lambda(s) = (q/pi)^{(s+p)/2} Gamma((s+p)/2) sum a_n n^{-s}. Side 'fwd' is B, side
# 'rev' is 1/B (density of -Y). The tilt B_t has the law of B tilted by e^{t Y^2/4} = e^{t (log B)^2/16}.
#
# Tests. kappa_k(lam) are the cumulants of the Esscher transform e^{-lam b} P(B in db)/E e^{-lam B}. B is infinitely
# divisible iff kappa_k(lam) >= 0 for all k >= 1 and lam > 0. If B is a GGC, the numbers kappa_{k+1}(lam)/k! are
# moments of a positive measure, so the Hankel matrices C_n = [kappa_{i+j+1}/(i+j)!] and C'_n = [kappa_{i+j+2}/(i+j+1)!],
# i, j < n, are positive semidefinite. Entries are divided by sc^{i+j} (C_n) and sc^{i+j+1} (C'_n), sc = kappa_2/kappa_1,
# a congruence that keeps the inertia; the least eigenvalues printed in the book are those of the scaled 7 x 7 matrices.
# The inertia of the leading blocks of order n <= 12 (n <= 10 in the survey) is read from the pivots of Gaussian
# elimination (Jacobi's rule).
#
# Method. Moments m_k(lam) = int e^{2ky} e^{-lam e^{2y}} e^{t y^2/4} e^{-+y/2} K(y) dy by the trapezoid rule on [-5, 5]
# (the kernels are analytic in |Im y| < pi/4 and below 10^-700 beyond |y| = 5), at two settings, 150 digits with step
# 1/64 and 200 digits with step 1/96, which must agree. Exact checks: at t = 0, B (xi) has Thorin measure
# 2 sum delta_{pi n^2} and the reflected chi_4 variable Z of Remark rem:dbeta:exit has Thorin measure
# sum delta_{pi (2m-1)^2/4}, so kappa_{k+1}(lam)/k! = 2 sum (lam + pi n^2)^{-k-1}, resp. sum (lam + pi(2m-1)^2/4)^{-k-1}.
# Independent check at lam = 0, with no kernel and no quadrature: E B^k = Lambda(2k)/Lambda(0), E B^{-k} = Lambda(1+2k)/Lambda(0).
# Survey: all fundamental discriminants 3 <= |d| <= 40 (positivity of the kernel, and the test to order 10 at
# lam in {0, 1/2, 2, 10, 20}, both sides, 80 digits, step 1/64) and the self dual family through L(s, chi_{-7}),
#   f_r = (1 - 2 r cos(phi/2)) L(s, chi_{-7}) + r e^{-i phi/2} L(s, chi) + r e^{i phi/2} L(s, chi-bar),  e^{i phi} = eps(chi),
# chi the character mod 7 with chi(3) = e^{2 pi i/6}, which has real coefficients, a_1 = 1 and Lambda(s) = Lambda(1 - s).
# Positivity of a kernel: theta(x) = e^{-pi x/q} (1 + r(x)) with r(x) = sum_{n>=2} n^p a_n e^{-pi (n^2-1) x/q}, and by the
# functional equation theta > 0 on (0, oo) iff 1 + r > 0 on [1, oo); the script reports S = sum_{n>=2} n^p |a_n| e^{-pi(n^2-1)/q}
# (S < 1 proves positivity, Lemma pk:lem:thetapos; certified in pakes_filter_cert.py) and min(1 + r) on [1, 60], step 1/200.
#
# Run on Hopper: bash misc/tools/hopper_run.sh -c 32 -g pakes_filter.json code/pakes_filter.py -- 32
# (python3 pakes_filter.py [ncpu]; writes pakes_filter.json in the current folder; about 15 minutes on 32 cores).
import json
import sys
import time
from multiprocessing import Pool

import mpmath as mp

SETTINGS = [(150, 64), (200, 96)]          # (digits, 1/h)
YMAX = 5
KMAX = 26                                  # cumulants kappa_1..kappa_26, a_0..a_25
NMAX = 12
LAMS = ['0', '0.5', '2', '10', '20', '50', '200']
TS = ['-20', '-10', '-5', '-2', '-1', '-0.5', '-0.2', '-0.1', '-0.05', '-0.02', '-0.01', '0',
      '0.01', '0.02', '0.05', '0.1', '0.2', '0.5', '1', '5']
TS_REV = ['-1', '0', '0.2']
NAMED = ['chi-3', 'chi-4', 'chi5', 'chi8', 'chi12', 'chi24', 'chi-7', 'chi-8', 'chi13', 'DH']
SURVEY_LAMS = ['0', '0.5', '2', '10', '20']
SURVEY_NMAX = 10
SURVEY_KMAX = 22
SURVEY_SETTING = (80, 64)
Q7_RS = ['-3', '-2', '-1.5', '-1', '-0.75', '-0.5', '-0.3', '-0.15', '0', '0.15', '0.3', '0.5', '0.75', '1',
         '1.5', '2', '3']


def kron(d, n):
    # Kronecker symbol (d/n), n >= 1
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


def fundamental(d):
    def squarefree(m):
        m = abs(m)
        p = 2
        while p * p <= m:
            if m % (p * p) == 0:
                return False
            p += 1
        return True
    if d % 4 == 1:
        return squarefree(d)
    if d % 4 == 0:
        return (d // 4) % 4 in (2, 3) and squarefree(d // 4)
    return False


SURVEY_DS = [d for d in range(-40, 41) if abs(d) >= 3 and fundamental(d)]


def q7_coeffs(r):
    # f_r above; chi(3^k) = w^k with w = e^{2 pi i/6}, eps = tau(chi)/(i sqrt 7), chi_{-7} = Legendre symbol mod 7
    w = mp.exp(2j * mp.pi / 6)
    chi = {int(pow(3, k, 7)): w ** k for k in range(6)}
    tau = mp.fsum(chi[a] * mp.exp(2j * mp.pi * a / 7) for a in range(1, 7))
    eps = tau / (1j * mp.sqrt(7))
    c = mp.mpf(r) * mp.exp(-1j * mp.arg(eps) / 2)
    a0 = 1 - 2 * mp.re(c)
    return {a: a0 * kron(-7, a) + 2 * mp.re(c * chi[a]) for a in range(1, 7)}


def case(name):
    # (q, {residue: coefficient}, parity) for a theta kernel, or None for xi
    if name == 'xi':
        return None
    if name == 'DH':
        kap = (mp.sqrt(10 - 2 * mp.sqrt(5)) - 2) / (mp.sqrt(5) - 1)
        return 5, {1: mp.mpf(1), 2: kap, 3: -kap, 4: mp.mpf(-1)}, 1
    if name.startswith('q7r'):
        a = q7_coeffs(name[3:])
        return 7, {k: a[k] for k in range(1, 7)}, 1
    d = int(name[3:])
    q = abs(d)
    return q, {r: mp.mpf(kron(d, r)) for r in range(1, q + 1) if kron(d, r) != 0}, (1 if d < 0 else 0)


def kernel(cs, y, dps):
    # the even kernel K at real y (by the functional equation K(-y) = K(y))
    y = abs(y)
    x = mp.exp(2 * y)
    s = mp.mpf(0)
    n = 1
    if cs is None:
        e9, e5 = mp.exp(9 * y / 2), mp.exp(5 * y / 2)
        while True:
            if n > 1 and mp.pi * (n * n - 1) * x > 2.4 * dps + 60:
                break
            s += (2 * mp.pi ** 2 * n ** 4 * e9 - 3 * mp.pi * n ** 2 * e5) * mp.exp(-mp.pi * n * n * x)
            n += 1
        return 2 * s
    q, a, p = cs
    while True:
        if n > 1 and mp.pi * (n * n - 1) * x / q > 2.4 * dps + 60:
            break
        r = n % q
        if r in a:
            s += (n ** p) * a[r] * mp.exp(-mp.pi * n * n * x / q)
        n += 1
    return 2 * mp.exp((p + mp.mpf(1) / 2) * y) * s


def cumulants(mu, K):
    c = [mp.mpf(0)] * (K + 1)
    for n in range(1, K + 1):
        c[n] = mu[n] - mp.fsum(mp.binomial(n - 1, i - 1) * c[i] * mu[n - i] for i in range(1, n))
    return c


def pivots(M, n):
    # Gaussian elimination without pivoting: the k-th pivot is D_k/D_{k-1}; by Jacobi's rule the number of negative
    # pivots among the first k is the number of negative eigenvalues of the leading k x k block
    M = [row[:] for row in M]
    piv = []
    for k in range(n):
        p = M[k][k]
        piv.append(p)
        if p == 0:
            break
        for i in range(k + 1, n):
            f = M[i][k] / p
            for j in range(k, n):
                M[i][j] -= f * M[k][j]
    return piv


def analyse(a, nmax):
    sc = a[1] / a[0]
    out = {}
    for shift, key in ((0, 'C'), (1, 'Cp')):
        M = [[a[i + j + shift] / sc ** (i + j + shift) for j in range(nmax)] for i in range(nmax)]
        piv = pivots(M, nmax)
        inert = [sum(1 for p in piv[:k] if p < 0) for k in range(1, len(piv) + 1)]
        out[key + '_inertia'] = inert
        out[key + '_first_fail'] = next((k + 1 for k, v in enumerate(inert) if v > 0), None)
        nf = out[key + '_first_fail']
        if nf is not None:
            E = mp.eigsy(mp.matrix([[M[i][j] for j in range(nf)] for i in range(nf)]), eigvals_only=True)
            out[key + '_mineig_first_fail'] = mp.nstr(min(E), 12)
        E = mp.eigsy(mp.matrix([[M[i][j] for j in range(nmax)] for i in range(nmax)]), eigvals_only=True)
        out[key + '_mineig_nmax'] = mp.nstr(min(E), 6)
        if nmax >= 7:
            E = mp.eigsy(mp.matrix([[M[i][j] for j in range(7)] for i in range(7)]), eigvals_only=True)
            out[key + '_mineig7'] = mp.nstr(min(E), 12)
    return out


def moments_rows(name, side, ts, lams, setting, kmax, nmax):
    dps, hinv = setting
    mp.mp.dps = dps
    cs = case(name)
    h = mp.mpf(1) / hinv
    N = 2 * YMAX * hinv
    ys = [-YMAX + j * h for j in range(N + 1)]
    sg = 1 if side == 'fwd' else -1
    g = [mp.exp(-sg * y / 2) * kernel(cs, y, dps) for y in ys]
    e2 = [mp.exp(2 * y) for y in ys]
    rows = []
    for tstr in ts:
        t = mp.mpf(tstr)
        gt = [g[j] * mp.exp(t * ys[j] ** 2 / 4) for j in range(N + 1)]
        for lamstr in lams:
            lam = mp.mpf(lamstr)
            w = [h * gt[j] * mp.exp(-lam * e2[j]) for j in range(N + 1)]
            m = []
            pw = [mp.mpf(1)] * (N + 1)
            for k in range(kmax + 1):
                m.append(mp.fsum(w[j] * pw[j] for j in range(N + 1)))
                pw = [pw[j] * e2[j] for j in range(N + 1)]
            mu = [x / m[0] for x in m]
            c = cumulants(mu, kmax)
            a = [c[k + 1] / mp.factorial(k) for k in range(kmax)]
            r = {'kernel': name, 'side': side, 't': tstr, 'lam': lamstr, 'dps': dps, 'hinv': hinv,
                 'negcum': [k for k in range(1, kmax + 1) if c[k] < 0], 'kappa1': mp.nstr(c[1], 15),
                 'end_weight': mp.nstr(max(abs(w[0]), abs(w[-1])) / abs(m[0]), 3)}
            r.update(analyse(a, nmax))
            r['_a'] = [mp.nstr(x, 140) for x in a[:2 * nmax + 1]]
            rows.append(r)
    return rows


def task(args):
    kind = args[0]
    if kind == 'main':
        _, name, side, ts, setting = args
        return moments_rows(name, side, ts, LAMS, setting, KMAX, NMAX)
    _, name, side = args
    return moments_rows(name, side, ['0'], SURVEY_LAMS, SURVEY_SETTING, SURVEY_KMAX, SURVEY_NMAX)


def positivity(name):
    mp.mp.dps = 30
    q, a, p = case(name)
    NT = int(mp.sqrt(250 * q / mp.pi)) + 3          # e^{-pi (n^2 - 1)/q} < e^{-240} beyond NT
    def r(x):
        return mp.fsum((n ** p) * a[n % q] * mp.exp(-mp.pi * (n * n - 1) * x / q) for n in range(2, NT) if n % q in a)
    S = mp.fsum((n ** p) * abs(a[n % q]) * mp.exp(-mp.pi * (n * n - 1) / q) for n in range(2, NT) if n % q in a)
    vals = [(1 + r(1 + mp.mpf(j) / 200), 1 + mp.mpf(j) / 200) for j in range(0, 200 * 59 + 1)]
    mn = min(vals)
    return {'S': mp.nstr(S, 6), 'min_1_plus_r': mp.nstr(mn[0], 6), 'argmin': mp.nstr(mn[1], 6)}


def exact_a(kind, lam, dps):
    # kappa_{k+1}(lam)/k! for k < KMAX from the Thorin measure: sum_z c_z (lam + z)^{-k-1}. At lam = 0 by zeta values,
    # 2 zeta(2k+2)/pi^{k+1} (xi) and (4/pi)^{k+1} (1 - 2^{-2k-2}) zeta(2k+2) (Z); for lam > 0 as (-1)^k S^{(k)}(lam)/k!,
    # S(x) = sum_z c_z/(x + z) in closed form, by mpmath's finite differences at raised working precision
    mp.mp.dps = dps + 40
    if lam == '0':
        if kind == 'xi':
            out = [2 * mp.zeta(2 * k + 2) / mp.pi ** (k + 1) for k in range(KMAX)]
        else:
            out = [(4 / mp.pi) ** (k + 1) * (1 - mp.mpf(2) ** (-2 * k - 2)) * mp.zeta(2 * k + 2) for k in range(KMAX)]
        mp.mp.dps = dps
        return out
    lam = mp.mpf(lam)
    if kind == 'xi':
        def S(x):            # 2 sum_n 1/(x + pi n^2) = (pi c coth(pi c) - 1)/(pi c^2), c^2 = x/pi
            c = mp.sqrt(x / mp.pi)
            return (mp.pi * c * mp.coth(mp.pi * c) - 1) / (mp.pi * c * c)
    else:
        def S(x):            # sum_m 1/(x + pi (2m-1)^2/4) = (d/dx) log cosh sqrt(pi x)
            b = mp.sqrt(mp.pi * x)
            return mp.tanh(b) * mp.pi / (2 * b)
    co = mp.taylor(S, lam, KMAX - 1)
    out = [(-1) ** k * co[k] for k in range(KMAX)]
    mp.mp.dps = dps
    return out


def exact_task(args):
    kind, lamstr = args
    return [mp.nstr(x, 140) for x in exact_a(kind, lamstr, SETTINGS[0][0])]


def Lambda_fn(name):
    # completed function with Lambda(s) = Lambda(1 - s); for xi, Lambda = 2 xi with Lambda(0) = 1
    if name == 'xi':
        def Lam(s):
            return s * (s - 1) * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)
        return Lam, mp.mpf(1)
    q, a, p = case(name)
    def L(s):
        return mp.mpf(q) ** (-s) * mp.fsum(a[r] * mp.zeta(s, mp.mpf(r) / q) for r in a)
    def Lam(s):
        return (mp.mpf(q) / mp.pi) ** ((s + p) / 2) * mp.gamma((s + p) / 2) * L(s)
    if p == 1:
        L0 = mp.fsum(a[r] * (mp.mpf(1) / 2 - mp.mpf(r) / q) for r in a)   # zeta(0, x) = 1/2 - x
        lam0 = mp.sqrt(q) * L0
    else:
        L1 = -mp.fsum(a[r] * mp.digamma(mp.mpf(r) / q) for r in a) / q      # L(1), sum a_r = 0
        lam0 = mp.sqrt(q) * L1                                                # Lambda(0) = Lambda(1)
    return Lam, lam0


def lvalue_check(name):
    mp.mp.dps = 80
    Lam, lam0 = Lambda_fn(name)
    K = 2 * 7 + 1
    out = {'Lambda0': mp.nstr(lam0, 20)}
    s = mp.mpc('0.31', '1.7')
    out['FE_rel_err'] = mp.nstr(abs(Lam(s) - Lam(1 - s)) / abs(Lam(s)), 3)
    for side in ('fwd', 'rev'):
        if side == 'fwd':
            m = [mp.mpf(1)] + [Lam(mp.mpf(2 * k)) / lam0 for k in range(1, K + 1)]
        else:
            m = [mp.mpf(1)] + [Lam(mp.mpf(1 + 2 * k)) / lam0 for k in range(1, K + 1)]
        c = cumulants(m, K)
        a = [c[k + 1] / mp.factorial(k) for k in range(K)]
        sc = a[1] / a[0]
        res = {}
        for shift, key in ((0, 'C'), (1, 'Cp')):
            ev = []
            for n in range(2, 8):
                M = mp.matrix([[a[i + j + shift] / sc ** (i + j + shift) for j in range(n)] for i in range(n)])
                ev.append(mp.nstr(min(mp.eigsy(M, eigvals_only=True)), 12))
            res[key + '_mineig_orders_2_to_7'] = ev
        res['negcum'] = [k for k in range(1, K + 1) if c[k] < 0]
        if side == 'fwd':
            res['kappa_2_3_4'] = [mp.nstr(c[k], 15) for k in (2, 3, 4)]
            res['3k3sq_minus_2k2k4'] = mp.nstr(3 * c[3] ** 2 - 2 * c[2] * c[4], 15)
        else:
            res['kappa_2_3_4'] = [mp.nstr(c[k], 15) for k in (2, 3, 4)]
            res['3k3sq_minus_2k2k4'] = mp.nstr(3 * c[3] ** 2 - 2 * c[2] * c[4], 15)
        out[side] = res
    return out


def q7_check(r):
    mp.mp.dps = 30
    q, a, p = case('q7r' + r)
    Lam, lam0 = Lambda_fn('q7r' + r)
    s = mp.mpc('0.3', '2.1')
    fe = abs(Lam(s) - Lam(1 - s)) / abs(Lam(s))
    return {'coeffs': [mp.nstr(a[k], 10) for k in range(1, 7)], 'FE_rel_err': mp.nstr(fe, 3)}


if __name__ == '__main__':
    ncpu = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    if 'quick' in sys.argv:          # a smoke test of the whole pipeline at low precision
        SETTINGS = [(40, 16), (50, 24)]
        TS = ['-5', '-1', '0', '0.2', '1', '5']
        LAMS = ['0', '2']
        NAMED = ['chi-4', 'DH']
        SURVEY_DS = [-19, -4, 5]
        Q7_RS = ['0', '3']
        SURVEY_SETTING = (30, 16)
    t0 = time.time()
    out = {'settings': SETTINGS, 'lams': LAMS, 'ts': TS, 'survey_lams': SURVEY_LAMS}
    tasks = []
    for st in SETTINGS:
        for t in TS:
            tasks.append(('main', 'xi', 'fwd', [t], st))
        tasks.append(('main', 'xi', 'rev', TS_REV, st))
        for nm in NAMED:
            for side in ('fwd', 'rev'):
                tasks.append(('main', nm, side, ['0'], st))
    survey_names = ['chi%d' % d for d in SURVEY_DS]
    q7_names = ['q7r' + r for r in Q7_RS]
    with Pool(ncpu) as P:
        pos = dict(zip(survey_names + q7_names, P.map(positivity, survey_names + q7_names)))
        stasks = [('survey', nm, side) for nm in survey_names + q7_names
                  if mp.mpf(pos[nm]['min_1_plus_r']) > 0 for side in ('fwd', 'rev')]
        res = P.map(task, tasks + stasks, chunksize=1)
        lv = dict(zip(['xi'] + NAMED, P.map(lvalue_check, ['xi'] + NAMED)))
        q7c = dict(zip(Q7_RS, P.map(q7_check, Q7_RS)))
        ekeys = [(kind, lamstr) for lamstr in LAMS for kind in ('xi', 'Z')]
        exact = dict(zip(ekeys, P.map(exact_task, ekeys)))
    rows = [r for o in res[:len(tasks)] for r in o]
    srows = [r for o in res[len(tasks):] for r in o]
    print('moments done, %.0f s' % (time.time() - t0), flush=True)

    # agreement of the two settings, and exact checks at t = 0
    def key(r):
        return (r['kernel'], r['side'], r['t'], r['lam'])
    by = {}
    for r in rows:
        by.setdefault(key(r), {})[(r['dps'], r['hinv'])] = r
    agree = []
    worst = mp.mpf(0)
    mp.mp.dps = SETTINGS[0][0]
    for k, d in by.items():
        r1, r2 = d[tuple(SETTINGS[0])], d[tuple(SETTINGS[1])]
        same = all(r1[x] == r2[x] for x in ('C_inertia', 'Cp_inertia', 'negcum'))
        a1 = [mp.mpf(x) for x in r1['_a']]
        a2 = [mp.mpf(x) for x in r2['_a']]
        rel = max(abs(x - y) / abs(y) for x, y in zip(a1, a2))
        worst = max(worst, rel)
        agree.append({'case': list(k), 'inertia_and_signs_agree': same, 'max_rel_diff_a': mp.nstr(rel, 3)})
    out['settings_agree_all'] = all(x['inertia_and_signs_agree'] for x in agree)
    out['settings_max_rel_diff_a'] = mp.nstr(worst, 3)
    checks = []
    for lamstr in LAMS:
        for kind, nm, side in (('xi', 'xi', 'fwd'), ('Z', 'chi-4', 'rev')):
            mp.mp.dps = 150
            ex = [mp.mpf(x) for x in exact[(kind, lamstr)]]
            r = by[(nm, side, '0', lamstr)][tuple(SETTINGS[0])]
            mp.mp.dps = 150
            a = [mp.mpf(x) for x in r['_a']]
            rel = max(abs(a[k] - ex[k]) / abs(ex[k]) for k in range(len(a)))
            checks.append({'kind': kind, 'lam': lamstr, 'max_rel_err_a': mp.nstr(rel, 3)})
    out['exact_checks'] = checks
    for r in rows + srows:
        r.pop('_a', None)
    out['rows'] = [r for r in rows if (r['dps'], r['hinv']) == tuple(SETTINGS[0])]
    out['rows_setting2'] = [r for r in rows if (r['dps'], r['hinv']) == tuple(SETTINGS[1])]
    out['positivity'] = pos
    out['survey_rows'] = srows
    out['lvalue_check'] = lv
    out['q7_family'] = q7c

    # digests for the text
    def row(nm, side, t, lam):
        return by[(nm, side, t, lam)][tuple(SETTINGS[0])]
    tab = []
    for (t, lam) in (('0', '2'), ('-1', '2'), ('-1', '0'), ('-5', '2'), ('0.2', '2'), ('0.2', '0'), ('1', '2'), ('5', '2')):
        r = row('xi', 'fwd', t, lam)
        tab.append({'t': t, 'lam': lam, 'C7_mineig': r['C_mineig7'], 'Cp7_mineig': r['Cp_mineig7'],
                    'C_first_fail': r['C_first_fail'], 'Cp_first_fail': r['Cp_first_fail']})
    out['table_tilt'] = tab
    ff = {}
    for t in TS:
        ff[t] = {lam: [row('xi', 'fwd', t, lam)['C_first_fail'], row('xi', 'fwd', t, lam)['Cp_first_fail']] for lam in LAMS}
    out['tilt_first_fail_C_Cp'] = ff
    chars = {}
    for nm in ['xi'] + NAMED:
        chars[nm] = {}
        for side in ('fwd', 'rev'):
            t = '0'
            rs = {lam: by[(nm, side, t, lam)][tuple(SETTINGS[0])] for lam in LAMS}
            chars[nm][side] = {'Cp7_mineig_lam0': rs['0']['Cp_mineig7'], 'Cp7_mineig_lam2': rs['2']['Cp_mineig7'],
                               'first_fail_C_Cp_by_lam': {lam: [rs[lam]['C_first_fail'], rs[lam]['Cp_first_fail']] for lam in LAMS}}
    out['table_chars'] = chars
    out['rows_without_negative_cumulant'] = sum(1 for r in out['rows'] if not r['negcum'])
    out['rows_total'] = len(out['rows'])
    surv = {}
    for nm in survey_names + q7_names:
        e = {'positivity': pos[nm]}
        for side in ('fwd', 'rev'):
            rs = sorted([r for r in srows if r['kernel'] == nm and r['side'] == side], key=lambda r: float(r['lam']))
            if rs:
                e[side] = [r['Cp_first_fail'] for r in rs]
                e[side + '_C'] = [r['C_first_fail'] for r in rs]
                e[side + '_negcum'] = sorted(set(k for r in rs for k in r['negcum']))
        surv[nm] = e
    out['survey_digest'] = surv
    out['time_s'] = round(time.time() - t0)
    with open('pakes_filter.json', 'w') as f:
        json.dump(out, f, indent=1)

    print('settings agree (inertia, cumulant signs):', out['settings_agree_all'], ' max rel diff of a_k:', out['settings_max_rel_diff_a'])
    for c in checks:
        print('exact check', c)
    print('rows without a negative cumulant: %d of %d' % (out['rows_without_negative_cumulant'], out['rows_total']))
    print('tilt table (t, lam, min eig C_7, C\'_7, first failing order C, C\'):')
    for x in tab:
        print('  ', x)
    print('tilt: first failing order of (C, C\') by lam', LAMS)
    for t in TS:
        print('  t = %6s: %s' % (t, ' '.join('%s/%s' % tuple('-' if v is None else str(v) for v in ff[t][lam]) for lam in LAMS)))
    print('characters, t = 0: min eig C\'_7 at lam = 0 and 2, first failing order of C\' by lam')
    for nm in ['xi'] + NAMED:
        for side in ('fwd', 'rev'):
            e = chars[nm][side]
            print('  %-6s %s %-20s %-20s %s' % (nm, side, e['Cp7_mineig_lam0'], e['Cp7_mineig_lam2'],
                                               ' '.join('-' if v[1] is None else str(v[1]) for v in e['first_fail_C_Cp_by_lam'].values())))
    print('L value check at lam = 0 (C\'_n, n = 2..7):')
    for nm, v in lv.items():
        print('  %-6s Lambda(0) = %s FE %s fwd %s | rev %s' % (nm, v['Lambda0'], v['FE_rel_err'], v['fwd']['Cp_mineig_orders_2_to_7'][-1], v['rev']['Cp_mineig_orders_2_to_7'][-1]))
    print('survey (first failing order of C\' at lam in %s, fwd | rev):' % SURVEY_LAMS)
    for nm in survey_names + q7_names:
        e = surv[nm]
        if 'fwd' in e:
            print('  %-8s S = %-9s min(1+r) = %-9s %s | %s' % (nm, e['positivity']['S'], e['positivity']['min_1_plus_r'],
                  ' '.join('-' if v is None else str(v) for v in e['fwd']), ' '.join('-' if v is None else str(v) for v in e['rev'])))
        else:
            print('  %-8s S = %-9s min(1+r) = %-9s kernel not positive' % (nm, e['positivity']['S'], e['positivity']['min_1_plus_r']))
    print('q7 family checks:', {r: v['FE_rel_err'] for r, v in q7c.items()})
    print('time %.0f s' % (time.time() - t0))
