# pakes_characters.py: the closed forms and checks behind Section 11.5.5 (label pk:sec:filter) that need no Hankel matrices:
# the Davenport and Heilbronn function (constant, root number, functional equation, two zeros off the line, the other self
# dual combination modulo 5), the theta functional equations and positivity bounds of Lemma pk:lem:thetapos, the Mellin
# identities E B^{s/2} = xi(s, a)/xi(0, a), the density of B = pi S_2/2, the reflected transforms of Proposition
# pk:prop:reflected [BPY01, (96)] with the zeros of their numerators, the arcsine mixture of Proposition pk:prop:arcsine, and the
# first order Thorin density of the tilted variable on the gaps (Proposition pk:prop:bpytilt), whose formula for Im psi is
# checked against quadrature along a rotated contour. Added for the book in October 2026, from Part D of the triage of the
# first author's chapter on Pakes (misc/reviews/nick-pakes/triage-D.md). The Hankel tests are in pakes_filter.py and the ball
# arithmetic certificates (two values of G', the bounds S_a < 1) in pakes_filter_cert.py.
#
# Conventions (book): for a real sequence a of period q and parity p (a_1 = 1), theta_a(x) = sum_{n>=1} n^p a_n e^{-pi n^2 x/q},
# xi(s, a) = (q/pi)^{(s+p)/2} Gamma((s+p)/2) L(s, a) = int_0^oo theta_a(x) x^{(s+p)/2} dx/x, kernel K_a(u) = 2 e^{(p+1/2)u}
# theta_a(e^{2u}) with xi(1/2 + s, a) = int K_a(u) e^{su} du; B_a > 0 with E B_a^{s/2} = xi(s, a)/xi(0, a) has the density
# theta_a(b) b^{p/2-1}/xi(0, a), and 1/B_a has the density theta_a(x) x^{(p-1)/2}/xi(0, a).
# Run on Hopper: bash misc/tools/hopper_run.sh -c 8 -g pakes_characters.json code/pakes_characters.py -- 8
# (python3 pakes_characters.py [ncpu]; writes pakes_characters.json in the current folder; a few minutes).
import json
import sys
import time
from multiprocessing import Pool

import mpmath as mp


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


def kappa_dh():
    return (mp.sqrt(10 - 2 * mp.sqrt(5)) - 2) / (mp.sqrt(5) - 1)


def case(name):
    if name == 'DH':
        k = kappa_dh()
        return 5, {1: mp.mpf(1), 2: k, 3: -k, 4: mp.mpf(-1)}, 1
    if name == 'DHminus':
        k = kappa_dh()
        return 5, {1: mp.mpf(1), 2: -1 / k, 3: 1 / k, 4: mp.mpf(-1)}, 1
    d = int(name[3:])
    q = abs(d)
    return q, {r: mp.mpf(kron(d, r)) for r in range(1, q + 1) if kron(d, r) != 0}, (1 if d < 0 else 0)


NAMED = ['chi-3', 'chi-4', 'chi5', 'chi8', 'chi12', 'chi24', 'chi-7', 'chi-8', 'chi13', 'DH']


def theta(cs, x):
    q, a, p = cs
    s = mp.mpf(0)
    n = 1
    while True:
        if n > 1 and mp.pi * (n * n - 1) * x / q > 2.4 * mp.mp.dps + 40:
            break
        r = n % q
        if r in a:
            s += (n ** p) * a[r] * mp.exp(-mp.pi * n * n * x / q)
        n += 1
    return s


def theta_any(cs, x):
    # theta at any x > 0, by the functional equation for x < 1 (sign w = +1)
    q, a, p = cs
    if x >= 1:
        return theta(cs, x)
    return x ** (-(p + mp.mpf(1) / 2)) * theta(cs, 1 / x)


def L(cs, s):
    q, a, p = cs
    return mp.mpf(q) ** (-s) * mp.fsum(a[r] * mp.zeta(s, mp.mpf(r) / q) for r in a)


def xi_a(cs, s):
    q, a, p = cs
    return (mp.mpf(q) / mp.pi) ** ((s + p) / 2) * mp.gamma((s + p) / 2) * L(cs, s)


def Phi(u):
    u = abs(u)
    x = mp.exp(2 * u)
    s = mp.mpf(0)
    n = 1
    while True:
        if n > 1 and mp.pi * (n * n - 1) * x > 2.4 * mp.mp.dps + 40:
            break
        s += (2 * mp.pi ** 2 * n ** 4 * mp.exp(9 * u / 2) - 3 * mp.pi * n ** 2 * mp.exp(5 * u / 2)) * mp.exp(-mp.pi * n * n * x)
        n += 1
    return 2 * s


def dPhi(u):
    # Phi'(u) for u >= 0
    x = mp.exp(2 * u)
    s = mp.mpf(0)
    n = 1
    while True:
        if n > 1 and mp.pi * (n * n - 1) * x > 2.4 * mp.mp.dps + 40:
            break
        A = 2 * mp.pi ** 2 * n ** 4 * mp.exp(9 * u / 2) - 3 * mp.pi * n ** 2 * mp.exp(5 * u / 2)
        dA = 9 * mp.pi ** 2 * n ** 4 * mp.exp(9 * u / 2) - mp.mpf(15) / 2 * mp.pi * n ** 2 * mp.exp(5 * u / 2)
        s += (dA - 2 * mp.pi * n * n * x * A) * mp.exp(-mp.pi * n * n * x)
        n += 1
    return 2 * s


def dh_part():
    mp.mp.dps = 40
    out = {}
    k = kappa_dh()
    chi = {1: mp.mpc(1, 0), 2: mp.mpc(0, 1), 3: mp.mpc(0, -1), 4: mp.mpc(-1, 0)}     # chi mod 5, chi(2) = i, odd
    tau = mp.fsum(chi[a] * mp.exp(2j * mp.pi * a / 5) for a in range(1, 5))
    eps = tau / (1j * mp.sqrt(5))
    c = (1 - 1j * k) / 2
    out['kappa'] = mp.nstr(k, 25)
    out['kappa_as_tan_half_arg_eps'] = mp.nstr(mp.tan(mp.arg(eps) / 2), 25)
    out['eps'] = mp.nstr(eps, 20)
    out['abs_eps'] = mp.nstr(abs(eps), 20)
    out['(1+ik)/(1-ik)'] = mp.nstr((1 + 1j * k) / (1 - 1j * k), 20)
    out['a_n_from_c'] = [mp.nstr(2 * mp.re(c * chi[n]), 15) for n in range(1, 5)]
    cs = case('DH')
    fe = {}
    for s in (mp.mpc('0.3', '2'), mp.mpc('0.8', '15.5')):
        fe[mp.nstr(s, 4)] = [mp.nstr(xi_a(cs, s), 18), mp.nstr(xi_a(cs, 1 - s), 18)]
    out['FE_xi(s)_xi(1-s)'] = fe
    mp.mp.dps = 25
    zs = []
    for z0 in (mp.mpc('0.8085', '85.6993'), mp.mpc('0.6508', '114.1633')):
        z = mp.findroot(lambda s: L(cs, s), z0)
        zs.append({'zero': mp.nstr(z, 15), 'abs_L': mp.nstr(abs(L(cs, z)), 3)})
    out['zeros_off_line'] = zs
    mp.mp.dps = 40
    # the other self dual combination mod 5 with a_1 = 1 (sign -1): theta(1/y) = -y^{3/2} theta(y), so theta(1) = 0
    cm = case('DHminus')
    out['other_combination_coeffs'] = [mp.nstr(cm[1][n], 12) for n in range(1, 5)]
    out['other_combination_FE_sign_check'] = [mp.nstr((theta(cm, 1 / y) + y ** 1.5 * theta(cm, y)) / theta(cm, y), 3)
                                              for y in (mp.mpf('1.3'), mp.mpf(2))]
    out['other_combination_theta(1)'] = mp.nstr(theta(cm, mp.mpf(1)), 3)
    out['L0_DH'] = mp.nstr(mp.fsum(cs[1][r] * (mp.mpf(1) / 2 - mp.mpf(r) / 5) for r in cs[1]), 20)
    out['xi0_DH'] = mp.nstr(mp.sqrt(5) * mp.fsum(cs[1][r] * (mp.mpf(1) / 2 - mp.mpf(r) / 5) for r in cs[1]), 20)
    return out


def kernel_part(name):
    mp.mp.dps = 30
    cs = case(name)
    q, a, p = cs
    out = {'q': q, 'parity': p}
    errs = []
    for y in (mp.mpf('0.45'), mp.mpf('0.8'), mp.mpf('1.7')):
        errs.append(abs(theta(cs, 1 / y) - y ** (p + mp.mpf(1) / 2) * theta(cs, y)) / abs(theta(cs, 1 / y)))
    out['FE_max_rel_err'] = mp.nstr(max(errs), 3)
    S = mp.fsum((n ** p) * abs(a[n % q]) * mp.exp(-mp.pi * (n * n - 1) / q) for n in range(2, 400) if n % q in a)
    out['S'] = mp.nstr(S, 6)
    # Mellin identity: int K(u) e^{(s-1/2)u} du = xi(s, a), and E B^{s/2} = xi(s, a)/xi(0, a) by the density theta(b) b^{p/2-1}
    def K(u):
        return 2 * mp.exp((p + mp.mpf(1) / 2) * abs(u)) * theta(cs, mp.exp(2 * abs(u)))
    rel = []
    for s in (mp.mpf(2), mp.mpf(3), mp.mpc('0.7', '3')):
        I = mp.quad(lambda u: K(u) * mp.exp((s - mp.mpf(1) / 2) * u), [-6, -2, -0.5, 0, 0.5, 2, 6])
        rel.append(abs(I - xi_a(cs, s)) / abs(xi_a(cs, s)))
    out['Mellin_max_rel_err'] = mp.nstr(max(rel), 3)
    return out


def bpy_part():
    mp.mp.dps = 30
    out = {}
    # density of B = pi S_2/2: f(b) = sum (4 pi^2 m^4 b - 6 pi m^2) e^{-pi m^2 b} = b^{-5/4} Phi(log(b)/2)
    def f(b):
        return mp.fsum((4 * mp.pi ** 2 * m ** 4 * b - 6 * mp.pi * m ** 2) * mp.exp(-mp.pi * m * m * b) for m in range(1, 80))
    out['density_vs_kernel_rel_err'] = mp.nstr(max(abs(f(b) - b ** mp.mpf(-1.25) * Phi(mp.log(b) / 2)) / f(b)
                                                   for b in (mp.mpf('0.5'), mp.mpf(1), mp.mpf(2))), 3)
    def fall(b):    # all b > 0 through the kernel
        return b ** mp.mpf(-1.25) * Phi(mp.log(b) / 2)
    tot = mp.quad(fall, [0, 0.1, 0.5, 1, 2, 5, 20, mp.inf])
    out['total_mass'] = mp.nstr(tot, 20)
    lt = []
    for lam in (mp.mpf(1), mp.mpf(3)):
        I = mp.quad(lambda b: mp.exp(-lam * b) * fall(b), [0, 0.1, 0.5, 1, 2, 5, 20, mp.inf])
        lt.append(abs(I - (mp.sqrt(mp.pi * lam) / mp.sinh(mp.sqrt(mp.pi * lam))) ** 2))
    out['laplace_vs_product_max_err'] = mp.nstr(max(lt), 3)
    # E B^{s/2} = 2 xi(s)
    xi = lambda s: s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)
    ms = []
    for s in (mp.mpf(2), mp.mpf(5), mp.mpc('0.5', '2')):
        I = mp.quad(lambda b: b ** (s / 2) * fall(b), [0, 0.1, 0.5, 1, 2, 5, 20, mp.inf])
        ms.append(abs(I - 2 * xi(s)) / abs(2 * xi(s)))
    out['mellin_2xi_max_rel_err'] = mp.nstr(max(ms), 3)
    return out


def R_num(cs, w):
    q, a, p = cs
    return mp.fsum(a[n] * mp.sinh((q - 2 * n) * w / q) for n in range(1, (q + 1) // 2) if n in a and 2 * n < q)


def reflected_part(name):
    mp.mp.dps = 30
    cs = case(name)
    q, a, p = cs
    out = {}
    R0 = mp.fsum(a[n] * mp.mpf(q - 2 * n) / q for n in range(1, (q + 1) // 2) if n in a and 2 * n < q)
    I0 = mp.quad(lambda x: theta_any(cs, x), [0, 0.25, 1, 4, 16, mp.inf])
    out['int_theta'] = mp.nstr(I0, 20)
    out['sqrtq_R0'] = mp.nstr(mp.sqrt(q) * R0, 20)
    out['xi(0,a)'] = mp.nstr(xi_a(cs, mp.mpf(0)), 20)          # odd a: xi(0, a) = sqrt(q) L(0, a)
    errs = []
    for lam in (mp.mpf('0.3'), mp.mpf(1), mp.mpf(3)):
        w = mp.sqrt(mp.pi * q * lam)
        I = mp.quad(lambda x: mp.exp(-lam * x) * theta_any(cs, x), [0, 0.25, 1, 4, 16, mp.inf])
        errs.append(abs(I / I0 - R_num(cs, w) / mp.sinh(w) / R0))
    out['closed_form_max_err'] = mp.nstr(max(errs), 3)
    # the transform on the negative axis: lam = -y^2/(pi q), w = i y; F(y) = num(iy)/(i sin y)/R0
    def num_i(y):
        return mp.im(R_num(cs, mp.mpc(0, y)))
    out['num_at_pi'] = mp.nstr(num_i(mp.pi), 15)
    out['num_at_2pi'] = mp.nstr(num_i(2 * mp.pi), 15)
    zs = []
    ys = [mp.pi * k / 400 for k in range(1, 1600)]
    for y0, y1 in zip(ys[:-1], ys[1:]):
        if num_i(y0) * num_i(y1) < 0:
            z = mp.findroot(num_i, (y0, y1), solver='bisect')
            zs.append({'y_over_pi': mp.nstr(z / mp.pi, 10), 'sin_y': mp.nstr(mp.sin(z), 6),
                       'lam0': mp.nstr(-z ** 2 / (mp.pi * q), 10)})
    out['numerator_sign_changes_y_in_(0,4pi)'] = zs
    return out


def special_part():
    mp.mp.dps = 30
    out = {}
    # chi_{-8}: (sinh(3w/4) + sinh(w/4))/sinh w = cosh(w/4)/cosh(w/2), zero at lam = -pi/2 (w = 2 pi i)
    errs = [abs((mp.sinh(3 * w / 4) + mp.sinh(w / 4)) / mp.sinh(w) - mp.cosh(w / 4) / mp.cosh(w / 2)) for w in (mp.mpf('0.7'), mp.mpc(1, 2), mp.mpf(3))]
    out['chi-8_identity_max_err'] = mp.nstr(max(errs), 3)
    w = mp.sqrt(mp.mpc(8 * mp.pi * (-mp.pi / 2)))
    out['chi-8_transform_at_-pi/2'] = mp.nstr(mp.cosh(w / 4) / mp.cosh(w / 2), 3)
    out['chi-8_product_check'] = mp.nstr(abs(mp.cosh(mp.sqrt(8 * mp.pi * 1) / 4) / mp.cosh(mp.sqrt(8 * mp.pi * 1) / 2)
                                         - mp.nprod(lambda k: (1 + 1 / (4 * mp.pi * (2 * k - 1) ** 2 / 8)) / (1 + 1 / (mp.pi * (2 * k - 1) ** 2 / 8)), [1, mp.inf])), 3)
    # chi_{-8} as pi T_1^{1/2}(R_1): cosh(a sqrt(2 lam'))/cosh(sqrt(2 lam')), lam' = pi lam, a = 1/2  [BPY01, (92)]
    lam = mp.mpf('0.7')
    wv = mp.sqrt(8 * mp.pi * lam)
    out['chi-8_bpy92_err'] = mp.nstr(abs(mp.cosh(wv / 4) / mp.cosh(wv / 2) - mp.cosh(mp.sqrt(2 * mp.pi * lam) / 2) / mp.cosh(mp.sqrt(2 * mp.pi * lam))), 3)
    # chi_{-3}: 3 sinh(w/3)/sinh w = prod_{3 notdivides m} (1 + 3 lam/(pi m^2))^{-1}, and = (3 pi/2) T_1^{1/3}(R_3)  [BPY01, (92)]
    lam = mp.mpf('1.3')
    w = mp.sqrt(3 * mp.pi * lam)
    lhs = 3 * mp.sinh(w / 3) / mp.sinh(w)
    logp = -mp.nsum(lambda m: mp.log(1 + 3 * lam / (mp.pi * m * m)) - mp.log(1 + 3 * lam / (mp.pi * (3 * m) ** 2)), [1, mp.inf])
    out['chi-3_product_err'] = mp.nstr(abs(lhs - mp.exp(logp)), 3)
    lp = 3 * mp.pi * lam / 2
    out['chi-3_bpy92_err'] = mp.nstr(abs(lhs - mp.sinh(mp.sqrt(2 * lp) / 3) / (mp.sinh(mp.sqrt(2 * lp)) / 3)), 3)
    # chi_4: sinh(w/2)/sinh(w) / R(0) = sech(w/2) = sech sqrt(pi lam), w = sqrt(4 pi lam)
    lam = mp.mpf('0.9')
    w = mp.sqrt(4 * mp.pi * lam)
    out['chi-4_sech_err'] = mp.nstr(abs(2 * mp.sinh(w / 2) / mp.sinh(w) - mp.sech(mp.sqrt(mp.pi * lam))), 3)
    # numerators on w = i y: chi_{-7}: sin(5y/7) + sin(3y/7) - sin(y/7) = +-sqrt(7)/2 at y = pi, 2 pi
    out['chi-7_num_pi_minus_sqrt7/2'] = mp.nstr(mp.sin(5 * mp.pi / 7) + mp.sin(3 * mp.pi / 7) - mp.sin(mp.pi / 7) - mp.sqrt(7) / 2, 3)
    out['chi-7_num_2pi_plus_sqrt7/2'] = mp.nstr(mp.sin(10 * mp.pi / 7) + mp.sin(6 * mp.pi / 7) - mp.sin(2 * mp.pi / 7) + mp.sqrt(7) / 2, 3)
    return out


def arcsine_part(t):
    mp.mp.dps = 20
    t = mp.mpf(t)
    Z = 2 * mp.quad(lambda u: mp.exp(t * u * u / 4) * Phi(u), [0, 0.5, 1, 2, 4])

    def p(u):
        return mp.exp(t * u * u / 4) * Phi(u) / Z

    def dp(u):
        return mp.exp(t * u * u / 4) * (t * u / 2 * Phi(u) + dPhi(u)) / Z

    def h(rho):
        # h_t(rho) = -2 rho int_rho^oo p'(r) dr/sqrt(r^2 - rho^2) = -2 rho int_0^oo p'(rho cosh w) dw
        W = mp.acosh(max(mp.mpf(5) / rho, mp.mpf(1)))
        return -2 * rho * mp.quad(lambda w: dp(rho * mp.cosh(w)), [0, W / 4, W / 2, W])
    tot = mp.quad(h, [0, 0.25, 0.5, 1, 1.5, 2, 3])
    rec = []
    for u in (mp.mpf('0.1'), mp.mpf('0.4'), mp.mpf('0.7'), mp.mpf(1)):
        V = mp.acosh(mp.mpf(3) / u)
        r = mp.quad(lambda v: h(u * mp.cosh(v)), [0, V / 4, V / 2, V]) / mp.pi
        rec.append(abs(r - p(u)) / p(u))
    hmin = min(h(mp.mpf(k) / 20) for k in range(1, 41))
    dmax = max(dp(mp.mpf(k) / 20) for k in range(1, 61))
    return {'t': mp.nstr(t, 3), 'int_h': mp.nstr(tot, 15), 'reconstruction_max_rel_err': mp.nstr(max(rec), 3),
            'min_h_on_grid_(0,2]': mp.nstr(hmin, 6), 'max_dp_on_grid_(0,3]': mp.nstr(dmax, 6)}


def ImPsi(x):
    M = int(mp.floor(mp.sqrt(x / mp.pi)))
    s = mp.mpf(0)
    for m in range(1, M + 1):
        sm = x - mp.pi * m * m
        s += 4 * mp.pi ** 2 * m ** 4 * (mp.log(sm) - 1 + mp.euler) / sm ** 2 + 6 * mp.pi * m ** 2 * (mp.log(sm) + mp.euler) / sm
    return mp.pi / 8 * s


def ImPsi_direct(x):
    # psi(lam) = (1/16) int e^{-lam b} f(b) (log b)^2 db at lam = -x + i0: [0, 1] on the real axis (real), and for b >= 1
    # each term with pi m^2 < x along b = 1 + r e^{-3 i pi/4}, the others on [1, oo)
    tot = mp.mpc(0)
    for m in range(1, 40):
        A, B = 4 * mp.pi ** 2 * m ** 4, -6 * mp.pi * m * m
        s = mp.pi * m * m - x
        if s > 0:
            continue
        w = mp.exp(-3j * mp.pi / 4)
        tot += mp.quad(lambda r: (A * (1 + r * w) + B) * mp.exp(-s * (1 + r * w)) * mp.log(1 + r * w) ** 2 * w,
                       [0, 1, 1 + 20 / abs(s), mp.inf])
    return mp.im(tot) / 16


def perturb_part():
    mp.mp.dps = 30
    out = {'check': {}}
    for xs in ('2', '3.3', '6', '12'):
        x = mp.pi * mp.mpf(xs)
        out['check'][xs + 'pi'] = [mp.nstr(ImPsi(x), 15), mp.nstr(ImPsi_direct(x), 15)]

    def phi0(x):
        return mp.pi * x / mp.sin(mp.sqrt(mp.pi * x)) ** 2

    def G(x):
        return ImPsi(x) / phi0(x)
    gaps = {}
    for M in (1, 2, 3, 4):
        lo, hi = mp.pi * M * M, mp.pi * (M + 1) ** 2
        xs = [lo + (hi - lo) * k / 400 for k in range(1, 400)]
        d = [mp.diff(G, x) for x in xs]
        neg = [x for x, v in zip(xs, d) if v < 0]
        gaps['(%d pi, %d pi)' % (M * M, (M + 1) ** 2)] = {
            'Gprime_negative_on_x/pi': [mp.nstr(neg[0] / mp.pi, 6), mp.nstr(neg[-1] / mp.pi, 6)] if neg else None,
            'Gprime_positive_points': sum(1 for v in d if v > 0), 'Gprime_negative_points': len(neg),
            'grid_negative_set_is_interval': (len(neg) == 0 or all(lo + (hi - lo) * k / 400 in neg for k in range(1, 400)
                                                                   if neg[0] <= lo + (hi - lo) * k / 400 <= neg[-1]))}
    out['gaps'] = gaps
    out['Gprime_1.015pi'] = mp.nstr(mp.diff(G, mp.pi * mp.mpf('1.015')), 10)
    out['Gprime_2.47pi'] = mp.nstr(mp.diff(G, mp.pi * mp.mpf('2.47')), 10)
    return out


def task(args):
    kind, x = args
    if kind == 'dh':
        return dh_part()
    if kind == 'kernel':
        return kernel_part(x)
    if kind == 'bpy':
        return bpy_part()
    if kind == 'reflected':
        return reflected_part(x)
    if kind == 'special':
        return special_part()
    if kind == 'arcsine':
        return arcsine_part(x)
    if kind == 'perturb':
        return perturb_part()


if __name__ == '__main__':
    ncpu = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    t0 = time.time()
    tasks = [('dh', None), ('bpy', None), ('special', None), ('perturb', None), ('arcsine', '0'), ('arcsine', '-1')]
    tasks += [('kernel', nm) for nm in NAMED]
    tasks += [('reflected', nm) for nm in ('chi-4', 'chi-3', 'chi-7', 'chi-8', 'DH')]
    with Pool(ncpu) as P:
        res = P.map(task, tasks, chunksize=1)
    out = {}
    for (kind, x), r in zip(tasks, res):
        out[kind if x is None else '%s %s' % (kind, x)] = r
    out['time_s'] = round(time.time() - t0)
    with open('pakes_characters.json', 'w') as f:
        json.dump(out, f, indent=1)
    for k, v in out.items():
        print(k, ':', json.dumps(v))
