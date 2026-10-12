# wigfree_gamma2.py: the Gamma(2) and Laguerre material of the free probability section (misc/reviews/nick-wigner/merge/wigfree.tex,
# label sec:wigfree:gamma2), including the Hankel table tab:wigfree:hankel. Added for the book in October 2026, from the first
# author's note on the semicircle law and the Thorin framework (misc/reviews/nick-wigner/source-v7/, version 7, Section 10),
# recomputed independently.
# 1. Williams' variable S_2 = sum_n 2 G_n/(pi^2 n^2), G_n Gamma(2,1): its moments from the cumulants 2 (m-1)! (2/pi^2)^m zeta(2m),
#    the identity 2 xi(s) = E (pi S_2/2)^{s/2} at s = 2, 4, 6 and at s = 1, through E X^{1/2} = (2 sqrt pi)^{-1}
#    int_0^oo lambda^{-3/2} (1 - E e^{-lambda X}) d lambda; and the one dimensional model S_1 = sum_n 2 e_n/(pi^2 n^2), whose
#    Mellin transform is E (pi S_1/2)^{s/2} = 2 (1 - 2^{1-s}) xi(s)/(s - 1), checked at s = 1, 2, 4, 6.
# 2. The Laguerre dictionary: the coefficients of the reversed Laguerre polynomial x^d L_d^{(a-1)}(1/x) in the Marcus Spielman
#    Srivastava normalisation are Gamma(d + a)/(d! Gamma(k + a)) (sympy, d <= 8, four values of a); the Bessel form of
#    sum_k x^k/(k! (a)_k) and its zeros -j_{a-1,m}^2/4; the example (1 + s)^2.
# 3. The Hankel test of Theorem thm:wigfree:nomix: for X = |U| and X = U^2, U with the Polya density Phi/Xi(0), and for X = S_2,
#    mu_k = E X^k/(a)_k must be a Stieltjes moment sequence if X is a scale mixture of Gamma(a) laws; the script records, for each
#    a, the least order N at which [mu_{i+j}]_{i,j<N} and the least at which [mu_{i+j+1}]_{i,j<N} is not positive definite
#    (Cholesky, 400 digits, repeated at 250 digits); the Stieltjes test fails at the smaller of the two. The moments of |U|
#    come from a tanh sinh rule (step 1/320) on [0, 3.6] (odd powers of |u| are not smooth at 0, so the trapezoid rule of
#    jacobi_zeros.py would lose accuracy), checked against the trapezoid values for even powers.
# 4. The 2 x 2 model W = sum_n c_n z_n z_n^*, c_n = 2/(pi^2 n^2), z_n standard complex Gaussian in C^2: the mixed characteristic
#    polynomial x^2 - 2 m_1 x + m_1^2 - m_2 (sympy), m_1 = 1/3, m_2 = 2/45, and a Monte Carlo check of E tr W and E det W.
# Run from any folder: python3 wigfree_gamma2.py; writes wigfree_gamma2.json there (a few minutes on one core).
import json
import time

import mpmath as mp
import numpy as np
import sympy as sp

t0 = time.time()
out = {}


def rec(k, v):
    out[k] = v
    print(k, v, round(time.time() - t0, 1), flush=True)


def s_(x, n=12):
    return mp.nstr(x, n)


# ---------------- 1. Williams' variable and the one dimensional model ----------------
mp.mp.dps = 40


def xi(s):
    # the completed zeta function, with its value 1/2 at s = 1, where (s - 1) zeta(s) -> 1
    if s == 1:
        return mp.mpf(1) / 2
    return s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)


def moments_from_cumulants(kap, K):
    # m_n = sum_{j=1}^n binom(n-1, j-1) kappa_j m_{n-j}
    m = [mp.mpf(1)]
    for n in range(1, K + 1):
        m.append(mp.fsum(mp.binomial(n - 1, j - 1) * kap[j] * m[n - j] for j in range(1, n + 1)))
    return m


K2 = [None] + [2 * mp.factorial(j - 1) * (2 / mp.pi ** 2) ** j * mp.zeta(2 * j) for j in range(1, 12)]
K1 = [None] + [mp.factorial(j - 1) * (2 / mp.pi ** 2) ** j * mp.zeta(2 * j) for j in range(1, 12)]
mS2 = moments_from_cumulants(K2, 6)
mS1 = moments_from_cumulants(K1, 6)
rec('E S_2, E S_2^2 (2/3, 8/15)', [s_(mS2[1], 20), s_(mS2[2], 20)])
rec('E (pi S_2/2)^{s/2} - 2 xi(s) at s = 2, 4, 6', [s_(mS2[1] * mp.pi / 2 - 2 * xi(mp.mpf(2)), 3), s_(mS2[2] * (mp.pi / 2) ** 2 - 2 * xi(mp.mpf(4)), 3),
                                                     s_(mS2[3] * (mp.pi / 2) ** 3 - 2 * xi(mp.mpf(6)), 3)])
# fractional moments through E X^{1/2} = (2 sqrt(pi))^{-1} int_0^oo lambda^{-3/2} (1 - L(lambda)) d lambda
L1 = lambda lam: mp.sqrt(2 * lam) / mp.sinh(mp.sqrt(2 * lam)) if lam > 0 else mp.mpf(1)
half = lambda L: mp.quad(lambda lam: lam ** mp.mpf(-1.5) * (1 - L(lam)), [0, 1, 10, 100, mp.inf]) / (2 * mp.sqrt(mp.pi))
ES2h = half(lambda lam: L1(lam) ** 2)
ES1h = half(L1)
rec('E (pi S_2/2)^{1/2} against 2 xi(1) = 1', [s_(mp.sqrt(mp.pi / 2) * ES2h, 20)])


def mellin_S1(s):
    # 2 (1 - 2^{1-s}) xi(s)/(s - 1), the limit 2 log 2 xi(1) at s = 1
    if s == 1:
        return 2 * mp.log(2) * xi(mp.mpf(1))
    return 2 * (1 - mp.mpf(2) ** (1 - s)) * xi(s) / (s - 1)


rec('E (pi S_1/2)^{s/2} against 2 (1 - 2^{1-s}) xi(s)/(s - 1) at s = 1, 2, 4, 6 (differences)',
    [s_(mp.sqrt(mp.pi / 2) * ES1h - mellin_S1(mp.mpf(1)), 3), s_(mS1[1] * mp.pi / 2 - mellin_S1(mp.mpf(2)), 3),
     s_(mS1[2] * (mp.pi / 2) ** 2 - mellin_S1(mp.mpf(4)), 3), s_(mS1[3] * (mp.pi / 2) ** 3 - mellin_S1(mp.mpf(6)), 3)])
rec('E S_1, E S_1^2 (1/3, 7/45), E (pi S_1/2)^{1/2} = 2 log 2 xi(1) = log 2', [s_(mS1[1], 15), s_(mS1[2], 15), s_(mp.sqrt(mp.pi / 2) * ES1h, 15)])
rec('zeros of 1 - 2^{1-s}: s = 1 + 2 pi i k/log 2, first ordinate', s_(2 * mp.pi / mp.log(2), 12))

# ---------------- 2. The Laguerre dictionary ----------------
a_, x_ = sp.symbols('a x', positive=True)
ok = True
for av in (sp.Rational(1, 2), 2, sp.Rational(7, 3), 6):
    for d in range(1, 9):
        q = sp.Poly(sp.expand(x_ ** d * sp.assoc_laguerre(d, av - 1, 1 / x_)), x_)
        for k in range(d + 1):
            # coefficient of x^{d-k} is (-1)^k binom(d,k) e_k with e_k = Gamma(d+a)/(d! Gamma(k+a))
            ek = sp.gamma(d + av) / (sp.factorial(d) * sp.gamma(k + av))
            if sp.simplify(q.coeff_monomial(x_ ** (d - k)) - (-1) ** k * sp.binomial(d, k) * ek) != 0:
                ok = False
rec('coefficients of x^d L_d^{(a-1)}(1/x) are (-1)^k binom(d,k) Gamma(d+a)/(d! Gamma(k+a)), d <= 8, a = 1/2, 2, 7/3, 6', ok)
# the reversed Laguerre polynomial has positive roots; an example of p boxtimes_d q_a for p with negative roots
d = 8
for av in (sp.Rational(1, 2), 2, 6):
    q = sp.Poly(sp.expand(x_ ** d * sp.assoc_laguerre(d, av - 1, 1 / x_)), x_)
    rts = [complex(r) for r in sp.Poly(q, x_).nroots(n=30)]
    rec('roots of x^%d L_%d^{(a-1)}(1/x), a = %s: least real part, largest |imaginary part|' % (d, d, av),
        ['%.6g' % min(r.real for r in rts), '%.2g' % max(abs(r.imag) for r in rts)])
mp.mp.dps = 30
for av in (mp.mpf('0.5'), mp.mpf(1), mp.mpf(2), mp.mpf('6.5')):
    xv = mp.mpf('3.7')
    lhs = mp.nsum(lambda k: xv ** k / (mp.factorial(k) * mp.rf(av, k)), [0, mp.inf])
    rhs = mp.gamma(av) * xv ** ((1 - av) / 2) * mp.besseli(av - 1, 2 * mp.sqrt(xv))
    j1 = mp.pi / 2 if av == mp.mpf('0.5') else mp.besseljzero(av - 1, 1)
    zv = mp.nsum(lambda k: (-j1 ** 2 / 4) ** k / (mp.factorial(k) * mp.rf(av, k)), [0, mp.inf])
    rec('a = %s: sum x^k/(k!(a)_k) minus Gamma(a) x^{(1-a)/2} I_{a-1}(2 sqrt x) at x = 3.7; value at -j_{a-1,1}^2/4; that point' % s_(av, 3),
        [s_(lhs - rhs, 3), s_(zv, 3), s_(-j1 ** 2 / 4, 10)])
s_sym = sp.symbols('s')
rec('(1+s)^2 inflated by (2)_k: 1 + 4 s + 6 s^2, discriminant', str(sp.discriminant(1 + 4 * s_sym + 6 * s_sym ** 2, s_sym)))

# ---------------- 3. The Hankel test ----------------


def phi_val(v):
    q = mp.exp(-mp.pi * mp.exp(2 * v))
    a9 = mp.exp(9 * v / 2)
    a5 = mp.exp(5 * v / 2)
    tot = mp.mpf(0)
    qn = q
    step = q ** 3
    n = 1
    while True:
        term = (2 * mp.pi ** 2 * n ** 4 * a9 - 3 * mp.pi * n ** 2 * a5) * qn
        tot += term
        if n > 3 and abs(term) < abs(tot) * mp.eps * mp.mpf(10) ** -20:
            break
        qn *= step
        step *= q * q
        n += 1
    return 2 * tot


def polya_moments(dps, K):
    # E|U|^k, k = 0..K, by a tanh sinh rule on [0, V]; and E U^{2k} by the trapezoid rule (step 1/600) for comparison
    mp.mp.dps = dps
    V = mp.mpf('3.6')
    h = mp.mpf(1) / 320
    tmax = mp.mpf('6.6')
    J = int(tmax / h)
    acc = [mp.mpf(0)] * (K + 1)
    for j in range(-J, J + 1):
        t = j * h
        sh = mp.pi / 2 * mp.sinh(t)
        x = V / 2 * (1 + mp.tanh(sh))
        if x <= 0 or x >= V:
            continue
        w = V / 2 * h * mp.pi / 2 * mp.cosh(t) / mp.cosh(sh) ** 2
        f = w * phi_val(x)
        p = mp.mpf(1)
        for k in range(K + 1):
            acc[k] += f * p
            p *= x
    # trapezoid even moments E U^{2k}, k = 0..K (the moments of X = U^2 up to order K)
    ht = mp.mpf(1) / 600
    ev = [mp.mpf(0)] * (K + 1)
    for j in range(int(V / ht) + 1):
        v = j * ht
        f = 2 * ht * phi_val(v) * (mp.mpf(1) / 2 if j == 0 else 1)
        p = mp.mpf(1)
        for k in range(K + 1):
            ev[k] += f * p
            p *= v * v
    M0 = 2 * acc[0]
    absm = [2 * acc[k] / M0 for k in range(K + 1)]
    sq = [ev[k] / ev[0] for k in range(K + 1)]
    diff = max(abs(absm[2 * k] - sq[k]) / sq[k] for k in range(K // 2 + 1))
    return absm, sq, diff, ev[0]


def first_fail(mom, a, Nmax):
    # the least order N at which [mu_{i+j}] (shift 0, the Hamburger condition) and [mu_{i+j+1}] (shift 1) fail to be positive
    # definite; a Stieltjes moment sequence needs both, so the Stieltjes test fails at the smaller of the two
    mu = [mom[k] / mp.rf(a, k) for k in range(len(mom))]
    res = []
    for shift in (0, 1):
        first = None
        for N in range(1, Nmax + 1):
            Hm = mp.matrix(N, N)
            for i in range(N):
                for j in range(N):
                    Hm[i, j] = mu[i + j + shift]
            try:
                mp.cholesky(Hm)
            except (ValueError, ZeroDivisionError):
                first = N
                break
        res.append(first)
    return res


AS = ['0.5', '1', '1.5', '2', '2.5', '3', '4', '6', '8', '12', '20']
NMAX = 20
res = {}
for dps in (400, 250):
    absm, sq, diff, M0 = polya_moments(dps, 2 * NMAX)
    if dps == 400:
        rec('tanh sinh even moments of |U| against the trapezoid values, largest relative difference', s_(diff, 3))
        rec('M_0 = Xi(0) against xi(1/2)', s_(M0 - xi(mp.mpf(1) / 2), 3))
        rec('E|U|, E U^2, E|U|^3', [s_(absm[1], 15), s_(absm[2], 15), s_(absm[3], 15)])
    mS2h = moments_from_cumulants([None] + [2 * mp.factorial(j - 1) * (2 / mp.pi ** 2) ** j * mp.zeta(2 * j) for j in range(1, 2 * NMAX + 1)], 2 * NMAX)
    for name, mom in (('|U|', absm), ('U^2', sq), ('S_2', mS2h)):
        for av in AS:
            res.setdefault(name, {}).setdefault(av, []).append(first_fail(mom, mp.mpf(av), NMAX))
rec('first failing order of [mu_{i+j}] and of [mu_{i+j+1}] at 400 and 250 digits, by variable and a (None: positive definite to order 20)', res)
table = {name: {av: min(v for v in res[name][av][0] if v is not None) for av in AS} for name in res}
rec('first failing order of the Stieltjes test (the smaller of the two), 400 digits', table)
table0 = {name: {av: res[name][av][0][0] for av in AS} for name in res}
rec('first failing order of the Hamburger test [mu_{i+j}] alone, 400 digits', table0)
rec('agreement of 400 and 250 digits', all(res[n][av][0] == res[n][av][1] for n in res for av in AS))

# ---------------- 4. The 2 x 2 model ----------------
cs = sp.symbols('c1:5')
zs = sp.symbols('z1:5')
X = sp.symbols('X')
f = (X + sum(c * z for c, z in zip(cs, zs))) ** 2
for z in zs:
    f = sp.expand(f - sp.diff(f, z))
mcp = sp.expand(f.subs({z: 0 for z in zs}))
m1 = sum(cs)
m2 = sum(c ** 2 for c in cs)
rec('mixed characteristic polynomial of c_i I_2 minus x^2 - 2 m_1 x + m_1^2 - m_2 (four terms, symbolic)', str(sp.simplify(mcp - (X ** 2 - 2 * m1 * X + m1 ** 2 - m2))))
mp.mp.dps = 30
m1v = 2 * mp.zeta(2) / mp.pi ** 2
m2v = 4 * mp.zeta(4) / mp.pi ** 4
rec('m_1, m_2, roots m_1 -+ sqrt(m_2)', [s_(m1v, 15), s_(m2v, 15), s_(m1v - mp.sqrt(m2v), 12), s_(m1v + mp.sqrt(m2v), 12)])
rng = np.random.default_rng(20261009)
Nt = 300
c = 2 / (np.pi ** 2 * np.arange(1, Nt + 1) ** 2)
trs, dets = [], []
for batch in range(40):
    Z = (rng.standard_normal((5000, Nt, 2)) + 1j * rng.standard_normal((5000, Nt, 2))) / np.sqrt(2)
    W = np.einsum('n,bni,bnj->bij', c, Z, Z.conj())
    trs.append(np.trace(W, axis1=1, axis2=2).real)
    dets.append(np.linalg.det(W).real)
trs = np.concatenate(trs)
dets = np.concatenate(dets)
m1t, m2t = c.sum(), (c ** 2).sum()
rec('Monte Carlo, 2 x 10^5 draws, n <= 300: mean tr W, its standard error, exact 2 m_1(300)', ['%.6f' % trs.mean(), '%.1g' % (trs.std() / np.sqrt(len(trs))), '%.6f' % (2 * m1t)])
rec('Monte Carlo: mean det W, its standard error, exact m_1(300)^2 - m_2(300)', ['%.6f' % dets.mean(), '%.1g' % (dets.std() / np.sqrt(len(dets))), '%.6f' % (m1t ** 2 - m2t)])
out['seconds'] = round(time.time() - t0, 1)
json.dump(out, open('wigfree_gamma2.json', 'w'), indent=1)
print('done', out['seconds'])
