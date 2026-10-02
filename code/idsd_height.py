# idsd_height.py: how low the verification height in Theorem thm:ch12:idsd can go, the remark after its proof in Section
# ch:heat (ch/ch12.tex). Added for the book in September 2026, from the checks of the authors' note on the theorem
# (misc/reviews/nick-thm511/, item 141 of the authors' report). As in the proof, c_k = (-1)^k W^(k) > 0 for k = 0, 1 follows
# on (0, tau] from the arithmetic side and on [tau, oo) from the zeros below a height R, which lie on the line. Three versions:
#  (a) the bounds of the proof, with tau = 1/8300 and R = 400;
#  (b) the archimedean error computed instead of bounded: with dev(r) = m(r) - log(r/2pi), I_k = int_0^5 r^{2k} |dev| dr and
#      r^2 |dev(r)| < C = 0.05 for r >= 5, it is at most I_k/(2 pi s_k(t)) + C/25; tau = 1/500 and R = 100;
#  (c) as (b), keeping also the polar term e^{t/4} in c_0, taking the prime terms with n <= 30 exactly, and taking the lower
#      bound for t >= tau from the first zero instead of a zero in (20, 50]; tau = 1/50 and R = 25.
# Every inequality is evaluated in ball arithmetic (python-flint, Arb) and counts only if it holds for the whole ball. Outside
# the certificate the script reports the least integer height at which the bound for t >= tau of each version closes, and the
# integrals I_k by mpmath quadrature. Run from any folder: python3 idsd_height.py; writes idsd_height.json there.
import json
import time

import mpmath as mp
from flint import arb, acb, ctx

ctx.prec = 128
pi, gE, log2 = arb.pi(), arb.const_euler(), arb(2).log()
half = arb(1) / 2
out, holds = {}, []
T0 = time.time()


def s(x, n=12):
    return x.str(n, radius=False) if isinstance(x, arb) else x


def check(name, cond, **vals):
    holds.append(bool(cond))
    out[name] = dict(holds=bool(cond), **{k: s(v) for k, v in vals.items()})
    print(name, bool(cond), round(time.time() - T0, 1), flush=True)


def ball(a, b):  # a ball containing [a, b], for floats a <= b
    return arb(a).union(arb(b))


def m(r):  # m(r) = Re psi(1/4 + i r/2) - log pi, for a real ball r
    return acb(arb(1) / 4, r / 2).digamma().real - pi.log()


def dev(r):
    return m(r) - (r / (2 * pi)).log()


def ddev(r):  # dev'(r) = -Im psi'(1/4 + i r/2)/2 - 1/r, with psi'(z) = zeta(2, z)
    return -acb(2).zeta(acb(arb(1) / 4, r / 2)).imag / 2 - 1 / r


def sup_abs(f, df, a, b):  # an upper bound for |f| on [a, b]: the smaller of the direct enclosure and the mean value form
    X, c, h = ball(a, b), (arb(a) + arb(b)) / 2, (arb(b) - arb(a)) / 2
    u, v = abs(f(X)).upper(), (abs(f(c)) + abs(df(X)) * h).upper()
    return u if u < v else v


# 1. The constants of the archimedean term. I_k is bounded above on [0, a0] by |dev| <= max|m| + log(2 pi/r), and on [a0, 5]
# by upper sums over pieces, geometric to 0.05 and of width 1e-4 beyond. C: r^2 |dev| < 0.05 on [5, 10] by pieces of width
# 0.001, bisected where needed, and for r >= 10 from Binet's formula
#   psi(z) = log z - 1/(2z) - 2 int_0^oo v dv/((v^2+z^2)(e^{2 pi v}-1))
# at z = 1/4 + ir/2, with 1/(v^2+z^2) = 1/z^2 - v^2/(z^2(v^2+z^2)): r^2 dev = T1 - T2 + T3 + T4, where
# T1 = (r^2/2) log(1 + 1/4r^2) is in [1/8 - 1/64r^2, 1/8], T2 = 2r^2/(1 + 4r^2) is in [1/2 - 1/8r^2, 1/2],
# T3 = -r^2 Re(1/12z^2) is in [1/3 - 1/4r^2 - 1/48r^4, 1/3], and |T4| = |2r^2 Re(K/z^2)| <= 8|K| with
# K = int_0^oo v^3 dv/((v^2+z^2)(e^{2 pi v}-1)). Since |v^2+z^2| >= (3r^2-1)/16 for v <= r/4 and >= r/4 beyond,
# |K| <= 1/(15(3r^2-1)) + (4/r) Tail(r/4), where Tail(a) = int_a^oo v^3 dv/(e^{2 pi v}-1) is at most
# e^{-ca}/(1-e^{-ca}) (a^3/c + 3a^2/c^2 + 6a/c^3 + 6/c^4), c = 2 pi. Each term of the bound decreases in r.
a0 = 2.0 ** -20
first = abs(m(ball(0.0, a0))).upper() * a0 + a0 * (1 + (2 * pi / a0).log())
I = [first, first * a0 ** 2]
pts = [a0]
while pts[-1] < 0.05:
    pts.append(pts[-1] * 1.001)
while pts[-1] < 5:
    pts.append(min(5.0, round(pts[-1] + 1e-4, 10)))
for a, b in zip(pts[:-1], pts[1:]):
    u = sup_abs(dev, ddev, a, b) * (arb(b) - arb(a))
    I[0] += u
    I[1] += u * arb(b) ** 2
I_STATED = [arb('0.836'), arb('0.404')]
check('archimedean_integrals', I[0] < I_STATED[0] and I[1] < I_STATED[1], I0_upper=I[0], I1_upper=I[1],
      I_stated='0.836, 0.404', pieces=len(pts) - 1)

def sup_cover(f, df, a, b, bound, depth=0):  # sup_abs on [a, b], bisected while it is not below bound; returns it and the pieces
    u = sup_abs(f, df, a, b)
    if u < bound or depth == 12:
        return u, 1
    c = (a + b) / 2
    (u1, n1), (u2, n2) = sup_cover(f, df, a, c, bound, depth + 1), sup_cover(f, df, c, b, bound, depth + 1)
    return (u1 if u1 > u2 else u2), n1 + n2


C = arb('0.05')
g = lambda r: r * r * dev(r)
dg = lambda r: 2 * r * dev(r) + r * r * ddev(r)
grid = [5 + i / 1000 for i in range(5001)]
supC, nC = arb(0), 0
for a, b in zip(grid[:-1], grid[1:]):
    u, n = sup_cover(g, dg, a, b, C)
    nC += n
    if u > supC:
        supC = u


def binet_tail(r):  # the bound for r^2 |dev(r)| when r >= 10, decreasing in r
    r = arb(r); c = 2 * pi; a = r / 4
    Ta = (-c * a).exp() / (1 - (-c * a).exp()) * (a ** 3 / c + 3 * a ** 2 / c ** 2 + 6 * a / c ** 3 + 6 / c ** 4)
    K = 1 / (15 * (3 * r * r - 1)) + 4 / r * Ta
    return arb(1) / 24 + 1 / (64 * r * r) + 1 / (4 * r * r) + 1 / (48 * r ** 4) + 8 * K


check('archimedean_sup', supC < C and binet_tail(10) < C, sup_on_5_10_upper=supC, pieces=nC, bound_beyond_10=binet_tail(10),
      r2dev_at_5=25 * dev(arb(5)), r2dev_at_10=100 * dev(arb(10)))

# 2. The arithmetic range. c_k(t) >= s_k(t) F_k(t) with F_k = ell + alpha_k - E_k - |P^(k)|/s_k + (-1/4)^k e^{t/4}/s_k, where
# ell = -log(t)/2, s_k(t) = Gamma(k+1/2) t^{-k-1/2}/(4 pi), alpha_k = psi(k+1/2)/2 - log 2pi, E_k bounds the archimedean error
# and P is the prime term.
alpha = [acb(k + half).digamma().real / 2 - (2 * pi).log() for k in (0, 1)]
kappa = log2 ** 2 / 8
LAM = {}  # the von Mangoldt function up to 30
for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]:
    n = p
    while n <= 30:
        LAM[n] = arb(p).log()
        n *= p


def s_k(k, t):
    return (k + half).gamma() * t ** (-(k + half)) / (4 * pi)


def E_book(k, t):  # the proof's bound, from int_0^20 |dev| < 136 and |dev| < 0.003 for r >= 20
    return 272 * arb(20) ** (2 * k) * t ** (k + half) / (k + half).gamma() + arb('0.003')


def E_comp(k, t):  # the computed bound, from I_k and C
    return I_STATED[k] / (2 * pi * s_k(k, t)) + C / 25


def P_book(k, t):  # the proof's bound t^{-k-1/2} e^{-kappa/t} for |P^(k)|, over s_k; valid for q = log 2/(8t) > 3
    return 4 * pi / (k + half).gamma() * (-kappa / t).exp()


def P_sharp(k, t):  # |P^(k)(t)| = (4 pi)^{-1/2} t^{-k-1/2} |sum_n Lambda(n) n^{-1/2} p_k(y_n) e^{-y_n}|, y_n = (log n)^2/4t,
    # p_0 = 1, p_1 = y - 1/2, over s_k; the terms n <= 30 exactly, and for n > 30 Lambda(n) n^{-1/2} < 1, |p_k(y)| e^{-y} <=
    # e^{-y/2} <= n^{-q1} with q1 = log 30/(8t), and sum_{n > 30} n^{-q1} <= 30^{1-q1}/(q1 - 1)
    tot = arb(0)
    for n, L in LAM.items():
        y = arb(n).log() ** 2 / (4 * t)
        tot += L / arb(n).sqrt() * abs(y - half if k else arb(1)) * (-y).exp()
    q1 = arb(30).log() / (8 * t)
    tot += arb(30) ** (1 - q1) / (q1 - 1)
    return 2 * pi.sqrt() / (k + half).gamma() * tot


def F(k, t, E, P, polar0):
    f = -t.log() / 2 + alpha[k] - E(k, t) - P(k, t)
    if k == 1:
        f -= (t / 4).exp() / 4 / s_k(1, t)
    elif polar0:
        f += (t / 4).exp() / s_k(0, t)
    return f


check('prime_constants', 2 * (arb(-3) / 4).exp() < 1 and 2 * (arb(-5) / 4).exp() < 1,
      note='(y + 1/2) e^{-y/2} <= 2 e^{-3/4} < 1, so |p_k(y)| e^{-y} <= e^{-y/2} for k = 0, 1')

# (a) and (b): ell decreases in t, and E_k, P_book and the polar term of c_1 increase, all over s_k, so F_k(t) >= F_k(tau) on
# (0, tau] with the polar term of c_0, which is positive, left out.
for name, tau, E in [('a_arithmetic', arb(1) / 8300, E_book), ('b_arithmetic', arb(1) / 500, E_comp)]:
    f = [F(k, tau, E, P_book, False) for k in (0, 1)]
    q = log2 / (8 * tau)
    check(name, q > 3 and f[0] > 0 and f[1] > 0, F0_at_tau=f[0], F1_at_tau=f[1], E0=E(0, tau), E1=E(1, tau), q=q)

# (c): on (0, t1] the argument of (b); on [t1, 0.02] F_k with the polar term of c_0 and P_sharp, on balls of relative width
# 0.0005 covering the interval, which contains (t1, 1/50].
t1f = 1e-3
t1 = arb(t1f)
f_low = [F(k, t1, E_comp, P_book, False) for k in (0, 1)]
a, nb, ok, low = t1f, 0, True, [None, None]
while a < 0.02:
    b = min(a * 1.0005, 0.02)
    T = ball(a, b)
    for k in (0, 1):
        f = F(k, T, E_comp, P_sharp, True)
        ok = ok and bool(f > 0)
        if low[k] is None or f.lower() < low[k][0]:
            low[k] = (f.lower(), a)
    nb += 1
    a = b
q1_end = arb(30).log() / (8 * arb(0.02))  # q1 of P_sharp decreases in t, so this is its least value on the cover
check('c_arithmetic', log2 / (8 * t1) > 3 and q1_end > 1 and f_low[0] > 0 and f_low[1] > 0 and ok, F0_at_t1=f_low[0],
      F1_at_t1=f_low[1], q1_at_end=q1_end, balls=nb, least_F0_lower=low[0][0], at_t=low[0][1], least_F1_lower=low[1][0],
      at_t_1=low[1][1])

# 3. The spectral range, t >= tau. The zeros below R lie on the line and contribute gamma^{2k} e^{-gamma^2 t} > 0 to c_k. A zero
# above R has Re a >= gamma^2 - 1/4 and |a| <= 2 gamma^2, and with N(r) <= r^2 for r >= R integration by parts bounds the
# modulus of their contributions by H_k(t) = 2^k e^{t/4} Gamma(k+2, t R^2)/t^{k+1}, which needs t R^2 >= k. Against a lower
# bound L_k(t) = g^{2k} e^{-g^2 t}, H_k/L_k is a sum of terms e^{-(R^2 - g^2 - 1/4) t} t^{-j}, so it decreases in t when
# R^2 > g^2 + 1/4 and H_k(tau) < L_k(tau) suffices.
def H(k, t, R):
    x = t * R * R
    G = [(-x).exp() * (1 + x), (-x).exp() * (2 + 2 * x + x * x)][k]
    return arb(2) ** k * (t / 4).exp() * G / t ** (k + 1)


def N_bounds(T):  # Trudgian's bound: |N(T) - (T/2pi) log(T/(2 pi e)) - 7/8| <= 0.112 log T + 0.278 log log T + 2.510 + 0.2/T
    T = arb(T)
    main = T / (2 * pi) * (T / (2 * pi)).log() - T / (2 * pi) + arb(7) / 8
    err = arb('0.112') * T.log() + arb('0.278') * T.log().log() + arb('2.510') + arb('0.2') / T
    return main - err, main + err


e1 = arb(1).exp()
ratio25 = 1 / (2 * pi * e1) + arb('0.39') / (e1 * 25) + arb('3.393') / 625
check('zero_counts', N_bounds(20)[1] < 5 and N_bounds(50)[0] > 6 and ratio25 < 1, N20_upper=N_bounds(20)[1],
      N50_lower=N_bounds(50)[0], N_over_r2_at_most=ratio25,
      note='a zero in (20, 50]; for r >= 25 the bound is at most (r/2pi + 0.39) log r + 3.393 <= r^2/(2 pi e) + 0.39 r/e + 3.393, '
           'so N(r)/r^2 is at most N_over_r2_at_most, which decreases in r')


def xi(s):  # xi(s) = s(s-1)/2 pi^{-s/2} Gamma(s/2) zeta(s), for a complex ball s
    return s * (s - 1) / 2 * (-(s / 2) * pi.log()).exp() * (s / 2).gamma() * s.zeta()


x1, x2 = xi(acb(half, arb('14.1347'))).real, xi(acb(half, arb('14.1348'))).real  # balls around the decimals
check('first_zero', x1 > 0 and x2 < 0, Xi_at_14_1347=x1, Xi_at_14_1348=x2,
      note='Xi(t) = xi(1/2 + it) is real, so it has a zero in between, on the line')
G1 = arb('14.1347').union(arb('14.1348'))


def least_height(tau, L, lo):  # outside the certificate: the least integer R >= lo at which H_k(tau, R) < L_k for k = 0, 1
    R = lo
    while not (tau * R * R > 1 and H(0, tau, arb(R)) < L[0] and H(1, tau, arb(R)) < L[1]):
        R += 1
    return R


least = {}
for name, tau, R in [('a', arb(1) / 8300, 400), ('b', arb(1) / 500, 100)]:  # a zero in (20, 50] contributes >= 400^k e^{-2500t}
    L = [(-2500 * tau).exp(), 400 * (-2500 * tau).exp()]
    h = [H(0, tau, arb(R)), H(1, tau, arb(R))]
    check(name + '_spectral', R >= 50 and arb(R * R) - half / 2 > 2500 and tau * R * R > 1 and h[0] < L[0] and h[1] < L[1],
          R=R, H0=h[0], L0=L[0], H1=h[1], L1=L[1])
    least[name] = least_height(tau, L, 51)
tau, R = arb(1) / 50, 25  # the first zero, in G1, contributes >= min over G1 of g^{2k} e^{-g^2 t}
L = [(-G1 * G1 * tau).exp(), G1 * G1 * (-G1 * G1 * tau).exp()]
h = [H(0, tau, arb(R)), H(1, tau, arb(R))]
check('c_spectral', arb(R * R) > G1 * G1 + half / 2 and tau * R * R > 1 and h[0] < L[0] and h[1] < L[1],
      R=R, H0=h[0], L0=L[0], H1=h[1], L1=L[1])
least['c'] = least_height(tau, L, 15)
out['least_heights_not_certified'] = least

# 4. Cross check, not part of the certificate: I_0 and I_1 by mpmath quadrature
mp.mp.dps = 20
dm = lambda r: mp.re(mp.digamma(mp.mpf(1) / 4 + 0.5j * r)) - mp.log(r / 2)
r0 = mp.findroot(dm, 0.03)
out['crosscheck'] = dict(sign_change_of_dev=mp.nstr(r0, 10),
                         I0=mp.nstr(mp.quad(lambda r: abs(dm(r)), [0, r0, 1, 5]), 12),
                         I1=mp.nstr(mp.quad(lambda r: r * r * abs(dm(r)), [0, r0, 1, 5]), 12),
                         note='mpmath quadrature, not part of the certificate')
out['all_hold'] = all(holds)
out['seconds'] = round(time.time() - T0, 1)
json.dump(out, open('idsd_height.json', 'w'), indent=1)
print('all inequalities hold:', all(holds), '| checks:', len(holds))
print(json.dumps(out, indent=1))
