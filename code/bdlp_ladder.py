# bdlp_ladder.py: the ladder of driving processes of the centre clock, Theorem pk:thm:ladder in Section pk:sec:ladder
# (ch/pakes2.tex). Added for the book in October 2026, for the merge of the first author's chapter on Pakes, the Lambert
# law and the Thorin condition (misc/reviews/nick-pakes/, triage-C.md; this script replaces the interval certificate of the
# chapter). With W^[j] = (-t d/dt)^j W, beta_0 = 1 and beta_{j+1}(u) = u beta_j(u) - u beta_j'(u), every zero contributes
# Re[beta_j(a t) e^{-a t}] to W^[j](t), and (-t d/dt)^j applied to the explicit formula eq:ch2:Warith gives
#   W^[j](t) = beta_j(-t/4) e^{t/4} + A_j(t) - P_j(t),  A_j(t) = (1/2pi) int_0^oo beta_j(r^2 t) e^{-r^2 t} m(r) dr,
#   P_j(t) = (2 sqrt(pi t))^{-1} sum_n Lambda(n) n^{-1/2} e^{-v_n} q_j(v_n),  v_n = (log n)^2/4t,
# with q_0 = 1, q_{j+1}(v) = (1/2 - v) q_j(v) + v q_j'(v). Replacing m(r) by log(r/2pi) in A_j gives the main term
# M_j(t) = 2^{-j} s_0(t) (ell + alpha_0 + j), written \bar A_j in the text, in the notation of the proof of Theorem
# thm:ch12:idsd, and |A_j - M_j| <= E_j(t) =
# (1/2pi)(L_j I_1 t + C sqrt(t) J_j/2) with L_j = sup |beta_j(u)/u| e^{-u}, J_j = int_0^oo |beta_j(u)| e^{-u} u^{-3/2} du,
# I_1 = int_0^5 r^2 |m(r) - log(r/2pi)| dr and r^2 |m(r) - log(r/2pi)| <= C for r >= 5. Certified here:
#  1. the polynomials: the Touchard and Stirling forms of beta_j, the Mellin identity int beta_j(u) e^{-u} u^{s-1} du =
#     s^j Gamma(s), the roots of beta_3, beta_4, beta_5, and that q_2, q_3, q_4 have no root in [9, oo);
#  2. the zeros: the first 1700 by Arb, all on the line, inside the balls of data/zeros_1700.txt, and N(2197.3) = 1700;
#  3. the constants: L_2 = L_3 = L_4 = 1, upper bounds for J_2, J_3, J_4, the archimedean constants I_1 < 0.404 and
#     C = 0.05 of idsd_height.py (recomputed here by the same method), and the bound C = 1.08 that Binet's formula gives
#     for every r >= 5 with no computation (Remark (b) after the proof);
#  4. the large t range: W^[3](t) > 0 for t >= tau_1 = (3 + sqrt 5)/(2 gamma_1^2), the zeros above height 2197.3 bounded
#     through the N(T) bound of Hasanalizade, Shen and Wong (J. Number Theory 235 (2022), Corollary 1.2), and the bound
#     of Remark (a) for the zeros above 3e12;
#  5. the small t range: W^[3](t) > 0 for t <= tau_2 = 1/75, from the explicit formula, with the computed constants
#     (0.404, 0.05) and with C = 1.08;
#  6. rung four: W^[4](e^{-4.15}) < 0, from the zeros and, with no zeros, from the explicit formula;
#  7. a second check, in three ranges: explicit formula cells of width 0.001 on y = -log t in [4.30, 8] for W^[2] and
#     W^[3], and a monotone bound for y >= 8.
# Every inequality is evaluated in ball arithmetic (python-flint, Arb) and counts only if it holds for the whole ball.
# Outside the certificate the script compares the explicit formula, with the archimedean integral by rigorous quadrature,
# with the zero sums at five points, gives point values of W^[j] for j <= 6, and scans W^[j], j <= 6, in double precision
# on y in [0, 10] (negative sets, with the sign at both ends of every bracket certified from the zeros, and the ratios
# W^[2]/W and W^[3]/W); M_j/M_0 at y = 32 shows the approach to the limits 1/4 and 1/8.
# Run on Hopper (about a minute, one core), from a folder containing zeros_1700.txt: python3 bdlp_ladder.py; writes
# bdlp_ladder.json there.
import json
import math
import time
from fractions import Fraction

import numpy as np
from flint import arb, acb, ctx

ctx.prec = 128
pi, half = arb.pi(), arb(1) / 2
out, holds = {}, []
T0 = time.time()


def s(x, n=12):
    return x.str(n, radius=False) if isinstance(x, arb) else x


def check(name, cond, **vals):
    holds.append(bool(cond))
    out[name] = dict(holds=bool(cond), **{k: s(v) for k, v in vals.items()})
    print(name, bool(cond), round(time.time() - T0, 1), flush=True)


def ball(a, b):  # a ball containing [a, b], for floats or arbs a <= b
    return arb(a).union(arb(b))


def A(c):
    return arb(c.numerator) / c.denominator if isinstance(c, Fraction) else arb(c)


def pv(c, x):  # the polynomial with coefficients c (low to high) at x
    r = x * 0
    for a in reversed(c):
        r = r * x + A(a)
    return r


def pv_float(c, x):
    r = 0.0
    for a in reversed(c):
        r = r * x + float(a)
    return r


# 1. The polynomials.
def beta_polys(J):
    B = [[1]]
    for _ in range(J):
        c, new = B[-1], [0] * (len(B[-1]) + 1)
        for k, a in enumerate(c):
            new[k + 1] += a
            new[k] -= k * a
        B.append(new)
    return B


def q_polys(J):
    Q = [[Fraction(1)]]
    for _ in range(J):
        c, new = Q[-1], [Fraction(0)] * (len(Q[-1]) + 1)
        for k, a in enumerate(c):
            new[k] += a / 2 + k * a
            new[k + 1] -= a
        Q.append(new)
    return Q


BETA, Q = beta_polys(8), q_polys(5)
BABS = [[abs(a) for a in c] for c in BETA]
stirling = [[1]]
for j in range(1, 9):
    prev = stirling[-1] + [0]
    stirling.append([0] + [i * prev[i] + prev[i - 1] for i in range(1, j + 1)])
touchard = all(BETA[j][i] == (-1) ** (j - i) * stirling[j][i] for j in range(9) for i in range(j + 1))


def rising(x, i):
    r = 1
    for k in range(i):
        r *= x + k
    return r


mellin = all(sum(BETA[j][i] * rising(x, i) for i in range(j + 1)) == x ** j for j in range(9) for x in range(1, 12))


def roots_in(c, lo, hi, step):  # sign changes of a polynomial on a grid (a grid point may be a root), each certified on a
    # ball of width 2e-12 by Arb
    f = lambda x: pv_float(c, x)
    rts, ok = [], True
    n = int(round((hi - lo) / step))
    for i in range(n):
        a, b = lo + i * step, lo + (i + 1) * step
        if f(a) == 0 and i > 0:
            r = a
        elif f(a) * f(b) < 0:
            for _ in range(200):
                m_ = (a + b) / 2
                if f(a) * f(m_) <= 0:
                    b = m_
                else:
                    a = m_
            r = (a + b) / 2
        else:
            continue
        ya, yb = pv(c, arb(r) - arb('1e-12')), pv(c, arb(r) + arb('1e-12'))
        ok = ok and bool((ya < 0 and yb > 0) or (ya > 0 and yb < 0))
        rts.append(r)
    return rts, ok


def no_root_beyond(c, v0):  # c(v0 + w), as a polynomial in w, has coefficients of one sign: no root in [v0, oo)
    n = len(c)
    sh = [Fraction(0)] * n
    for k, a in enumerate(c):
        for i in range(k + 1):
            sh[i] += Fraction(a) * math.comb(k, i) * Fraction(v0) ** (k - i)
    signs = {(x > 0) - (x < 0) for x in sh if x != 0}
    return len(signs) == 1


rb = {}
rok = True
for j in (3, 4, 5):
    r, ok = roots_in(BETA[j][1:], 0.0, 20.0, 0.001)
    rb[j], rok = r, rok and ok and len(r) == j - 1
out['beta'] = {j: BETA[j] for j in range(7)}
out['q'] = {j: [str(a) for a in Q[j]] for j in range(5)}
check('polynomials', touchard and mellin and rok and all(no_root_beyond(Q[j], 9) for j in (2, 3, 4)),
      roots_beta3=[round(x, 6) for x in rb[3]], roots_beta4=[round(x, 6) for x in rb[4]],
      roots_beta5=[round(x, 4) for x in rb[5]], roots_beta3_exact='(3 -+ sqrt 5)/2',
      note='beta_j(u) = sum_i (-1)^(j-i) S(j,i) u^i = (-1)^j B_j(-u) (Touchard), int beta_j e^-u u^(s-1) du = s^j Gamma(s) '
           'checked exactly at s = 1, ..., 11 for j <= 8; q_2, q_3, q_4 have no root in [9, oo)')

# 2. The zeros.
NZ = 1700
Z = acb.zeta_zeros(1, NZ)
tab = [arb(l.strip(), '1e-66') for l in open('zeros_1700.txt') if l.strip()]
on_line = all(float((z.real - half).mid()) == 0.0 and float(z.real.rad()) == 0.0 for z in Z)
inside = all(z.imag.overlaps(g) for z, g in zip(Z, tab))
GAM = arb('2197.3')
NG = GAM.zeta_nzeros()
gam = [z.imag for z in Z]
g2 = [g * g for g in gam]
check('zeros', on_line and inside and NG.overlaps(arb(NZ)) and float(NG.rad()) == 0.0 and gam[-1] < GAM,
      count=NZ, N_of_2197_3=NG, gamma_1=gam[0], gamma_2=gam[1], gamma_1700=gam[-1],
      max_abs_arb_minus_table=max(abs(float((z.imag - g).mid())) for z, g in zip(Z, tab)))

# 3. The constants. L_j = sup |beta_j(u)/u| e^{-u} = 1 for j = 2, 3, 4: the function p(u) = (beta_j(u)/u) e^{-u} has
# |p(0)| = 1, tends to 0, and is monotone between the zeros of its derivative, where |p| < 1 (certified on balls).
Lok, crit = True, {}
for j in (2, 3, 4):
    p = BETA[j][1:]
    dp = [(k + 1) * p[k + 1] for k in range(len(p) - 1)] + [0]
    r_ = [a - b for a, b in zip(dp, p)]  # p' - p, so that (p e^{-u})' = (p' - p) e^{-u}
    while r_ and r_[-1] == 0:
        r_.pop()
    cr, ok = roots_in(r_, 0.0, 60.0, 0.001)
    cauchy = 1 + max(abs(a) for a in r_[:-1]) / abs(r_[-1])  # every root of p' - p is below this
    vals = [abs(pv(p, ball(x - 1e-12, x + 1e-12)) * (-ball(x - 1e-12, x + 1e-12)).exp()) for x in cr]
    Lok = Lok and ok and len(cr) == len(r_) - 1 and cauchy < 60 and abs(p[0]) == 1 and all(v < 1 for v in vals)
    crit[j] = [round(x, 4) for x in cr], [s(v, 6) for v in vals]


def J_upper(j):  # int_0^oo |beta_j(u)| e^{-u} u^{-3/2} du = 2 int_0^oo |p_j(w^2)| e^{-w^2} dw, split at the roots
    p = BETA[j][1:]
    pa = [abs(a) for a in p]
    rts, ok = roots_in(p, 0.0, 20.0, 0.001)
    assert ok and len(rts) == len(p) - 1  # every root of p is real, positive and below 20 (here j <= 4)
    eps = 1e-12
    ws = [arb(0)] + [arb(r).sqrt() for r in rts] + [arb(12)]
    tot = arb(0)
    for a, b in zip(ws[:-1], ws[1:]):
        tot += abs(acb.integral(lambda w, an: 2 * pv(p, w * w) * (-(w * w)).exp(), a, b).real)
    slack = arb(0)
    for r in rts:  # the root lies within eps of r: splitting at sqrt(r) errs by at most twice the integral over the gap
        wa, wb = (arb(r) - eps).sqrt(), (arb(r) + eps).sqrt()
        W_ = wa.union(wb)
        slack += 2 * (wb - wa) * abs(2 * pv(p, W_ * W_) * (-(W_ * W_)).exp())
    tail = 2 * pv(pa, arb(144)) * arb(-144).exp()  # w >= 12: |p(w^2)| e^{-w^2} <= pa(144) e^{-144} e^{-(w-12)}
    return arb((tot + slack + tail).upper())


J = {j: J_upper(j) for j in (2, 3, 4)}


def m(r):  # m(r) = Re psi(1/4 + i r/2) - log pi, for a real ball r
    return acb(arb(1) / 4, r / 2).digamma().real - pi.log()


def dev(r):  # m(r) - log(r/2pi)
    return m(r) - (r / (2 * pi)).log()


def ddev(r):  # dev'(r) = -Im psi'(1/4 + i r/2)/2 - 1/r, with psi'(z) = zeta(2, z)
    return -acb(2).zeta(acb(arb(1) / 4, r / 2)).imag / 2 - 1 / r


def sup_abs(f, df, a, b):  # an upper bound for |f| on [a, b]: the smaller of the direct enclosure and the mean value form
    X, c, h = ball(a, b), (arb(a) + arb(b)) / 2, (arb(b) - arb(a)) / 2
    u, v = abs(f(X)).upper(), (abs(f(c)) + abs(df(X)) * h).upper()
    return u if u < v else v


# I_1 = int_0^5 r^2 |dev| dr: on [0, a0] by |dev| <= max|m| + log(2pi/r), then upper sums over pieces, geometric to 0.05
# and of width 1e-4 beyond (the method of idsd_height.py).
a0 = 2.0 ** -20
I1 = (abs(m(ball(0.0, a0))).upper() * a0 + a0 * (1 + (2 * pi / a0).log())) * a0 ** 2
pts = [a0]
while pts[-1] < 0.05:
    pts.append(pts[-1] * 1.001)
while pts[-1] < 5:
    pts.append(min(5.0, round(pts[-1] + 1e-4, 10)))
for a, b in zip(pts[:-1], pts[1:]):
    I1 += sup_abs(dev, ddev, a, b) * (arb(b) - arb(a)) * arb(b) ** 2


def sup_cover(f, df, a, b, bound, depth=0):  # sup_abs on [a, b], bisected while it is not below bound
    u = sup_abs(f, df, a, b)
    if u < bound or depth == 12:
        return u
    c = (a + b) / 2
    u1, u2 = sup_cover(f, df, a, c, bound, depth + 1), sup_cover(f, df, c, b, bound, depth + 1)
    return u1 if u1 > u2 else u2


C = arb('0.05')
supC = arb(0)
for i in range(5000):
    u = sup_cover(lambda r: r * r * dev(r), lambda r: 2 * r * dev(r) + r * r * ddev(r), 5 + i / 1000, 5 + (i + 1) / 1000, C)
    if u > supC:
        supC = u


def tail_v(a, k):  # int_a^oo v^k dv/(e^{2 pi v} - 1) <= e^{-ca}/(1 - e^{-ca}) int_a^oo v^k e^{-c(v-a)} dv, c = 2pi; k = 1, 3
    c = 2 * pi
    e = (-c * a).exp() / (1 - (-c * a).exp())
    if k == 1:
        return e * (a / c + 1 / c ** 2)
    return e * (a ** 3 / c + 3 * a ** 2 / c ** 2 + 6 * a / c ** 3 + 6 / c ** 4)


def binet_tail(r):  # idsd_height.py: the bound for r^2 |dev(r)| when r >= 10, using the cancellation 1/8 - 1/2 + 1/3
    r = arb(r)
    K = 1 / (15 * (3 * r * r - 1)) + 4 / r * tail_v(r / 4, 3)
    return arb(1) / 24 + 1 / (64 * r * r) + 1 / (4 * r * r) + 1 / (48 * r ** 4) + 8 * K


# The analytic bound for every r >= 5 (Remark (b)): with z = 1/4 + ir/2, Binet's formula gives dev(r) =
# (1/2) log(1 + 1/4r^2) - 2/(1 + 4r^2) - 2 Re int_0^oo v dv/((v^2 + z^2)(e^{2 pi v} - 1)), and |v^2 + z^2| is at least
# (3r^2 - 1)/16 for v <= r/4 and at least r/4 beyond, so r^2 |dev| <= 1/8 + 1/2 + 4r^2/(3(3r^2 - 1)) + 8r T_1(r/4),
# T_1(a) = int_a^oo v dv/(e^{2 pi v} - 1). Each term is nonincreasing in r (the last because a^2 > a/c + 1/c^2 for
# a = r/4 >= 1.25), so its value at r = 5 bounds r^2 |dev| for every r >= 5.
r5 = arb(5)
CD_pieces = [arb(1) / 8, half, 4 * r5 ** 2 / (3 * (3 * r5 ** 2 - 1)), 8 * r5 * tail_v(r5 / 4, 1)]
CD5 = sum(CD_pieces)
CD = arb('1.08')
a5 = r5 / 4
c2p = 2 * pi
check('constants', Lok and I1 < arb('0.404') and supC < C and binet_tail(10) < C and CD5 < CD
      and a5 ** 2 > a5 / c2p + 1 / c2p ** 2,
      L_234='1 (attained at u = 0)', critical_points_and_values=crit, J2_upper=J[2], J3_upper=J[3], J4_upper=J[4],
      I1_upper=arb(I1.upper()), I1_stated='0.404', pieces=len(pts) - 1, sup_r2dev_on_5_10=supC,
      binet_tail_at_10=binet_tail(10), C_stated='0.05', analytic_pieces_at_5=[s(x, 6) for x in CD_pieces],
      analytic_bound_at_5=CD5, analytic_stated='1.08')
I1b = arb('0.404')


def E_bound(j, t, Cc=C):  # |A_j - M_j| <= (1/2pi)(L_j I_1 t + C sqrt(t) J_j/2), L_j = 1, nondecreasing in t
    return (I1b * t + Cc * t.sqrt() / 2 * J[j]) / (2 * pi)


# The explicit formula.
alpha0 = -arb.const_euler() / 2 - (4 * pi).log()


def M(j, t):  # 2^{-j} s_0(t) (ell + alpha_0 + j)
    return arb(2) ** (-j) / (4 * (pi * t).sqrt()) * (-(t.log()) / 2 + alpha0 + j)


def polar(j, t):
    return pv(BETA[j], -t / 4) * (t / 4).exp()


def polar_abs(j, t):  # |beta_j(-t/4)| e^{t/4} <= btilde_j(t/4) e^{t/4}, nondecreasing in t
    return pv(BABS[j], t / 4) * (t / 4).exp()


def vonmangoldt(n):
    for p in range(2, n + 1):
        if n % p == 0:
            k = n
            while k % p == 0:
                k //= p
            return arb(p).log() if k == 1 else None
    return None


NP = 200
LAM = {n: vonmangoldt(n) for n in range(2, NP + 1)}
LAM = {n: v for n, v in LAM.items() if v is not None}
QABS = [[abs(a) for a in c] for c in Q]


def prime_tail(j, t):  # for n > NP: |q_j(v)| e^{-v} <= e^{-v/2} (v >= 40) and e^{-(log n)^2/8t} <= n^{-sigma}, sigma =
    # log NP/(8t); with Lambda(n) <= log n, the sum over n > NP is at most (2 sqrt(pi t))^{-1} int_NP^oo log x x^{1/2-sigma-1}
    # dx, which is below the expression returned (it has NP^{1-sigma} and sigma - 1 for NP^{1/2-sigma} and sigma - 1/2).
    tl = arb(t.lower())
    vmin = arb(NP + 1).log() ** 2 / (4 * arb(t.upper()))
    assert vmin > 40 and pv(QABS[j], arb(40)) < arb(20).exp()
    sig = arb(NP).log() / (8 * arb(t.upper()))
    return (2 * (pi * tl).sqrt()) ** -1 * arb(NP) ** (1 - sig) * (arb(NP).log() / (sig - 1) + 1 / (sig - 1) ** 2)


def P(j, t):  # P_j(t), with the tail as a ball
    acc = t * 0
    for n, L in LAM.items():
        v = arb(n).log() ** 2 / (4 * t)
        acc += L / arb(n).sqrt() * (-v).exp() * pv(Q[j], v)
    acc = acc / (2 * (pi * t).sqrt())
    return acc + arb(0, prime_tail(j, t).upper())


def P_abs(j, t):  # the sum of the moduli of the terms of P_j(t), n <= NP exactly, n > NP by prime_tail
    acc = arb(0)
    for n, L in LAM.items():
        v = arb(n).log() ** 2 / (4 * t)
        acc += L / arb(n).sqrt() * (-v).exp() * abs(pv(Q[j], v))
    return acc / (2 * (pi * t).sqrt()) + prime_tail(j, t)


# The zero sums: the 1700 zeros on the line, and a bound for the zeros above GAM wherever they lie. Such a zero has
# Re a >= gamma^2 - 1/4 and |a| <= gamma^2 + 1/4, so its term is at most f(gamma) = btilde_j((gamma^2 + 1/4) t)
# e^{-(gamma^2 - 1/4) t} in modulus, btilde_j with the moduli of the coefficients; f decreases for (gamma^2 + 1/4) t > j.
# Partial summation against N(T) <= N+(T) = (T/2pi) log(T/2pi e) + 0.1038 log T + 0.2573 log log T + 9.3675 (HSW,
# Corollary 1.2), whose smooth part has derivative at most (1/2pi) log(T/2pi) + 0.4/T, gives the sum over the zeros above GAM
# at most f(GAM)(N+(GAM) - 1700) + int_GAM^oo ((1/2pi) log(T/2pi) + 0.4/T) f(T) dT. For (GAM^2 + 1/4) t >= 2j,
# f(T) <= f(GAM) e^{-(T^2 - GAM^2) t/2} <= f(GAM) e^{-GAM (T - GAM) t}, and log(T/2pi) <= log(GAM/2pi) + (T - GAM)/GAM, so
# the integral is at most f(GAM)(L/(GAM t) + 1/(2 pi GAM^3 t^2)), L = log(GAM/2pi)/2pi + 0.4/GAM.
def N_upper(T):
    return T / (2 * pi) * (T / (2 * pi * arb(1).exp())).log() + arb('0.1038') * T.log() + arb('0.2573') * T.log().log() \
        + arb('9.3675')


def zero_tail(j, t):  # t > 0 an arb, the lower end of the range of t
    v = (GAM ** 2 + half / 2) * t
    assert v > 2 * max(j, 1)
    f = pv(BABS[j], v) * (-(GAM ** 2 - half / 2) * t).exp()
    L = (GAM / (2 * pi)).log() / (2 * pi) + arb('0.4') / GAM
    return f * ((N_upper(GAM) - NZ) + L / (GAM * t) + 1 / (2 * pi * GAM ** 3 * t ** 2))


def W_zeros(j, t):  # W^[j](t) from the 1700 zeros, with the tail as a ball; t a ball
    acc = t * 0
    for a in g2:
        u = a * t
        acc += pv(BETA[j], u) * (-u).exp()
    return acc + arb(0, zero_tail(j, arb(t.lower())).upper())


def xb(y0, y1=None):  # the ball of t = e^{-y} for y in [y0, y1]
    if y1 is None:
        return (-arb(y0)).exp()
    return ball(arb((-arb(y1)).exp().lower()), arb((-arb(y0)).exp().upper()))


# 4. Large t. tau_1 = (3 + sqrt 5)/(2 gamma_1^2): for t >= tau_1 every one of the 1700 zeros contributes beta_3(gamma^2 t)
# e^{-gamma^2 t} >= 0, and the second at least beta_3(gamma_2^2 tau_1) e^{-gamma_2^2 t}, since beta_3 increases beyond its
# largest root (3 + sqrt 5)/2. The tail B(t) = zero_tail(3, t) has e^{gamma_2^2 t} B(t) decreasing for t >= tau_1: its
# factors are K(t), decreasing, and btilde_3(v) e^{-(GAM^2 - 1/4 - gamma_2^2) t}, whose logarithmic derivative is at most
# 3/t - (GAM^2 - 1/4 - gamma_2^2) < 0.
tau1 = (3 + arb(5).sqrt()) / (2 * g2[0])
sec = pv(BETA[3], g2[1] * tau1)
B1 = zero_tail(3, tau1)
dec = 3 / tau1 < GAM ** 2 - half / 2 - g2[1]
lead = (g2[1] * tau1).exp() * B1
check('large_t', sec > arb('99.3') and lead < arb(10) ** -27458 and dec and g2[1] * tau1 > (3 + arb(5).sqrt()) / 2,
      tau_1=tau1, y_I=-tau1.log(), gamma2sq_tau1=g2[1] * tau1, beta3_at_gamma2sq_tau1=sec,
      tail_at_tau1_log10=B1.log() / arb(10).log(), e_gamma2sq_tau1_times_tail_log10=lead.log() / arb(10).log(),
      W3_at_tau1_zeros=W_zeros(3, tau1))

# Remark (a) after the proof: with the zeros below TPT = 3e12 on the line (Platt and Trudgian) and N(r) <= r^2 for r >= 200
# (Backlund), the zeros above TPT contribute at most f(TPT)(TPT^2 + 2/t) in modulus for t >= tau_1, with f as in zero_tail
# for j = 3: partial summation bounds their sum by -int_TPT^oo f'(r) r^2 dr = f(TPT) TPT^2 + int_TPT^oo 2r f(r) dr, and
# f(r) <= f(TPT) e^{-(r^2 - TPT^2) t/2}. The bound, and e^{gamma_2^2 t} times it, decrease for t >= tau_1 (as for B(t)).
TPT = 3 * arb(10) ** 12
remA = (pv(BABS[3], (TPT ** 2 + half / 2) * tau1)).log() - (TPT ** 2 - half / 2) * tau1 + (TPT ** 2 + 2 / tau1).log()
check('remark_a', remA + g2[1] * tau1 < -arb(10) ** 23 and 3 / tau1 < TPT ** 2 - half / 2 - g2[1],
      log_bound_at_tau1=remA, log_e_gamma2sq_tau1_times_bound=remA + g2[1] * tau1)

# 5. Small t. For t <= tau_2 = 1/75: M_3 decreases (ell + alpha_0 + 4 > 0), the polar modulus and E_3 increase, and every
# prime term of P_3 increases in modulus, because v_n >= v_2(tau_2) > 9, where q_3 < 0 < q_4, and -t d/dt of the n-th term of
# P_3 is the n-th term of P_4. So W^[3](t) >= M_3(tau_2) - |polar| - E_3 - sum |prime terms|, all at tau_2.
tau2 = arb(1) / 75
v2 = arb(2).log() ** 2 / (4 * tau2)
mono = (-(tau2.log()) / 2 + alpha0 + 4 > 0) and v2 > 9 and tau2 > tau1
pieces = dict(M3=M(3, tau2), polar=polar_abs(3, tau2), E3=E_bound(3, tau2), P3=P_abs(3, tau2),
              E3_analytic=E_bound(3, tau2, CD))
low = pieces['M3'] - pieces['polar'] - pieces['E3'] - pieces['P3']
low_an = pieces['M3'] - pieces['polar'] - pieces['E3_analytic'] - pieces['P3']
check('small_t', mono and low > arb('0.29') and low_an > arb('0.28'), tau_2=tau2, y_2=-tau2.log(), v_2=v2,
      **{k: v for k, v in pieces.items()}, W3_lower=low, W3_lower_analytic=low_an, W3_at_tau2_zeros=W_zeros(3, tau2))

# 6. Rung four at t_0 = e^{-4.15}.
t0 = xb('4.15')
W4z = W_zeros(4, t0)
first5 = [pv(BETA[4], a * t0) * (-(a * t0)).exp() for a in g2[:5]]
ef = polar(4, t0) + M(4, t0) - P(4, t0)
n2 = arb(2).log() / arb(2).sqrt() * (-(arb(2).log() ** 2 / (4 * t0))).exp() * pv(Q[4], arb(2).log() ** 2 / (4 * t0)) \
    / (2 * (pi * t0).sqrt())
up, up_an = ef + E_bound(4, t0), ef + E_bound(4, t0, CD)
r4 = rb[4]
check('rung_four', W4z < 0 and up < arb('-0.109') and up_an < arb('-0.090') and arb(r4[1]) < g2[0] * t0 < arb(r4[2]),
      W4_zeros=W4z.str(22), W_zeros_at_t0=W_zeros(0, t0), gamma1sq_t0=g2[0] * t0,
      first_five_terms=[s(x, 8) for x in first5], polar=polar(4, t0), M4=M(4, t0), P4=P(4, t0), P4_n2=n2,
      centre=ef, E4=E_bound(4, t0), E4_analytic=E_bound(4, t0, CD), W4_upper=up, W4_upper_analytic=up_an)

# 7. The second check, in three ranges. Range I is that of part 4 (y <= y_I = 4.3348...). Range II: cells of width 0.001
# on y in [4.30, 8], the explicit formula evaluated on the ball of t of each cell, E_j at the upper end. Range III, y >= 8:
# the monotone bound of part 5 at t = e^{-8} (q_2 > 0 > q_3 on [9, oo) gives the monotonicity of the terms of P_2).
cells = {}
for j in (2, 3):
    worst, fails = None, 0
    for i in range(3700):
        X = xb(arb(4300 + i) / 1000, arb(4301 + i) / 1000)
        v = polar(j, X) + M(j, X) - P(j, X) - E_bound(j, arb(X.upper()))
        if not v > 0:
            fails += 1
        if worst is None or v.lower() < worst[0].lower():
            worst = (arb(v.lower()), 4.3 + i / 1000)
    cells[j] = dict(cells=3700, fails=fails, least_lower=s(worst[0], 8), at_y=round(worst[1], 3))
x8 = xb('8')
r3 = {}
for j in (2, 3):
    r3[j] = M(j, x8) - polar_abs(j, x8) - E_bound(j, x8) - P_abs(j, x8)
check('three_range', all(c['fails'] == 0 for c in cells.values()) and all(v > 0 for v in r3.values())
      and arb(2).log() ** 2 / (4 * x8) > 9 and -(x8.log()) / 2 + alpha0 + 3 > 0,
      range_II=cells, range_III_lower_W2=r3[2], range_III_lower_W3=r3[3])
out['all_hold'] = all(holds)


# Outside the certificate. (i) The explicit formula with the archimedean integral by rigorous quadrature, against the
# zero sums.
def A_quad(j, t):
    R = (arb(1000) / t).sqrt()

    def f(r, an):
        u = r * r * t
        mm = ((acb(arb(1) / 4) + acb(0, half) * r).digamma() + (acb(arb(1) / 4) - acb(0, half) * r).digamma()) / 2 \
            - pi.log()
        return pv(BETA[j], u) * (-u).exp() * mm
    I = acb.integral(f, 0, R, eval_limit=10 ** 7, depth_limit=10 ** 5).real / (2 * pi)
    tail = (R.log() + 3) * (-(R * R * t) / 2).exp() / (R * t) / (2 * pi)
    return I + arb(0, tail.upper())


xc = {}
for ys in ('4.15', '4.30', '5', '6', '8'):
    t = xb(ys)
    row = {}
    for j in (0, 2, 3, 4):
        aq = A_quad(j, t)
        kz, ke = W_zeros(j, t), polar(j, t) + aq - P(j, t)
        row[j] = dict(zeros=kz.str(15), explicit=ke.str(15), agree=bool(kz.overlaps(ke)), A_minus_M=s(aq - M(j, t), 6),
                      E_bound=s(E_bound(j, t), 6) if j else None)
    xc[ys] = row
    print('crosscheck', ys, round(time.time() - T0, 1), flush=True)
out['crosscheck_explicit_vs_zeros'] = xc
out['point_values_zeros'] = {ys: {j: W_zeros(j, xb(ys)).str(14) for j in range(7)}
                             for ys in ('4.15', '4.30', '4.3349', '5', '6', '7', '8')}

# (ii) Double precision scan of W^[j], j <= 6, on y in [0, 10], from the 1700 zeros (the zeros above contribute below
# 1e-90 there); negative sets, with the sign at the grid points on either side of every change certified from the zeros.
g2f = np.array([float(a.mid()) for a in g2])
ys = np.round(np.arange(0, 10001) / 1000.0, 3)
Wv = np.zeros((7, len(ys)))
for i0 in range(0, len(ys), 500):
    xs = np.exp(-ys[i0:i0 + 500])
    U = np.outer(xs, g2f)
    E = np.exp(-U)
    for j in range(7):
        pj = np.zeros_like(U)
        for c in reversed(BETA[j]):
            pj = pj * U + c
        Wv[j, i0:i0 + 500] = (pj * E).sum(axis=1)
scan = {}
for j in range(7):
    neg = Wv[j] < 0
    ints, brackets_ok = [], True
    i = 0
    while i < len(ys):
        if neg[i]:
            k = i
            while k + 1 < len(ys) and neg[k + 1]:
                k += 1
            ints.append((float(ys[i]), float(ys[k + 1]) if k + 1 < len(ys) else float(ys[k])))
            for yy, sg in ((ys[max(i - 1, 0)], 1 if i else -1), (ys[i], -1), (ys[k], -1),
                           (ys[min(k + 1, len(ys) - 1)], 1 if k + 1 < len(ys) else -1)):
                v = W_zeros(j, xb(arb(int(round(yy * 1000))) / 1000))
                brackets_ok = brackets_ok and bool(v > 0 if sg > 0 else v < 0)
            i = k + 1
        else:
            i += 1
    scan[j] = dict(negative_on_y=ints, brackets_certified=brackets_ok)
idx = {y: int(round(y * 1000)) for y in (4.3, 5, 6, 7, 8, 9, 10)}
ratios = {str(y): [round(Wv[2, k] / Wv[0, k], 6), round(Wv[3, k] / Wv[0, k], 6)] for y, k in idx.items()}
seg = slice(4300, 10001)
decr = bool(np.all(np.diff(Wv[2, seg] / Wv[0, seg]) < 0) and np.all(np.diff(Wv[3, seg] / Wv[0, seg]) < 0))
x32 = xb('32')
out['scan'] = dict(negative_sets=scan, ratios_W2_W3_over_W=ratios, ratios_decrease_on_4_30_to_10=decr,
                   main_term_ratios_at_y32=[s(M(2, x32) / M(0, x32), 8), s(M(3, x32) / M(0, x32), 8)],
                   min_W2_W3_on_4_30_to_8=[round(float(Wv[2, 4300:8001].min()), 6), round(float(Wv[3, 4300:8001].min()), 6)])
out['seconds'] = round(time.time() - T0, 1)
json.dump(out, open('bdlp_ladder.json', 'w'), indent=1)
print(json.dumps({k: out[k] for k in ('large_t', 'small_t', 'rung_four', 'three_range', 'scan', 'all_hold', 'seconds')},
                 indent=1))
