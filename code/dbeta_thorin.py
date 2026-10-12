# dbeta_thorin.py: the numbers of Section sec:dbeta:beta (ch/ch16.tex), the Dirichlet beta function L(s, chi_4) in Thorin
# form, and of its tables tab:dbeta:levy, tab:dbeta:kron and tab:dbeta:qi. Added for the book in October 2026, from the first
# author's note on the Dirichlet beta function in Thorin form (misc/reviews/nick-dbeta/), recomputed independently in ball
# arithmetic with Arb (python-flint) and, for the truncated prime sums, in double precision with numpy.
# xi(s, chi_4) = (4/pi)^{(s+1)/2} Gamma((s+1)/2) L(s, chi_4) (the normalisation of Section sec:ch16:eta, which is 2/sqrt(pi)
# times the note's Lambda), b_alpha = (xi'/xi)(alpha, chi_4), v_alpha(t) = Re (xi'/xi)(alpha+it, chi_4), V = (log xi)''(1/2).
# Computed: (1) constants: xi(1/2), L(1/2), K_4(0) = 2 eta(i)^3 = Gamma(1/4)^3/(4 pi^{9/4}), b_alpha, V; (2) the kernel: the
# Mellin form xi(1/2+w) = int K_4(u) e^{wu} du by the trapezoid rule (step 1/64 on |u| <= 3, spectrally accurate since K_4 is
# analytic in |Im u| < pi/4), the total mass int eta(ix)^3 dx = 1 and its Laplace transform sech(sqrt(pi u)); (3) the Levy
# form and the sine squared form against the closed forms, with prime powers up to 2e5 and 1e7; (4) the termwise thresholds
# alpha_c for chi_4 and chi_12 and the Kronecker table S(alpha), t_0(alpha), by bisection with certified signs; (5) the
# decomposition at alpha = 3/2 on [0,30]; (6) the Q(i) factorisation xi(s, chi_4) = 4 xi_{Q(i)}(s)/xi(s), with the note's
# xi_{Q(i)}(s) = s(s-1)/2 pi^{-s} Gamma(s) zeta_{Q(i)}(s), a quarter of the xi_K of Section sec:sqfree:nf (there xi_K = xi xi(., chi_4)),
# and the table of the three brackets at alpha = 3/2; (7) the characteristic function of the primal law at t = 5.9 and
# 6.1; (8) the argument form: (1/pi) arg xi(sigma+it, chi_4) by certified continuation against (1/pi) int_0^t v_sigma,
# and the asymptotic formula;
# (9) the reciprocal eta product identity; (10) signs of the prime power weights up to 3000; (11) the Thorin representation
# log(xi(alpha+sqrt u)/xi(alpha)) = (1/pi) int_0^oo log(1+u/t^2) v_alpha(t) dt at (alpha, u) = (2, 1) and (1, 1), by
# Gauss-Legendre with 24 nodes on unit intervals (on (0,1] after t = e^{-x}, which removes the logarithmic singularity) up to
# t = 150 and t = 500 (the piece on (0,1] at 256 bits), with the tail (u/2pi)(log(2T/pi)+1)/T from v_alpha(t) ~ (1/2) log(2t/pi).
# (12) the first passage identity behind the sine squared lemma (Remark rem:dbeta:ig) and (13) the half kernel against the
# incomplete Mellin transform, both with mpmath at 30 digits.
# Run on Hopper (a few minutes): python3 dbeta_thorin.py; writes dbeta_thorin.json in the current folder.
import json
import math
import time

import numpy as np
from flint import acb, acb_series, arb, ctx, dirichlet_char

t0 = time.time()
out = {}


def rec(k, v):
    out[k] = v
    print(k, v, round(time.time() - t0, 1), flush=True)


ctx.prec = 128
CHI4 = dirichlet_char(4, 3)
CHI12 = dirichlet_char(12, 11)
PI = arb.pi()
LOGPI = PI.log()
LOG4PI = (arb(4) / PI).log()
LOG12PI = (arb(12) / PI).log()
LOG2 = arb(2).log()
LOG3 = arb(3).log()


def S(x, n=10):
    return x.str(n, radius=False)


def xi4(s):
    return ((s + 1) / 2 * LOG4PI).exp() * ((s + 1) / 2).gamma() * acb.dirichlet_l(s, CHI4)


def xiz(s):
    return s * (s - 1) / 2 * (-s / 2 * LOGPI).exp() * (s / 2).gamma() * s.zeta()


def xiK(s):
    # xi_{Q(i)}(s) = s(s-1)/2 pi^{-s} Gamma(s) zeta(s) L(s, chi_4)
    return s * (s - 1) / 2 * (-s * LOGPI).exp() * s.gamma() * s.zeta() * acb.dirichlet_l(s, CHI4)


def dlogL(s, chi):
    c = acb_series.dirichlet_l(acb_series([s, 1], prec=2), chi).coeffs()
    return c[1] / c[0]


def dlogzeta(s):
    c = acb_series.zeta(acb_series([s, 1], prec=2)).coeffs()
    return c[1] / c[0]


def LL4(s):
    return LOG4PI / 2 + ((s + 1) / 2).digamma() / 2 + dlogL(s, CHI4)


def LLxi(s):
    return 1 / s + 1 / (s - 1) - LOGPI / 2 + (s / 2).digamma() / 2 + dlogzeta(s)


def LLK(s):
    return 1 / s + 1 / (s - 1) - LOGPI + s.digamma() + dlogzeta(s) + dlogL(s, CHI4)


def A4(a, t):
    return LOG4PI / 2 + ((acb(a, t) + 1) / 2).digamma().real / 2


def S4(a):
    a = acb(a)
    return -dlogzeta(a).real - LOG2 / (arb(2) ** a.real - 1)


def A12(a, t):
    return LOG12PI / 2 + (acb(a, t) / 2).digamma().real / 2


def S12(a):
    return S4(a) - LOG3 / (arb(3) ** arb(a) - 1)


def bisect(F, a, b, tol=1e-13):
    # root of F in [a, b] with certified signs at the ends; returns the final bracket
    fa = F(arb(a))
    fb = F(arb(b))
    assert (fa < 0 and fb > 0) or (fa > 0 and fb < 0), (a, b, fa, fb)
    sa = 1 if fa > 0 else -1
    while b - a > tol:
        m = (a + b) / 2
        fm = F(arb(m))
        if not (fm > 0 or fm < 0):
            break
        if (1 if fm > 0 else -1) == sa:
            a = m
        else:
            b = m
    return a, b


# (1) constants
half = acb(0.5)
xh = xi4(half)
K0 = arb(0.25).gamma() ** 3 / (4 * PI ** arb(2.25))
K0s = 2 * sum((-1) ** n * (2 * n + 1) * (-PI * (2 * n + 1) ** 2 / 4).exp() for n in range(14))
b1 = LL4(acb(1)).real
ser = acb_series([half, 1], prec=3)
logxi = (ser + 1) / 2 * acb(LOG4PI) + ((ser + 1) / 2).lgamma() + acb_series.dirichlet_l(ser, CHI4).log()
V = 2 * logxi.coeffs()[2].real
rec('constants', {'xi(1/2,chi4)': S(xh.real, 15), 'L(1/2,chi4)': S(acb.dirichlet_l(half, CHI4).real, 15),
                  'K4(0)=Gamma(1/4)^3/(4pi^(9/4))': S(K0, 15), 'K4(0) series': S(K0s, 15),
                  'xi(1,chi4)': S(xi4(acb(1)).real, 15), 'xi(3,chi4)': S(xi4(acb(3)).real, 15), 'pi/2': S(PI / 2, 15),
                  'b_1': S(b1, 12), 'b_1/sqrt(pi)': S(b1 / PI.sqrt(), 8), 'V': S(V, 12), 'p1=V/2': S(V / 2, 12),
                  "L'(1,chi4)": S(acb_series.dirichlet_l(acb_series([acb(1), 1], prec=2), CHI4).coeffs()[1].real, 12),
                  "(pi/4)(gamma+2log2+3logpi-4logGamma(1/4))": S(PI / 4 * (arb.const_euler() + 2 * LOG2 + 3 * LOGPI
                                                                           - 4 * arb(0.25).gamma().log()), 12),
                  'b_1 closed: 3log2+(5/2)logpi+gamma/2-4logGamma(1/4)': S(3 * LOG2 + 5 * LOGPI / 2 + arb.const_euler() / 2
                                                                          - 4 * arb(0.25).gamma().log(), 12)})
bt = {}
for a in ['0.55', '0.6', '0.75', '1', '1.5', '2', '2.5', '3.15']:
    b = LL4(acb(arb(a))).real
    bt[a] = {'b': S(b, 10), 'b/(a-1/2)': S(b / (arb(a) - 0.5), 8)}
rec('b_alpha', bt)


# (2) the kernel, the Mellin form and the sech law
def K4(u):
    u = abs(u)
    x = (2 * u).exp()
    return 2 * (1.5 * u).exp() * sum((-1) ** n * (2 * n + 1) * (-PI * (2 * n + 1) ** 2 * x / 4).exp() for n in range(14))


h = arb(1) / 64
us = [arb(k) / 64 for k in range(-192, 193)]
Ks = [K4(u) for u in us]
mel = {}
for w in [acb(0), acb(0.5), acb(0.1, 5), acb(0.1, 20)]:
    trap = h * sum(K * (w * u).exp() for K, u in zip(Ks, us))
    mel[w.str(3, radius=False)] = {'trapezoid': trap.str(14, radius=False), 'direct': xi4(half + w).str(14, radius=False)}
# with x = e^{2u}: eta(ix)^3 dx = K_4(u) e^{u/2} du
mass = h * sum(K * (u / 2).exp() for K, u in zip(Ks, us))
lap = h * sum((-(2 * u).exp()).exp() * K * (u / 2).exp() for K, u in zip(Ks, us))
mel['int eta(ix)^3 dx'] = S(mass, 14)
mel['int e^{-x} eta(ix)^3 dx'] = S(lap, 14)
mel['sech(sqrt(pi))'] = S(1 / PI.sqrt().cosh(), 14)
rec('kernel', mel)


# (3) Levy form and sine squared form against the closed forms
def prime_powers(N):
    isp = np.ones(N + 1, dtype=bool)
    isp[:2] = False
    for p in range(2, int(N ** 0.5) + 1):
        if isp[p]:
            isp[p * p::p] = False
    ns, lp = [], []
    for p in np.nonzero(isp)[0]:
        p = int(p)
        q = p
        while q <= N:
            ns.append(q)
            lp.append(math.log(p))
            q *= p
    ns = np.array(ns, dtype=np.float64)
    lp = np.array(lp)
    nn = ns.astype(np.int64)
    chi = np.where(nn % 2 == 0, 0, np.where(nn % 4 == 1, 1, -1)).astype(np.float64)
    return ns, lp, chi


PPS = {X: prime_powers(X) for X in (200000, 10000000)}
lev = []
for a, s in [(2, 2), (3, 2), (1.5, 1)]:
    A_, S_ = acb(a), acb(s)
    ba = LL4(A_).real
    exact = (xi4(A_ + S_) / xi4(A_)).real.log() - ba * s
    gam = ((A_ + S_ + 1) / 2).lgamma().real - ((A_ + 1) / 2).lgamma().real - S_.real / 2 * ((A_ + 1) / 2).digamma().real
    for X, (ns, lp, chi) in PPS.items():
        ln = np.log(ns)
        pr = float(np.sum(chi * lp / ln * ns ** (-a) * (ns ** (-s) - 1 + s * ln)))
        lev.append({'alpha': a, 's': s, 'X': X, 'b s': S(ba * s, 8), 'gamma piece': S(gam, 10), 'primes': '%.10f' % pr,
                    'total': '%.10f' % (float(gam.mid()) + pr), 'exact': S(exact, 10)})
rec('levy', lev)
sq = []
for a, t in [(2, 1), (2, 7), (1.5, 4), (1.5, 6.020948905)]:
    A_ = acb(a)
    ba = LL4(A_).real
    G = (((A_ + 1 + acb(0, t)) / 2).digamma().real - ((A_ + 1) / 2).digamma()) / 2
    closed = LL4(acb(a, t)).real
    for X, (ns, lp, chi) in PPS.items():
        P = float(np.sum(chi * lp * ns ** (-a) * (1 - np.cos(t * np.log(ns)))))
        sq.append({'alpha': a, 't': t, 'X': X, 'b+G': S(ba + G.real, 10), 'primes': '%.10f' % P,
                   'sum': '%.10f' % (float((ba + G.real).mid()) + P), 'closed v_alpha(t)': S(closed, 10)})
rec('sin2', sq)

# (4) termwise thresholds and the Kronecker table
a4 = bisect(lambda a: A4(a, 0) - S4(a), 2.2, 2.4)
a12 = bisect(lambda a: A12(a, 0) - S12(a), 1.7, 1.95)
rec('alpha_c', {'chi4': [a4[0], a4[1]], 'chi12': [a12[0], a12[1]]})
kr = {}
for a in ['2', '1.75', '1.5', '1.25', '1.1', '1.05']:
    al = arb(a)
    Sa = S4(al)
    lo, hi = bisect(lambda x: A4(al, x.exp()) - Sa, 0.0, 60.0, 1e-12)
    kr[a] = {'S': S(Sa, 8), 't0': '%.6g' % math.exp(lo), 'log t0': '%.6f' % lo,
             'log t0 - 2/(alpha-1)': '%.4f' % (lo - 2 / (float(a) - 1))}
rec('kronecker', kr)

# (5) decomposition at alpha = 3/2 on [0,30], step 0.01, at 64 bits
ctx.prec = 64
a = acb(1.5)
bb = dlogL(a, CHI4).real
ba = LL4(a).real
rows = []
for k in range(0, 3001):
    t = k / 100
    s = acb(1.5, t)
    P = dlogL(s, CHI4).real - bb
    cg = A4(1.5, t) + bb
    rows.append((t, float(P.mid()), float(cg.mid())))
kmin = min(range(len(rows)), key=lambda k: rows[k][1])
kmax = max(range(len(rows)), key=lambda k: rows[k][1])
vmin = min(range(1, len(rows)), key=lambda k: rows[k][1] + rows[k][2])
rec('decomposition_3/2', {'b_3/2': S(ba, 8), 'min P': '%.5f' % rows[kmin][1], 'at t': rows[kmin][0],
                          'b+G there': '%.5f' % rows[kmin][2], 'v there': '%.5f' % (rows[kmin][1] + rows[kmin][2]),
                          'max P': '%.5f' % rows[kmax][1], 'at t (max)': rows[kmax][0],
                          'min v on (0,30]': '%.6f' % (rows[vmin][1] + rows[vmin][2]), 'at t (min v)': rows[vmin][0],
                          'fraction of grid with P<0': '%.3f' % (sum(1 for r in rows[1:] if r[1] < 0) / (len(rows) - 1))})
ctx.prec = 128

# (6) the Q(i) factorisation and the three brackets at alpha = 3/2
fac = {}
for s in [acb(0.3, 4.1), acb(2.5)]:
    fac[s.str(3, radius=False)] = {'xi(s,chi4)': xi4(s).str(14, radius=False),
                                   '4 xi_Q(i)/xi': (4 * xiK(s) / xiz(s)).str(14, radius=False)}
tab = {}
for t in [0, 2, 5, 10, 20]:
    s = acb(1.5, t)
    tab[str(t)] = {'K': S(LLK(s).real, 6), 'xi': S(LLxi(s).real, 6), 'chi4': S(LL4(s).real, 6),
                   'K-xi-chi4': (LLK(s).real - LLxi(s).real - LL4(s).real).str(3)}
fac['table_alpha_3/2'] = tab
rec('factorisation', fac)

# (7) the characteristic function of the primal law across the first ordinate
rec('charfn', {'t=5.9': S((xi4(acb(0.5, 5.9)) / xh).real, 6), 't=6.1': S((xi4(acb(0.5, 6.1)) / xh).real, 6)})


# (8) the argument form, by certified continuation as in dbeta_zeros.py
def ball(z0, z1):
    return acb(arb((z0[0] + z1[0]) / 2, abs(z1[0] - z0[0]) / 2), arb((z0[1] + z1[1]) / 2, abs(z1[1] - z0[1]) / 2))


def pt(z):
    return acb(arb(z[0]), arb(z[1]))


def seg_arg(f, z0, z1, depth=0):
    W = f(ball(z0, z1))
    m = W.mid()
    r = (W.real.rad() ** 2 + W.imag.rad() ** 2).sqrt()
    if abs(m) > 1.1 * r:
        return (f(pt(z1)) / f(pt(z0))).arg()
    if depth > 40:
        raise RuntimeError('piece not resolved near %r' % (z0,))
    zm = ((z0[0] + z1[0]) / 2, (z0[1] + z1[1]) / 2)
    return seg_arg(f, z0, zm, depth + 1) + seg_arg(f, zm, z1, depth + 1)


def vert_arg(f, sg, t, h=1 / 16):
    n = int(round(t / h))
    tot = arb(0)
    for k in range(n):
        tot += seg_arg(f, (sg, k * t / n), (sg, (k + 1) * t / n))
    return tot


GLX, GLW = np.polynomial.legendre.leggauss(24)


def int_v(sg, t):
    tot = arb(0)
    for k in range(int(t)):
        for x, w in zip(GLX, GLW):
            tot += float(w) / 2 * LL4(acb(sg, k + (float(x) + 1) / 2)).real
    return tot


ar = {}
for sg, t in [(1.5, 10), (1.5, 25), (2.0, 30)]:
    ar['(%s,%s)' % (sg, t)] = {'arg/pi': S(vert_arg(xi4, sg, t) / PI, 12), 'int v/pi': S(int_v(sg, t) / PI, 12)}
L4 = lambda s: acb.dirichlet_l(s, CHI4)
for t in [50, 100, 200]:
    sg = 0.75
    U = vert_arg(xi4, sg, t) / PI
    aL = vert_arg(L4, sg, t) / PI
    T_ = arb(t)
    main = T_ / (2 * PI) * (2 * T_ / (PI * arb(1).exp())).log() + arb(sg) / 4
    ar['asymptotic t=%d' % t] = {'U': S(U, 10), 'main + arg L/pi': S(main + aL, 10), 'difference': (U - main - aL).str(3),
                                 't*difference': (T_ * (U - main - aL)).str(3)}
rec('argument', ar)

# (9) the reciprocal eta product: e^{-pi x/4}/eta(ix)^3 = prod (1-q^m)^{-3} = exp(3 sum sigma_{-1}(N) q^N), q = e^{-2 pi x}
x = arb(0.3)
q = (-2 * PI * x).exp()
prod = arb(1)
for m in range(1, 400):
    prod *= (1 - q ** m) ** (-3)
sm = arb(0)
for N in range(1, 400):
    sm += sum(arb(1) / d for d in range(1, N + 1) if N % d == 0) * q ** N
eta3x = sum((-1) ** n * (2 * n + 1) * (-PI * (2 * n + 1) ** 2 / (4 * x)).exp() for n in range(14)) / x ** 1.5
rec('recip_product', {'prod': S(prod, 15), 'exp(3 sum)': S((3 * sm).exp(), 15),
                      'e^{-pi x/4}/eta(ix)^3': S((-PI * x / 4).exp() / eta3x, 15)})

# (10) signs of the prime power weights chi_4(n) Lambda(n), n <= 3000
ns, lp, chi = prime_powers(3000)
rec('weights_3000', {'prime powers': int(len(ns)), 'negative': int(np.sum(chi < 0)), 'zero (powers of 2)': int(np.sum(chi == 0)),
                     'positive': int(np.sum(chi > 0))})

# (11) the Thorin representation on the lines alpha = 2 and alpha = 1, at u = 1
ctx.prec = 64
lf = {}
for a, u in [(2, 1), (1, 1)]:
    ua = arb(u)
    # (0,1]: t = e^{-x}, x in [0, 60]; at 256 bits, since L(1+it, chi_4) for t near e^{-60} is a difference of two Hurwitz
    # zeta values of size 1/t and loses about log2(1/t) bits
    ctx.prec = 256
    head = arb(0)
    for k in range(60):
        for x, w in zip(GLX, GLW):
            xx = k + (float(x) + 1) / 2
            t = arb(-xx).exp()
            head += float(w) / 2 * (1 + ua / t ** 2).log() * LL4(acb(arb(a), t)).real * t
    ctx.prec = 64
    body = arb(0)
    res = {}
    for k in range(1, 500):
        for x, w in zip(GLX, GLW):
            t = k + (float(x) + 1) / 2
            body += float(w) / 2 * (1 + ua / arb(t) ** 2).log() * LL4(acb(a, t)).real
        if k + 1 in (150, 500):
            T_ = arb(k + 1)
            tail = ua / (2 * PI) * ((2 * T_ / PI).log() + 1) / T_
            res['T=%d' % (k + 1)] = S((head + body) / PI + tail, 10)
    exact = (xi4(acb(a) + ua.sqrt()) / xi4(acb(a))).real.log()
    lf['(alpha,u)=(%d,%d)' % (a, u)] = dict(res, exact=S(exact, 10))
rec('thorin_representation', lf)
ctx.prec = 128

# (12) the first passage identity behind the sine squared lemma: for s, x > 0,
# e^{-sx} + sx - 1 = int_0^oo (1 - e^{-s^2 t/2})(1 - e^{-x^2/(2t)}) x (2 pi t^3)^{-1/2} dt, at (s, x) = (1.3, 0.8), by mpmath's
# tanh-sinh quadrature at 30 digits.
import mpmath  # noqa: E402

mpmath.mp.dps = 30
s_, x_ = mpmath.mpf('1.3'), mpmath.mpf('0.8')
lhs = mpmath.exp(-s_ * x_) + s_ * x_ - 1
rhs = mpmath.quad(lambda t: mpmath.expm1(-s_ ** 2 * t / 2) * mpmath.expm1(-x_ ** 2 / (2 * t)) * x_
                  / mpmath.sqrt(2 * mpmath.pi * t ** 3), [0, 0.1, 1, 10, 100, mpmath.inf])
rec('first_passage_identity', {'(s,x)': '(1.3,0.8)', 'e^{-sx}+sx-1': mpmath.nstr(lhs, 15),
                               'integral': mpmath.nstr(rhs, 15)})

# (13) the half kernel int_0^oo e^{-wu} K_4(u) du at w = 0.3, by Gauss-Legendre (24 nodes on panels of width 1/16 on [0,3],
# beyond which K_4 < e^{-300}), against the incomplete Mellin transform A(s) = int_1^oo eta(ix)^3 x^{(s+1)/2} dx/x at
# s = 1/2 - w, summed over the Jacobi series with the upper incomplete gamma function:
# A(s) = sum_{n>=0} (-1)^n (2n+1) (pi (2n+1)^2/4)^{-(s+1)/2} Gamma((s+1)/2, pi (2n+1)^2/4); and xi(s, chi_4) = A(s) + A(1-s).
ctx.prec = 128
w_ = arb('0.3')
hk = arb(0)
for k in range(48):
    for xg, wg in zip(GLX, GLW):
        u = (arb(k) + (float(xg) + 1) / 2) / 16
        hk += float(wg) / 32 * K4(u) * (-w_ * u).exp()


def A_inc(s):
    a = (s + 1) / 2
    return mpmath.fsum((-1) ** n * (2 * n + 1) * (mpmath.pi * (2 * n + 1) ** 2 / 4) ** (-a)
                       * mpmath.gammainc(a, mpmath.pi * (2 * n + 1) ** 2 / 4) for n in range(12))


sA = mpmath.mpf('0.5') - mpmath.mpf('0.3')
rec('half_kernel', {'w': '0.3', 'int_0^oo e^{-wu} K_4(u) du': S(hk, 14), 'A(1/2-w)': mpmath.nstr(A_inc(sA), 14),
                    'A(s)+A(1-s) at s=0.2': mpmath.nstr(A_inc(sA) + A_inc(1 - sA), 14),
                    'xi(0.2, chi_4)': S(xi4(acb(arb('0.2'))).real, 14)})
json.dump(out, open('dbeta_thorin.json', 'w'), indent=1)
print('done', round(time.time() - t0, 1))
