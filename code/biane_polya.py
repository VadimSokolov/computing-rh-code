# biane_polya.py: Polya's function P(t) = K_{it/2}(2 pi), its zeros, its hitting time, its index side clock, and the other
# numbers of Section ch:biane (ch/biane.tex, Sections sec:biane:hit to sec:biane:kg) and of Table tab:biane:zeros. Added for
# the book in October 2026, from the first author's note on Biane's paper on Polya (misc/reviews/nick-biane/), recomputed
# independently. K_mu is the standard Macdonald function, K_mu(x) = int_0^oo exp(-x cosh w) cosh(mu w) dw.
# Zeros of P. P(t) = K_{i nu}(2 pi) with nu = t/2, and K_{i nu}(x) = pi |I_{i nu}(x)| sin(Psi(nu)) / sinh(pi nu) with the
# continuous phase Psi(nu) = Im log Gamma(1 + i nu) - nu log(x/2) - arg 0F1(; 1 + i nu; x^2/4), Psi(0) = 0. So the zeros are
# exactly the solutions nu_n of Psi(nu_n) = n pi, n = 1, 2, ..., once Psi is increasing (checked below for nu >= 8, and
# Psi < pi before the first zero). The script solves Psi = n pi by the secant method at 40 digits, checks the phase formula
# against Arb's K_{i nu}(2 pi), encloses the first 40 zeros by sign changes of Arb balls, and counts the sign changes of
# K_{i nu}(2 pi) on a grid of step 0.02 in nu <= 150 as an independent count.
# Zeros of Xi*(t) = 4 pi^2 (K_{9/4+it/2}(2 pi) + K_{9/4-it/2}(2 pi)) = 8 pi^2 Re K_{9/4+it/2}(2 pi) by sign changes of Arb
# balls on a grid of step 0.02 in (0, 1000] and bisection in ball arithmetic; Polya proved that they are all real and simple.
# Zeros of Xi from Arb (acb.zeta_zeros). Counts below 60, 100, 240; products over the zeros against the closed forms on the
# index side (Theorem thm:biane:indexP, the Xi* and Xi analogues), with a tail from the smooth counting function; the mean
# and variance of the clock tau_oo and of the clock from the finite start R_0 = 2 pi e^{2.5}; and the identities of the
# hitting time sections (Mellin transforms, Levy and Thorin densities, heat kernels on H^3 and S^3, the curvature family,
# the Klein Gordon equation, the sine squared images of the two clocks), each computed two ways.
# Run from any folder: python3 biane_polya.py; writes biane_polya.json there (about seven minutes on one core).
import json
import random
import time

import mpmath as mp
import sympy as sp
from flint import acb, arb, ctx

t0 = time.time()
out = {}


def rec(k, v):
    out[k] = v
    print(k, v, round(time.time() - t0, 1), flush=True)


def s(x, n=12):
    return mp.nstr(x, n)


mp.mp.dps = 40
PI = mp.pi
ctx.prec = 200
API = arb.pi()
Z2PI = acb(2 * API)


def A(x):
    return arb(mp.nstr(mp.mpf(x), 45))


def K_arb(nu_re, nu_im, prec=200):
    old = ctx.prec
    ctx.prec = prec
    v = acb(2 * arb.pi()).bessel_k(acb(A(nu_re), A(nu_im)))
    ctx.prec = old
    return v


def sign_arb(v):
    if v > 0:
        return 1
    if v < 0:
        return -1
    return 0


def KiP(nu):                     # K_{i nu}(2 pi), Arb ball, real part, with a precision retry
    for prec in (200, 600, 1500):
        v = K_arb(0, nu, prec).real
        if sign_arb(v) != 0:
            return v
    return v


def Xs_arb(t):                   # Xi*(t) = 8 pi^2 Re K_{9/4+it/2}(2 pi), Arb ball
    for prec in (200, 600, 1500):
        v = K_arb(mp.mpf(9) / 4, mp.mpf(t) / 2, prec).real
        if sign_arb(v) != 0:
            return 8 * API ** 2 * v
    return 8 * API ** 2 * v


# ---------------------------------------------------------------- the phase and the zeros of P
def F0(nu):
    return mp.hyp0f1(1 + 1j * nu, PI ** 2)


def Psi(nu):
    return mp.im(mp.loggamma(1 + 1j * nu)) - nu * mp.log(PI) - mp.arg(F0(nu))


# continuity of the branch: unwrap Psi on [0, 12] with step 0.005 and compare with the principal formula
prev = mp.mpf(0)
offset = 0
worst_jump = 0
nu = mp.mpf(0)
mins = []
while nu < 12:
    nu += mp.mpf('0.005')
    v = Psi(nu) + 2 * PI * offset
    while v - prev > PI:
        offset -= 1
        v -= 2 * PI
    while v - prev < -PI:
        offset += 1
        v += 2 * PI
    worst_jump = max(worst_jump, abs(v - prev))
    prev = v
rec('phase: branch offset of the principal formula on [0, 12] (multiples of 2 pi)', offset)
rec('phase: largest step of the unwrapped phase on the grid', s(worst_jump, 4))
minReF = min(mp.re(F0(mp.mpf(x) / 4)) for x in range(32, 4001))
rec('phase: least Re 0F1(;1+i nu;pi^2) on [8, 1000] (positive: principal arg is continuous there)', s(minReF, 6))
# check the representation K = pi |I| sin Psi / sinh(pi nu) against Arb at some points
chk = []
for nu in [mp.mpf('0.5'), mp.mpf(3), mp.mpf('9.5'), mp.mpf(25), mp.mpf('101.3'), mp.mpf(700)]:
    I = (PI ** (1j * nu)) / mp.gamma(1 + 1j * nu) * F0(nu)
    lhs = PI * abs(I) * mp.sin(Psi(nu) + 2 * PI * offset) / mp.sinh(PI * nu)
    ref = mp.mpf(KiP(nu).mid().str(40, radius=False))
    chk.append([s(nu, 6), s(ref, 10), s(abs(lhs / ref - 1), 3)])
rec('phase formula vs Arb K_{i nu}(2 pi): nu, K, relative difference', chk)


def Psi_c(nu):
    return Psi(nu) + 2 * PI * offset


def nu_zero(n, guess):
    return mp.findroot(lambda x: Psi_c(x) - n * PI, guess)


NUMAX = 1000
nus = []
g = mp.mpf('9.77')
n = 1
while True:
    z = nu_zero(n, g)
    if z > NUMAX:
        break
    nus.append(z)
    n += 1
    # next guess from the local slope log(nu/pi)
    g = z + PI / mp.log(z / PI)
rec('P: number of zeros nu_n <= 1000 found by the phase', len(nus))
rec('P: first five nu_n', [s(z, 16) for z in nus[:5]])
mon = all(Psi_c((nus[k] + nus[k + 1]) / 2) > (k + 1) * PI and Psi_c((nus[k] + nus[k + 1]) / 2) < (k + 2) * PI
          for k in range(len(nus) - 1))
dmin = min(mp.diff(Psi_c, mp.mpf(x) / 10) for x in range(80, 10001, 7))
rec('P: phase between consecutive zeros in (n pi, (n+1) pi) for all n', bool(mon))
rec('P: least Psi\'(nu) sampled on [8, 1000] (step 0.7)', s(dmin, 6))
rec('P: Psi(8) / pi (below 1, so no zero below 8)', s(Psi_c(mp.mpf(8)) / PI, 8))
# certify the first 40 zeros by sign changes of Arb balls
cert = 0
for z in nus[:40]:
    a, b = KiP(z - mp.mpf('1e-25')), KiP(z + mp.mpf('1e-25'))
    if sign_arb(a) * sign_arb(b) == -1:
        cert += 1
rec('P: first 40 zeros enclosed to 1e-25 by sign changes of Arb balls', cert)
# independent count by sign changes of K_{i nu}(2 pi) on a grid of step 0.02, nu in (0, 150]
sg = [sign_arb(KiP(mp.mpf(k) / 50)) for k in range(1, 7501)]
nsc = sum(1 for k in range(len(sg) - 1) if sg[k] * sg[k + 1] < 0)
rec('P: sign changes of K_{i nu}(2 pi) on (0, 150], step 0.02', nsc)
rec('P: zeros nu_n < 150 from the phase', sum(1 for z in nus if z < 150))
rec('P: zeros nu_n < 120 from the phase', sum(1 for z in nus if z < 120))
rec('P: zeros nu_n < 50 (t = 2 nu < 100)', sum(1 for z in nus if z < 50))
tP = [2 * z for z in nus]
rec('P: zeros t_n = 2 nu_n in (0, 60)', [s(t, 10) for t in tP if t < 60])
# asymptotics of the phase: Psi(nu) - [nu log(nu/pi) - nu + pi/4] vs (pi^2 - 1/12)/nu, and the Debye phase
asy = []
for nu in [mp.mpf(20), mp.mpf(100), mp.mpf(500)]:
    main = nu * mp.log(nu / PI) - nu + PI / 4
    debye = nu * mp.acosh(nu / (2 * PI)) - mp.sqrt(nu ** 2 - 4 * PI ** 2) + PI / 4
    asy.append([s(nu, 4), s((Psi_c(nu) - main) * nu, 8), s((Psi_c(nu) - debye) * nu, 8)])
rec('phase: nu (Psi - main) -> pi^2 - 1/12, nu (Psi - Debye) -> -1/12', asy)
rec('pi^2 - 1/12', s(PI ** 2 - mp.mpf(1) / 12, 8))

# ---------------------------------------------------------------- zeros of Xi*
TX = 1000
grid = [mp.mpf(k) / 50 for k in range(1, 50 * TX + 1)]
sx = [sign_arb(Xs_arb(t)) for t in grid]
rec('Xi*: grid points with undetermined sign', sum(1 for v in sx if v == 0))
xz = []
for k in range(len(grid) - 1):
    if sx[k] * sx[k + 1] < 0:
        lo, hi, slo = grid[k], grid[k + 1], sx[k]
        for _ in range(48):
            mid = (lo + hi) / 2
            sm = sign_arb(Xs_arb(mid))
            if sm == 0:
                break
            if sm == slo:
                lo = mid
            else:
                hi = mid
        xz.append((lo + hi) / 2)
rec('Xi*: zeros in (0, 1000]', len(xz))
rec('Xi*: zeros below 60, 100, 240', [sum(1 for z in xz if z < T) for T in (60, 100, 240)])
rec('Xi*: zeros in (0, 60)', [s(z, 10) for z in xz if z < 60])

# ---------------------------------------------------------------- zeros of Xi
ctx.prec = 128
GAM = [mp.mpf(z.imag.mid().str(40, radius=False)) for z in acb.zeta_zeros(1, 1700)]
ctx.prec = 200
rec('Xi: zeros below 60, 100, 240', [sum(1 for g in GAM if g < T) for T in (60, 100, 240)])
rec('Xi: zeros in (0, 60)', [s(g, 10) for g in GAM if g < 60])


def N0(T):
    return T / (2 * PI) * mp.log(T / (2 * PI * mp.e))


cnt = []
for T in [60, 100, 240, 1000]:
    cnt.append([T, sum(1 for g in GAM if g < T), sum(1 for z in xz if z < T), sum(1 for t in tP if t < T),
                s(N0(mp.mpf(T)), 6)])
rec('counts N(T) for Xi, Xi*, P and N_0(T) = (T/2pi) log(T/(2 pi e))', cnt)
# counting functions for the figure (exact zeros up to 100)
out['fig_counting'] = {'xi': [float(g) for g in GAM if g < 100], 'xistar': [float(z) for z in xz if z < 100],
                       'P': [float(t) for t in tP if t < 100]}

# ---------------------------------------------------------------- products on the index side
def logK(mu, x=2 * PI):
    return mp.log(mp.besselk(mu, x))


def tailP(f, V):
    # int_V^oo f(nu) dN_P(nu) with dN_P = Psi'(nu)/pi dnu, Psi' = log(nu/pi) - (pi^2 - 1/12)/nu^2, minus f(V)/2
    dens = lambda nu: (mp.log(nu / PI) - (PI ** 2 - mp.mpf(1) / 12) / nu ** 2) / PI
    return mp.quad(lambda nu: f(nu) * dens(nu), [V, 2 * V, 10 * V, mp.inf]) - f(V) / 2


def tailW(f, V):
    # int_V^oo f(T) dN(T) with the Riemann von Mangoldt density log(T/2pi)/(2pi), minus f(V)/2
    return mp.quad(lambda T: f(T) * mp.log(T / (2 * PI)) / (2 * PI), [V, 2 * V, 10 * V, mp.inf]) - f(V) / 2


rows = []
for th in [mp.mpf('0.5'), mp.mpf(2), mp.mpf(5), mp.mpf('2.89'), mp.mpf(100)]:
    exact = logK(mp.sqrt(th)) - logK(0)
    f = lambda nu: mp.log(1 + th / nu ** 2)
    prod = mp.fsum(f(z) for z in nus) + tailP(f, nus[-1])
    z137 = [z for z in nus if z < 150]
    nick = mp.fsum(f(z) for z in z137) + mp.quad(lambda nu: f(nu) * mp.log(nu / PI) / PI, [z137[-1], mp.inf])
    rows.append([s(th, 4), s(exact, 12), s(prod, 12), s(prod - exact, 3), s(nick - exact, 3)])
rec('index side of P: theta, log(K_sqrt(theta)/K_0)(2pi), sum over nu_n<=1000 + tail, difference, '
    'difference with the 137 zeros below 150 and the plain Weyl tail', rows)
lam = mp.mpf('1.7')
exact17 = logK(lam) - logK(0)
z101 = [z for z in nus if z < 120]
f17 = lambda nu: mp.log(1 + lam ** 2 / nu ** 2)
nick17 = mp.fsum(f17(z) for z in z101) + mp.quad(lambda nu: f17(nu) * mp.log(nu / PI) / PI, [z101[-1], mp.inf])
rec('wave trace identity at lambda = 1.7: log(K_1.7/K_0)(2pi), 101 zeros below 120 + plain Weyl tail',
    [s(exact17, 10), s(nick17, 10), len(z101)])

# Xi*: closed form against the product
def logXs(th):
    r = mp.sqrt(th) / 2
    return mp.log(mp.besselk(mp.mpf(9) / 4 - r, 2 * PI) + mp.besselk(mp.mpf(9) / 4 + r, 2 * PI)) - \
        mp.log(2 * mp.besselk(mp.mpf(9) / 4, 2 * PI))


rows = []
x240 = [z for z in xz if z < 240]
for th in [mp.mpf(10), mp.mpf(60), mp.mpf(200)]:
    f = lambda T: mp.log(1 + th / T ** 2)
    exact = logXs(th)
    prod = mp.fsum(f(z) for z in xz) + tailW(f, xz[-1])
    nick = mp.fsum(f(z) for z in x240) + mp.quad(lambda T: f(T) * mp.log(T / (2 * PI)) / (2 * PI), [x240[-1], mp.inf])
    rows.append([s(th, 4), s(exact, 10), s(prod, 10), s(prod - exact, 3), s(nick, 10)])
rec('Xi*: theta, log(Xi*(i sqrt theta)/Xi*(0)), product over zeros <= 1000 + tail, difference, '
    'product over the zeros below 240 + plain Weyl tail', rows)
rec('Xi*: number of zeros below 240 used in the last column', len(x240))
g = lambda u: (mp.exp(4.5 * u) + mp.exp(-4.5 * u)) * mp.exp(-2 * PI * mp.cosh(2 * u))
Mg60 = mp.quad(lambda u: g(u) * mp.cosh(mp.sqrt(60) * u), [-3, -1, 0, 1, 3]) / mp.quad(g, [-3, -1, 0, 1, 3])
rec('Xi*: M_g(60) from the kernel, exp of the closed form', [s(Mg60, 10), s(mp.exp(logXs(mp.mpf(60))), 10)])


def xi(sv):
    return sv * (sv - 1) / 2 * PI ** (-sv / 2) * mp.gamma(sv / 2) * mp.zeta(sv)


rows = []
for th in [mp.mpf(10), mp.mpf(60), mp.mpf(200)]:
    f = lambda T: mp.log(1 + th / T ** 2)
    exact = mp.log(xi(mp.mpf(1) / 2 + mp.sqrt(th)) / xi(mp.mpf(1) / 2))
    prod = mp.fsum(f(g) for g in GAM) + tailW(f, GAM[-1])
    nick = mp.fsum(f(g) for g in GAM[:300]) + mp.quad(lambda T: f(T) * mp.log(T / (2 * PI)) / (2 * PI), [GAM[299], mp.inf])
    rows.append([s(th, 4), s(exact, 10), s(prod, 10), s(prod - exact, 3), s(nick, 10)])
rec('Xi: theta, log(xi(1/2+sqrt theta)/xi(1/2)), product over 1700 zeros + tail, difference, 300 zeros + plain tail', rows)

# ---------------------------------------------------------------- the index side clock
def kmom(j, x):
    # e^x k_{2j}(x); the integrand is below e^{-400} beyond W, and a finite range keeps mpmath away from exp(-x cosh w)
    # at the huge nodes of an infinite interval (where argument reduction would need about cosh(w) bits)
    W = mp.acosh(1 + 400 / x)
    return mp.quad(lambda w: w ** (2 * j) * mp.exp(-x * (mp.cosh(w) - 1)), [0, W / 16, W / 8, W / 4, W / 2, W])


def a12(x):
    k0, k2, k4 = kmom(0, x), kmom(1, x), kmom(2, x)
    return k2 / (2 * k0), k4 / (24 * k0) - k2 ** 2 / (8 * k0 ** 2)


b = 2 * PI
a1b, a2b = a12(b)
Etau = a1b
Vtau = -2 * a2b
sum2 = mp.fsum(1 / z ** 2 for z in nus) + tailP(lambda nu: 1 / nu ** 2, nus[-1])
sum4 = mp.fsum(1 / z ** 4 for z in nus) + tailP(lambda nu: 1 / nu ** 4, nus[-1])
rec('clock: E tau_oo = K\'\'_0(2pi)/(2 K_0(2pi)), sum nu_n^-2 (zeros + tail)', [s(Etau, 12), s(sum2, 12)])
rec('clock: Var tau_oo, sum nu_n^-4', [s(Vtau, 12), s(sum4, 12)])
R0 = 2 * PI * mp.e ** mp.mpf('2.5')
a1r, a2r = a12(R0)
m0 = a1b - a1r
v0 = 2 * (a2r - a2b)
rec('clock from R_0 = 2 pi e^2.5: R_0, mean, sd', [s(R0, 8), s(m0, 10), s(mp.sqrt(v0), 10)])
rec('clock: nu_1^2, nu_2^2, nu_3^2', [s(z ** 2, 8) for z in nus[:3]])

# ---------------------------------------------------------------- the hitting time identities
# Proposition prop:biane:PT at t = 7
t7 = mp.mpf(7)
P7 = mp.besselk(3.5j, 2 * PI)
dens_T = lambda tau: tau ** mp.mpf(-1.5) * mp.exp(2 * PI - PI * (tau + 1 / tau))
ET = mp.quad(lambda tau: tau ** (mp.mpf(1) / 2 + 3.5j) * dens_T(tau), [0, mp.mpf('0.1'), 1, 10, mp.inf])
rec('PT at t=7: 2 K_{3.5i}(2pi), e^{-2pi} E T^{1/2+3.5i}, normalisation of the density',
    [s(2 * mp.re(P7), 10), s(mp.exp(-2 * PI) * ET, 10), s(mp.quad(dens_T, [0, 1, mp.inf]), 12)])
# Mellin transform of W_1 against 2 (2/pi)^s xi(2s)
fW = lambda x: mp.nsum(lambda n: (PI ** 4 * n ** 4 * x - 3 * PI ** 2 * n ** 2) * mp.exp(-PI ** 2 * n ** 2 * x / 2), [1, mp.inf])
rows = []
for sv in [mp.mpf(1), mp.mpf('0.7'), mp.mpc('0.5', 3)]:
    lhs = mp.quad(lambda x: x ** sv * fW(x), [mp.mpf('0.02'), mp.mpf('0.2'), 1, 5, 40])
    rhs = 2 * (2 / PI) ** sv * xi(2 * sv)
    rows.append([s(sv, 4), s(lhs, 10), s(rhs, 10)])
rec('Mellin of W_1: s, int x^s f_W, 2 (2/pi)^s xi(2s)', rows)
rec('Mellin of W_1: mass of the density below 0.02', s(mp.quad(fW, [mp.mpf('0.005'), mp.mpf('0.02')]), 3))
PhiN = lambda u: mp.nsum(lambda n: (2 * PI ** 2 * n ** 4 * mp.exp(4.5 * u) - 3 * PI * n ** 2 * mp.exp(2.5 * u)) *
                         mp.exp(-PI * n ** 2 * mp.exp(2 * u)), [1, mp.inf])
rows = []
for u in [mp.mpf(0), mp.mpf('0.3'), mp.mpf('-0.4')]:
    flogY = fW(2 / PI * mp.exp(2 * u)) * 4 / PI * mp.exp(2 * u)
    rows.append([s(u, 3), s(PhiN(u), 10), s(mp.exp(u / 2) * flogY / 4, 10)])
rec('Phi_note(u) = Phi(u)/2 against e^{u/2} f_{log Y}(u)/4', rows)
rec('asymptote 2 pi^2 (e^{9u/2}+e^{-9u/2}) e^{-2 pi cosh 2u} at u = 0, and Phi_note(1.9), P kernel at 1.9',
    [s(4 * PI ** 2 * mp.exp(-2 * PI), 6), s(PhiN(mp.mpf('1.9')), 4), s(mp.exp(-2 * PI * mp.cosh(mp.mpf('3.8'))), 4)])
# Levy and Thorin measures of T_x, a = x = 1
tt = mp.mpf('0.8')
lev = mp.exp(-tt / 2) / mp.sqrt(2 * PI * tt ** 3)
lev2 = mp.quad(lambda c: mp.exp(-c * tt) / (PI * mp.sqrt(2 * c - 1)), [mp.mpf(1) / 2, 1, 5, mp.inf]) / tt
inv = mp.sqrt(5) / PI
inv2 = mp.im(mp.sqrt(1 + 2 * mp.mpc(-3, mp.mpf('1e-30')))) / PI
rec('T_x (a=x=1): Levy density at t=0.8 two ways; Thorin mass of (0,3) two ways',
    [s(lev, 10), s(lev2, 10), s(inv, 10), s(inv2, 10)])
# heat kernels
rows = []
for t in [mp.mpf('0.3'), mp.mpf(1)]:
    pH = lambda r: mp.exp(-t / 2) / (2 * PI * t) ** 1.5 * r / mp.sinh(r) * mp.exp(-r ** 2 / (2 * t))
    pS = lambda th: mp.exp(t / 2) / (2 * PI * t) ** 1.5 * mp.fsum((th + 2 * PI * k) / mp.sin(th) *
                                                                   mp.exp(-(th + 2 * PI * k) ** 2 / (2 * t)) for k in range(-12, 13))
    IH = mp.quad(lambda r: pH(r) * 4 * PI * mp.sinh(r) ** 2, [0, 1, 4, 20, mp.inf])
    IS = mp.quad(lambda th: pS(th) * 4 * PI * mp.sin(th) ** 2, [0, PI / 2, PI])
    rows.append([s(t, 3), s(IH, 15), s(IS, 15)])
rec('heat kernels of H^3 and S^3: t, total masses', rows)
# Klein Gordon residual
F = lambda x, w: mp.exp(-mp.exp(x) * mp.cosh(w))
x0, w0 = mp.mpf('0.7'), mp.mpf('0.4')
res = mp.diff(lambda x: F(x, w0), x0, 2) - mp.diff(lambda w: F(x0, w), w0, 2) - mp.exp(2 * x0) * F(x0, w0)
rec('Klein Gordon residual at (0.7, 0.4)', s(res, 3))
# the sine squared images of the two clocks
lam, a, ell = mp.mpf('1.7'), mp.mpf(1), mp.mpf('0.8')
l1 = 2 * mp.log(mp.sinh(lam * a) / (lam * a))
r1 = mp.quad(lambda v: (1 - mp.cos(lam * v)) * 4 / (v * mp.expm1(PI * v / a)), mp.linspace(0, 40, 41) + [mp.inf])
l2 = ell * (mp.sqrt(a ** 2 + lam ** 2) - a)
r2 = mp.quad(lambda v: (1 - mp.cos(lam * v)) * 2 * ell * a / PI * mp.besselk(1, a * v) / v, mp.linspace(0, 60, 61) + [mp.inf])
rec('sine squared images at lambda=1.7, a=1, ell=0.8: Planck two ways, Macdonald two ways',
    [s(l1, 16), s(r1, 16), s(l2, 16), s(r2, 16)])

# ---------------------------------------------------------------- the curvature family
def sk(k, z):
    if k > 0:
        return mp.sin(mp.sqrt(k) * z) / mp.sqrt(k)
    if k < 0:
        return mp.sinh(mp.sqrt(-k) * z) / mp.sqrt(-k)
    return z


rows = []
for k in [mp.mpf(1), mp.mpf(-1), mp.mpf(1) / 3]:
    th, aa, xx = mp.mpf('0.7'), mp.mpf('0.9'), mp.mpf('0.55')
    q = mp.sqrt(2 * th - k + 0j)
    f = lambda x: mp.re(mp.sinh(x * q) / q) / sk(k, x)
    ode = mp.diff(f, xx, 2) / 2 + (mp.diff(lambda x: sk(k, x), xx) / sk(k, xx)) * mp.diff(f, xx) - th * f(xx)
    L = lambda thv: mp.re(mp.sqrt(2 * thv - k + 0j) / mp.sinh(aa * mp.sqrt(2 * thv - k + 0j))) * sk(k, aa)
    L0 = lambda thv: mp.re(mp.sqrt(2 * thv + 0j) * aa / mp.sinh(aa * mp.sqrt(2 * thv + 0j))) if thv != 0 else mp.mpf(1)
    tilt = L0(th - k / 2) / L0(-k / 2)
    cn = lambda n: PI ** 2 * n ** 2 / (2 * aa ** 2) - k / 2
    prodv = mp.exp(-mp.nsum(lambda n: mp.log(1 + th / cn(n)), [1, mp.inf]))
    rows.append([s(k, 4), s(ode, 3), s(L(th), 14), s(tilt, 14), s(prodv, 14)])
rec('curvature family (theta=0.7, a=0.9): kappa, ODE residual at x=0.55, Laplace transform, tilt formula, Thorin product', rows)
# Mellin asymmetry of the tilted doubled clock (a = 1, kappa = 1) against the symmetry at kappa = 0
fW1 = lambda x, k: fW(x) * mp.exp(k * x / 2)
rows = []
for k in [mp.mpf(0), mp.mpf(1)]:
    Z = mp.quad(lambda x: fW1(x, k), [mp.mpf('0.02'), 1, 5, 60])
    M = lambda sv: (PI / 2) ** sv * mp.quad(lambda x: x ** sv * fW1(x, k), [mp.mpf('0.02'), 1, 5, 60]) / Z
    rows.append([s(k, 3), s(M(mp.mpf('0.1')), 10), s(M(mp.mpf('0.4')), 10)])
rec('Mellin symmetry (pi/2)^s E W^s at s=0.1 and s=0.4 (equal iff symmetric under s -> 1/2 - s), kappa = 0, 1', rows)
# Delta h_kappa = -kappa |rho|^2 h_kappa: the Laplacian taken symbolically by sympy for n = 2, 3, 4 and kappa = -1, 0, 1,
# the ratio Delta h/h evaluated at three random rational points with 60 digits
res_sym = []
random.seed(2)
for nn in [2, 3, 4]:
    y = sp.symbols('y1:%d' % (nn + 1))
    for k, fn in [(-1, sp.sinh), (0, lambda z: z), (1, sp.sin)]:
        h = sp.Mul(*[fn(y[i] - y[j]) for i in range(nn) for j in range(i + 1, nn)])
        lap = sum(sp.diff(h, v, 2) for v in y)
        rho2 = sp.Rational(nn * (nn * nn - 1), 3)
        worst = 0
        for _ in range(3):
            pt = {v: sp.Rational(random.randint(-90, 90), 61) for v in y}
            r = sp.N((lap / h).subs(pt), 60)
            worst = max(worst, abs(float(r + k * rho2)))
        res_sym.append([nn, k, str(rho2), worst])
rec('Delta h_kappa / h_kappa + kappa |rho|^2 at random points: n, kappa, |rho|^2, largest absolute value', res_sym)
# Biane's Proposition 2.6: h(theta) = prod sin((theta_j - theta_k)/2) has eigenvalue (n - n^3)/24 for Delta/2
rows = []
random.seed(1)
for nn in [2, 3, 4]:
    th = [mp.mpf(random.randint(1, 60)) / 97 for _ in range(nn - 1)]
    th.append(-mp.fsum(th))
    def hfun(*v):
        return mp.fprod(mp.sin((v[j] - v[k]) / 2) for j in range(nn) for k in range(j + 1, nn))
    lap = mp.fsum(mp.diff(hfun, th, tuple(2 if i == m else 0 for i in range(nn))) for m in range(nn))
    rows.append([nn, s(lap / 2 / hfun(*th), 12), s(mp.mpf(nn - nn ** 3) / 24, 12)])
rec('Biane Prop 2.6 normalisation: n, (Delta h/2)/h in the angles theta, (n - n^3)/24', rows)

out['seconds'] = round(time.time() - t0, 1)
json.dump(out, open('biane_polya.json', 'w'), indent=1)
print('done', out['seconds'])
