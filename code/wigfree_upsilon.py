# wigfree_upsilon.py: the Upsilon image and the free counterpart of the centre clock, for the subsections on the Upsilon
# transform and on the free clock in the free probability section (misc/reviews/nick-wigner/merge/wigfree.tex, labels
# sec:wigfree:ups and sec:wigfree:clock). Added for the book in October 2026, from the first author's note on the semicircle law
# and the Thorin framework (misc/reviews/nick-wigner/source-v7/, version 7, Sections 3 and 4 and Table 1), recomputed and extended.
# Notation of the book: G(u) = xi(1/2 + sqrt u)/xi(1/2), W(t) = sum_j exp(-a_j t), a_j = gamma_j^2, p_n = sum_j a_j^{-n}.
# Part A, the table: psi(s) = int_0^oo log G(s x) e^{-x} dx, which is the Laplace exponent of the Upsilon image of the clock and
# minus the free cumulant transform of its free counterpart at -s, by quadrature from xi alone (30 digits), against
# sum_j phi(a_j/s), phi(z) = e^z E_1(z), over the first 300 and 1700 zeros, with the first order tail s(log(T/2 pi) + 1)/(2 pi T)
# of the note, T the last ordinate used, and with the full smooth tail int_T^oo phi(g^2/s) dN_0(g), dN_0 = log(g/2 pi) dg/(2 pi),
# plus the boundary term phi(T^2/s)(N_0(T) - N), N_0(T) = (T/2 pi) log(T/2 pi e) + 7/8. At s = 10^3 and 10^4 the same comparison
# measures the accuracy of the expansion used in Part C.
# Part B: the Thorin measure of the Upsilon image has distribution function x -> W(1/x); in the chart x = theta^2 this is
# W(theta^{-2}) = sum_j exp(-gamma_j^2/theta^2), recorded with the zero count N(theta) on a grid for the figure.
# Part C: densities of three laws with the Levy measure W(t) dt/t: the clock T (fixed Talbot inversion of 1/G, mpmath, 25
# digits), its Upsilon image (fixed Talbot inversion of exp(-psi(u)), double precision, with psi(u) = sum over the 1700 zeros of
# phi(a_j/u) plus u tau_1 - u^2 tau_2 + 2 u^3 tau_3, tau_k the smooth tail of sum a_j^{-k} with the boundary term), and the
# free counterpart Lambda(T), whose Cauchy transform w = G(x + i0) solves 1 + x s = psi(s), s = -w, so that the density is
# Im(s)/pi on the support [c, oo), c = max_{s > 0} (psi(s) - 1)/s. The densities are checked against the cumulants: T has mean
# p_1 and variance p_2, the Upsilon image p_1 and 2 p_2, and Lambda(T) p_1 and p_2, with third central moment 2 p_3 for T and
# Lambda(T). Reads zeros_1700.txt; run from its folder: python3 wigfree_upsilon.py; writes wigfree_upsilon.json (a few minutes).
import json
import time

import mpmath as mp
import numpy as np

trap = getattr(np, 'trapezoid', None) or np.trapz

t0 = time.time()
out = {}


def rec(k, v):
    out[k] = v
    print(k, v, round(time.time() - t0, 1), flush=True)


def s_(x, n=12):
    return mp.nstr(x, n)


GAM = [mp.mpf(l.split()[-1]) for l in open('zeros_1700.txt') if l.strip() and not l.startswith('#')]
NZ = len(GAM)
rec('zeros read, last ordinate', [NZ, s_(GAM[-1], 15)])

# ---------------- Part A ----------------
mp.mp.dps = 30
HALF = mp.mpf(1) / 2


def xi(s):
    return s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)


XI0 = xi(HALF)


def logG(u):
    if u == 0:
        return mp.mpf(0)
    return mp.log(xi(HALF + mp.sqrt(u)) / XI0)


def phi(z):
    return mp.exp(z) * mp.e1(z)


def N0(T):
    return T / (2 * mp.pi) * mp.log(T / (2 * mp.pi * mp.e)) + mp.mpf(7) / 8


def tail_first(s, T):
    return s * (mp.log(T / (2 * mp.pi)) + 1) / (2 * mp.pi * T)


def tail_full(s, T, N):
    integral = mp.quad(lambda g: phi(g * g / s) * mp.log(g / (2 * mp.pi)) / (2 * mp.pi), [T, 2 * T, 8 * T, mp.inf])
    return integral + phi(T * T / s) * (N0(T) - N)


tabA = {}
for sv in (1, 10, 100, 1000, 10000):
    s = mp.mpf(sv)
    pts = sorted(set([mp.mpf(0)] + [k / s for k in (1, 10, 100, 1000) if k / s < 1] + [mp.mpf(v) for v in (1, 5, 20, 60, 120)]))
    I = mp.quad(lambda x: logG(s * x) * mp.exp(-x), pts)
    terms = [phi(g * g / s) for g in GAM]
    Z300 = mp.fsum(terms[:300])
    Z1700 = mp.fsum(terms)
    T300, T1700 = GAM[299], GAM[-1]
    row = {'integral from xi': s_(I, 15), '300 zeros': s_(Z300, 12), '300 zeros + first order tail': s_(Z300 + tail_first(s, T300), 12),
           '300 zeros + full tail': s_(Z300 + tail_full(s, T300, 300), 12),
           '1700 zeros + first order tail': s_(Z1700 + tail_first(s, T1700), 12),
           '1700 zeros + full tail': s_(Z1700 + tail_full(s, T1700, NZ), 15)}
    row['differences from the integral: 300, 300+first, 300+full, 1700+first, 1700+full'] = [
        s_(Z300 - I, 3), s_(Z300 + tail_first(s, T300) - I, 3), s_(Z300 + tail_full(s, T300, 300) - I, 3),
        s_(Z1700 + tail_first(s, T1700) - I, 3), s_(Z1700 + tail_full(s, T1700, NZ) - I, 3)]
    tabA[sv] = row
    rec('free transform at s = %d' % sv, row)

# ---------------- Part B ----------------
gam = np.array([float(g) for g in GAM])
A = gam ** 2
th = np.linspace(0.5, 100, 400)
Wth = np.array([np.exp(-A / t ** 2).sum() for t in th])
Nth = np.array([(gam <= t).sum() for t in th])
out['thorin of the Upsilon image, theta grid'] = {'theta': [round(float(t), 4) for t in th], 'W(1/theta^2)': [float('%.8g' % w) for w in Wth],
                                                  'N(theta)': [int(n) for n in Nth]}
ratio = {}
for t in (50.0, 100.0, 200.0, 400.0):
    ratio[t] = [float('%.6g' % np.exp(-A / t ** 2).sum()), int((gam <= t).sum())]
rec('W(1/theta^2) and N(theta) at theta = 50, 100, 200, 400 (Gaussian smoothing of the count; ratio tends to sqrt(pi)/2)', ratio)

# ---------------- Part C ----------------
# p_1, p_2, p_3 from the Taylor coefficients (as in wigfree_jacobi.py, 200 digits suffice here)
mp.mp.dps = 200
hh = mp.mpf(1) / 200


def phi_nodes(h, V):
    nodes = []
    for j in range(int(mp.nint(V / h)) + 1):
        v = j * h
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
        nodes.append((v, 2 * tot))
    return nodes


M = [mp.mpf(0)] * 5
for v, ph in phi_nodes(hh, mp.mpf(3)):
    w = 2 * hh * ph * (mp.mpf(1) / 2 if v == 0 else 1)
    for k in range(5):
        M[k] += w * v ** (2 * k)
a = [M[k] / (mp.factorial(2 * k) * M[0]) for k in range(5)]
psi_ = [mp.mpf(0)] * 5
for n in range(1, 5):
    psi_[n] = (n * a[n] - mp.fsum(j * psi_[j] * a[n - j] for j in range(1, n))) / n
P = [None] + [float((-1) ** (n + 1) * n * psi_[n]) for n in range(1, 5)]
rec('p_1, p_2, p_3 (Taylor side)', ['%.12g' % P[1], '%.12g' % P[2], '%.12g' % P[3]])
mp.mp.dps = 25

# tails tau_k = sum_{gamma > T} gamma^{-2k} by the smooth density with the boundary term, T = gamma_1700
T = float(GAM[-1])
N0T = T / (2 * np.pi) * np.log(T / (2 * np.pi * np.e)) + 7 / 8


def tau(k):
    m = 2 * k - 1
    integral = T ** (-m) / (2 * np.pi * m) * (np.log(T / (2 * np.pi)) + 1 / m)
    return integral + T ** (-2 * k) * (N0T - NZ)


TAU = [None, tau(1), tau(2), tau(3)]
rec('tail sums tau_1, tau_2, tau_3 beyond gamma_1700', ['%.6g' % TAU[1], '%.6g' % TAU[2], '%.6g' % TAU[3]])


def phiv(z):
    # e^z E_1(z), principal branch, elementwise, in double precision: the asymptotic series for |z| >= 30 (error below 10^-12),
    # the power series of E_1 for |z| < 2 and for |arg z| >= 3 pi/4 (little cancellation there), and the continued
    # fraction e^z E_1(z) = 1/(z + 1 - 1/(z + 3 - 4/(z + 5 - ...))), evaluated backwards from depth 160, for the rest
    z = np.asarray(z, dtype=complex)
    res = np.empty_like(z)
    big = np.abs(z) >= 30
    ser = (~big) & ((np.abs(z) < 2) | ((z.real < 0) & (np.abs(z.imag) <= -z.real)))
    cf = (~big) & (~ser)
    zb = z[big]
    acc = np.zeros_like(zb)
    term = 1 / zb
    for n in range(16):
        acc += term
        term = term * (-(n + 1)) / zb
    res[big] = acc
    zs = z[ser]
    tot = np.zeros_like(zs)
    term = np.ones_like(zs)
    for k in range(1, 121):
        term = term * (-zs) / k
        tot += term / k
    res[ser] = np.exp(zs) * (-np.euler_gamma - np.log(zs) - tot)
    zc = z[cf]
    tt = zc + 2 * 160 + 1
    for n in range(159, -1, -1):
        tt = zc + 2 * n + 1 - (n + 1) ** 2 / tt
    res[cf] = 1 / tt
    return res


ztest = np.array([0.5 + 0.5j, 3 - 4j, 20 - 20j, -10 + 25j, -25 - 3j, 29 - 1j, 2 + 0j, -5 - 0.1j, 40 - 10j, -35 - 1j, 2.5 - 25j, -2.1 - 2.2j])
zref = np.array([complex(mp.exp(mp.mpc(z.real, z.imag)) * mp.e1(mp.mpc(z.real, z.imag))) for z in ztest])
rec('scaled exponential integral against mpmath at 12 test points, largest relative error', '%.2g' % np.max(np.abs(phiv(ztest) / zref - 1)))


def psi(u):
    # Laplace exponent of the Upsilon image at complex u (array), off (-oo, 0]
    u = np.atleast_1d(np.asarray(u, dtype=complex))
    z = A[None, :] / u[:, None]
    return phiv(z).sum(axis=1) + u * TAU[1] - u ** 2 * TAU[2] + 2 * u ** 3 * TAU[3]


def dpsi(u):
    u = np.atleast_1d(np.asarray(u, dtype=complex))
    z = A[None, :] / u[:, None]
    return ((phiv(z) - 1 / z) * (-z / u[:, None])).sum(axis=1) + TAU[1] - 2 * u * TAU[2] + 6 * u ** 2 * TAU[3]


# accuracy of psi against Part A at real s
chk = {}
for sv in (1, 10, 100, 1000, 10000):
    chk[sv] = '%.3g' % (psi(float(sv))[0].real - float(mp.mpf(tabA[sv]['integral from xi'])))
rec('psi (1700 zeros + three term tail) minus the integral, s = 1 .. 10^4', chk)

# C1: density of T by mpmath Talbot inversion of 1/G
xT = np.round(np.arange(0.002, 0.1500001, 0.0005), 6)
fT = []
for x in xT:
    fT.append(float(mp.invertlaplace(lambda u: XI0 / xi(HALF + mp.sqrt(u)), mp.mpf(str(x)), method='talbot')))
fT = np.array(fT)
rec('density of T computed at %d points' % len(xT), round(time.time() - t0, 1))


# C2: density of the Upsilon image by the fixed Talbot method in double precision
def talbot_upsilon(x, Mt=18):
    r = 2 * Mt / (5 * x)
    k = np.arange(1, Mt)
    thk = k * np.pi / Mt
    cot = np.cos(thk) / np.sin(thk)
    S = r * thk * (cot + 1j)
    sig = thk + (thk * cot - 1) * cot
    F0 = np.exp(-psi(r)[0].real)
    FS = np.exp(-psi(S))
    return r / Mt * (0.5 * F0 * np.exp(r * x) + np.sum((np.exp(x * S) * FS * (1 + 1j * sig)).real))


# Along the left part of the contour 1/E e^{-u Upsilon(T)} grows like exp(0.7 sqrt|u|), so at small x the double precision sum
# cancels badly; below x = 0.0125 the inversion runs in mpmath (Talbot, 30 digits and more inside), with the same zero sum.
AM = [g * g for g in GAM]
TAUm = [None] + [mp.mpf(float(v)) for v in TAU[1:]]


def F_upsilon_mp(u):
    tot = mp.fsum(mp.exp(a / u) * mp.e1(a / u) for a in AM)
    return mp.exp(-(tot + u * TAUm[1] - u ** 2 * TAUm[2] + 2 * u ** 3 * TAUm[3]))


xUs = np.round(np.arange(0.004, 0.01249, 0.0005), 6)
mp.mp.dps = 30
fUs = np.array([float(mp.invertlaplace(F_upsilon_mp, mp.mpf(str(x)), method='talbot')) for x in xUs])
over = {}
for x in (0.0125, 0.015, 0.02, 0.03):
    over[x] = ['%.10g' % float(mp.invertlaplace(F_upsilon_mp, mp.mpf(str(x)), method='talbot')), '%.10g' % talbot_upsilon(x)]
rec('Upsilon image: density in mpmath and in double precision at x = 0.0125, 0.015, 0.02, 0.03', over)
mp.mp.dps = 25
xUb = np.round(np.concatenate([np.arange(0.0125, 0.15, 0.0005), np.arange(0.15, 0.6, 0.002), np.arange(0.6, 1.2001, 0.01)]), 6)
fUb = np.array([talbot_upsilon(x) for x in xUb])
fUb16 = np.array([talbot_upsilon(x, 14) for x in xUb])
xU = np.concatenate([xUs, xUb])
fU = np.concatenate([fUs, fUb])
fU16 = np.concatenate([fUs, fUb16])
rec('Upsilon image: density at x = 0.004, 0.006, 0.008, 0.010 (mpmath)', ['%.4g' % fUs[i] for i in (0, 4, 8, 12)])
rec('Upsilon image: largest change of the density between 13 and 17 Talbot nodes, relative to its maximum', '%.2g' % (np.abs(fU - fU16).max() / fU.max()))
# the right tail: the Levy density of the Upsilon image is (2/sqrt x) sum_j gamma_j K_1(2 gamma_j sqrt x), stretched exponential
mp.mp.dps = 20
tailchk = {}
for xv in (0.2, 0.4, 0.7, 1.0):
    lev = 2 / mp.sqrt(xv) * mp.fsum(g * mp.besselk(1, 2 * g * mp.sqrt(xv)) for g in GAM[:40])
    i = int(np.argmin(np.abs(xU - xv)))
    tailchk[xv] = ['%.6g' % fU[i], '%.6g' % float(lev), '%.4g' % (fU[i] / float(lev))]
rec('Upsilon image, right tail: density, Levy density (2/sqrt x) sum gamma_j K_1(2 gamma_j sqrt x), ratio, at x = 0.2, 0.4, 0.7, 1', tailchk)
mp.mp.dps = 25

# C3: density of Lambda(T) by subordination
def H(s, x):
    return 1 + x * s - psi(s)[0]


# left edge: q(s) = s psi'(s) - psi(s) + 1 decreases from 1; its root s* gives c = psi'(s*)
lo, hi = 1.0, 10.0
while (hi * dpsi(hi)[0] - psi(hi)[0] + 1).real > 0:
    hi *= 2
for _ in range(200):
    mid = (lo + hi) / 2
    if (mid * dpsi(mid)[0] - psi(mid)[0] + 1).real > 0:
        lo = mid
    else:
        hi = mid
sstar = (lo + hi) / 2
c = dpsi(sstar)[0].real
h2 = 1e-3 * sstar
psi2 = ((psi(sstar + h2)[0] - 2 * psi(sstar)[0] + psi(sstar - h2)[0]) / h2 ** 2).real
rec('free clock: left edge c of the support, s* and psi\'\'(s*)', ['%.10g' % c, '%.10g' % sstar, '%.6g' % psi2])
K = 500
xmax = 0.15
xL = c + (xmax - c) * (np.arange(1, K + 1) / K) ** 2
sv = []
s = sstar + 1j * np.sqrt(2 * (xL[0] - c) * sstar / abs(psi2))
for x in xL:
    for it in range(60):
        Hs = H(s, x)
        ds = Hs / (x - dpsi(s)[0])
        s = s - ds
        if abs(ds) < 1e-14 * abs(s):
            break
    if s.imag <= 0:
        s = complex(s.real, abs(s.imag))
    sv.append(s)
sv = np.array(sv)
fL = sv.imag / np.pi
resid = max(abs(H(s, x)) for s, x in zip(sv, xL))
rec('free clock: largest residual |1 + x s - psi(s)| on the grid, least Im s', ['%.2g' % resid, '%.3g' % sv.imag.min()])


def moments(x, f, edge=None):
    # trapezoid moments; for the free clock add the edge piece int_c^{x_1} f ~ (2/3) f(x_1) (x_1 - c)
    m0 = trap(f, x)
    m1 = trap(x * f, x)
    if edge is not None:
        m0 += 2 / 3 * f[0] * (x[0] - edge)
        m1 += 2 / 3 * f[0] * (x[0] - edge) * (x[0] + edge) / 2
    mean = m1 / m0
    var = trap((x - mean) ** 2 * f, x) / m0
    k3 = trap((x - mean) ** 3 * f, x) / m0
    return m0, mean, var, k3


for name, x, f, e in (('T', xT, fT, None), ('Upsilon(T)', xU, fU, None), ('Lambda(T)', xL, fL, c)):
    m0, mean, var, k3 = moments(x, f, e)
    i = int(np.argmax(f))
    rec('%s: mass, mean, variance, third central moment, mode, density at the mode' % name,
        ['%.8f' % m0, '%.8g' % mean, '%.6g' % var, '%.5g' % k3, '%.5g' % x[i], '%.6g' % f[i]])
rec('targets: p_1, p_2, 2 p_2, 2 p_3, 12 p_3', ['%.8g' % P[1], '%.6g' % P[2], '%.6g' % (2 * P[2]), '%.5g' % (2 * P[3]), '%.5g' % (12 * P[3])])
# unimodality on the grids: the number of sign changes of the discrete derivative
for name, f in (('T', fT), ('Upsilon(T)', fU), ('Lambda(T)', fL)):
    f = f[f > 1e-8 * f.max()]
    d = np.sign(np.diff(f))
    d = d[d != 0]
    rec('%s: sign changes of the discrete derivative where the density exceeds 1e-8 of its maximum' % name, int((d[1:] != d[:-1]).sum()))
out['densities'] = {'T': {'x': xT.tolist(), 'f': [float('%.8g' % v) for v in fT]},
                    'Upsilon(T)': {'x': xU.tolist(), 'f': [float('%.8g' % v) for v in fU]},
                    'Lambda(T)': {'x': [float('%.8g' % v) for v in xL], 'f': [float('%.8g' % v) for v in fL]}}
out['seconds'] = round(time.time() - t0, 1)
json.dump(out, open('wigfree_upsilon.json', 'w'), indent=1)
print('done', out['seconds'])
