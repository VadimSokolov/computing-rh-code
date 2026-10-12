# sqfree_constants.py: the constants of Section sec:sqfree (ch/sqfree.tex), with Table tab:sqfree:zerosum, on the squarefree Dirichlet series
# D(s) = zeta(s)/zeta(2s). Added for the book in October 2026, from Nick Polson's note on the squarefree series
# (misc/reviews/nick-squarefree/), recomputed independently. It computes
#  (1) the constant C of sum_{n<=x} mu(n)^2/n = (6/pi^2)(log x + C) + O(x^{-1/2}) in its Glaisher form 24 log A - gamma - 2 log 2pi,
#      its zeta form gamma - 2 zeta'(2)/zeta(2) and its Xi form 2 + log 2 - (1/2) log pi - gamma/2 + b_1 - 2 b_2, b_a = xi'/xi(a),
#      and the value F(1) = -2 zeta'(2)/zeta(2)^2 of the second answer of the thread;
#  (2) the sum over the zeros sum_gamma (4/(1 + 4 gamma^2) - 24/(9 + 4 gamma^2)) = b_1 - 2 b_2 over the first 60, 1000 and 10000
#      zeros (Arb, acb.zeta_zeros), with the tail from the smooth zero density, and the mean of the centre clock
#      E T = sum 1/gamma^2 = xi''(1/2)/(2 xi(1/2)) both ways;
#  (3) the constant C_K of the Gaussian field Q(i) as gamma_K - 2 zeta_K'(2)/zeta_K(2) and as A_K plus its zero part, the zero part
#      split between zeta and L(s, chi_{-4}), and the L part checked against the zeros of L(s, chi_{-4}) below height 400 (Arb's
#      Hardy Z function of the character, sign changes refined by bisection) with the smooth tail;
#  (4) the factor Gamma_sq(s) = ((2s-1)/(s-1)) 2^s pi^{-(s+1)/2} Gamma((s+1)/2) of zeta(s)/zeta(2s) = Gamma_sq(s) xi(s)/xi(2s) at four
#      points, the Polya form at s = 3/4 with E e^{wX} = int Phi(u) e^{wu} du / xi(1/2) by the trapezoid rule (step 1/400 on [0, 4],
#      Phi of eq:ch2:Phi), the Polya integral at w = 0 and 0.2 + i, and the sign of zeta(a)/zeta(2a) at negative a;
#  (5) the residues of the squarefree series at the first ten zeros: R = zeta(rho/2)/(rho zeta'(rho)) (sharp count),
#      A = Gamma(rho/2) zeta(rho/2)/(2 zeta'(rho)) (exponential smoothing), B = (1 + rho/2) Gamma(rho/2) zeta(rho/2) 2^{-rho/2}/(2 zeta'(rho))
#      (Gamma(2) smoothing), R1 = zeta(rho/2)/((rho - 2) zeta'(rho)) (the alpha = 1 sum), and the constants c_1, c_2 of the double pole
#      at s = -1 of Gamma(s) D(s) x^s in closed form.
# Run from any folder: python3 sqfree_constants.py; writes sqfree_constants.json there (a few minutes, most of it the zeros).
import json
import time

import mpmath as mp
from flint import acb, ctx, dirichlet_char

t0 = time.time()
out = {}


def rec(k, v):
    out[k] = v
    print(k, v, round(time.time() - t0, 1), flush=True)


def s(x, n=15):
    return mp.nstr(x, n)


def cs(z, n=12):
    return [mp.nstr(mp.re(z), n), mp.nstr(mp.im(z), n)]


mp.mp.dps = 40
g = mp.euler
lpi = mp.log(mp.pi)
z2 = mp.zeta(2)
dz2 = mp.zeta(2, derivative=1)
six = 6 / mp.pi ** 2

# (1) the alpha = 1 constant
C_gl = 24 * mp.log(mp.glaisher) - g - 2 * mp.log(2 * mp.pi)
C_z = g - 2 * dz2 / z2
b1 = 1 + g / 2 - mp.log(4 * mp.pi) / 2                       # xi'/xi(1)
b2 = mp.mpf(1) / 2 + 1 - lpi / 2 + mp.digamma(1) / 2 + dz2 / z2  # xi'/xi(2) = 1/2 + 1/(2-1) - log(pi)/2 + psi(1)/2 + zeta'/zeta(2)
C_xi = 2 + mp.log(2) - lpi / 2 - g / 2 + b1 - 2 * b2


def xi(z):
    return z * (z - 1) / 2 * mp.pi ** (-z / 2) * mp.gamma(z / 2) * mp.zeta(z)


b2_num = mp.diff(xi, 2) / xi(2)
b1_num = mp.diff(xi, mp.mpf(1) + mp.mpf('1e-30')) / xi(mp.mpf(1) + mp.mpf('1e-30'))
rec('C: Glaisher form, zeta form, Xi form', [s(C_gl, 25), s(C_z, 25), s(C_xi, 25)])
rec('(6/pi^2) C', s(six * C_gl, 20))
rec('F(1) = -2 zeta\'(2)/zeta(2)^2 (second answer of the thread)', s(-2 * dz2 / z2 ** 2, 20))
rec('int_1^oo Delta(t) t^-2 dt = (6/pi^2)(C - 1)', s(six * (C_gl - 1), 20))
rec('b_1 closed form, numerical derivative; b_2 closed form, numerical derivative', [s(b1, 20), s(b1_num, 20), s(b2, 20), s(b2_num, 20)])
rec('b_1/2 and b_2/2 (the tilted means E_1[T], E_3[T] of the note)', [s(b1 / 2, 12), s(b2 / 2, 12)])
rec('b_1 - 2 b_2 (= G\'/G(1/4) - 6 G\'/G(9/4), the zero part of C)', s(b1 - 2 * b2, 20))
rec('elementary part 2 + log 2 - log(pi)/2 - gamma/2', s(2 + mp.log(2) - lpi / 2 - g / 2, 20))

# (2) sums over the zeros of zeta
ctx.prec = 80
NZ = 10000
GAM = [mp.mpf(z.imag.mid().str(25, radius=False)) for z in acb.zeta_zeros(1, NZ + 1)]
rec('zeros computed', [NZ + 1, s(GAM[0], 15), s(GAM[NZ], 15)])


def fz(t):
    return 4 / (1 + 4 * t * t) - 24 / (9 + 4 * t * t)


def nsmooth(T, q=1, c0=mp.mpf(7) / 8):
    return T / (2 * mp.pi) * mp.log(q * T / (2 * mp.pi * mp.e)) + c0


def tail(f, T, k, q=1, c0=mp.mpf(7) / 8):
    # sum_{gamma > T} f(gamma) ~ int_T^oo f(t) (1/2pi) log(q t/2pi) dt + f(T) (N_smooth(T) - k), N(T) = k
    integ = mp.quad(lambda t: f(t) * mp.log(q * t / (2 * mp.pi)) / (2 * mp.pi), [T, 2 * T, 10 * T, mp.inf])
    return integ + f(T) * (nsmooth(T, q, c0) - k)


zs = {}
for K in (60, 1000, 10000):
    part = mp.fsum(fz(t) for t in GAM[:K])
    T = (GAM[K - 1] + GAM[K]) / 2
    tl = tail(fz, T, K)
    zs[K] = {'partial sum': s(part, 12), 'smooth tail': s(tl, 8), 'total': s(part + tl, 12),
             'total minus (b_1 - 2 b_2)': s(part + tl - (b1 - 2 * b2), 3)}
rec('sum over zeros of 4/(1+4g^2) - 24/(9+4g^2)', zs)
inv2_60 = mp.fsum(1 / t ** 2 for t in GAM[:60])
ET = mp.diff(xi, mp.mpf(1) / 2, 2) / (2 * xi(mp.mpf(1) / 2))
rec('E T = xi\'\'(1/2)/(2 xi(1/2))', s(ET, 25))
inv2 = mp.fsum(1 / t ** 2 for t in GAM[:NZ])
T = (GAM[NZ - 1] + GAM[NZ]) / 2
rec('sum 1/gamma^2 over 10000 zeros plus smooth tail, minus E T', [s(inv2, 15), s(inv2 + tail(lambda t: 1 / t ** 2, T, NZ) - ET, 3)])
rec('the note\'s tail estimate -5 sum_{k>60} 1/gamma_k^2 = -5 (E T - sum_{k<=60} 1/gamma_k^2)', s(-5 * (ET - inv2_60), 6))
# empirical check of the constant 7/8 in N(T) from the zeros: mean of (k - 1/2) - N_smooth(gamma_k) without the constant
c0z = mp.fsum((k + 1 - mp.mpf(1) / 2) - nsmooth(GAM[k], 1, 0) for k in range(100, NZ)) / (NZ - 100)
rec('mean of N(T) - (T/2pi) log(T/(2 pi e)) at the zeros, against 7/8', s(c0z, 6))

# (3) the Gaussian field Q(i)
chi4 = [0, 1, 0, -1]
L1 = mp.dirichlet(1, chi4)
dL1 = mp.dirichlet(1, chi4, derivative=1)
L2 = mp.dirichlet(2, chi4)
dL2 = mp.dirichlet(2, chi4, derivative=1)
LL1_closed = g + 2 * mp.log(2) + 3 * lpi - 4 * mp.log(mp.gamma(mp.mpf(1) / 4))
gK = g + dL1 / L1
CK = gK - 2 * (dz2 / z2 + dL2 / L2)
AK = 4 - g - lpi                                   # 2 + (1/2) log 4 + (2 - gamma - log 2pi)
xchi = lambda a, LL: mp.log(4 / mp.pi) / 2 + mp.digamma((a + 1) / mp.mpf(2)) / 2 + LL   # xi_chi'/xi_chi, xi_chi = (4/pi)^{(s+1)/2} Gamma((s+1)/2) L
Lpart = xchi(1, dL1 / L1) - 2 * xchi(2, dL2 / L2)
rec('Q(i): L(1) - pi/4, L(2) - Catalan, L\'/L(1) and its closed form', [s(L1 - mp.pi / 4, 3), s(L2 - mp.catalan, 3), s(dL1 / L1, 20), s(LL1_closed, 20)])
rec('Q(i): gamma_K, zeta_K\'/zeta_K(2), C_K = gamma_K - 2 zeta_K\'(2)/zeta_K(2)', [s(gK, 20), s(dz2 / z2 + dL2 / L2, 20), s(CK, 22)])
rec('Q(i): residue over zeta_K(2) = 3/(2 pi Catalan), and (pi/4)/(zeta(2) L(2))', [s(3 / (2 * mp.pi * mp.catalan), 15), s(L1 / (z2 * L2), 15)])
rec('Q(i): A_K = 4 - gamma - log pi; zeta part b_1 - 2 b_2; L part; A_K + both - C_K',
    [s(AK, 20), s(b1 - 2 * b2, 15), s(Lpart, 15), s(AK + b1 - 2 * b2 + Lpart - CK, 3)])
ctx.prec = 64
chi = dirichlet_char(4, 3)
rec('character mod 4 index 3: values at 1, 3; primitive; parity', [str(chi(1)), str(chi(3)), chi.is_primitive(), chi.parity()])


def Z(t):
    return float(chi.hardy_z(acb(t)).real.mid())


TL = 400.0
step = 0.02
zeros_L = []
t_prev, z_prev = 0.5, Z(0.5)
t = t_prev
while t < TL:
    t = t_prev + step
    z = Z(t)
    if z == 0 or (z > 0) != (z_prev > 0):
        a, b, fa = t_prev, t, z_prev
        for _ in range(45):
            m = (a + b) / 2
            fm = Z(m)
            if (fm > 0) == (fa > 0):
                a, fa = m, fm
            else:
                b = m
        zeros_L.append(mp.mpf((a + b) / 2))
    t_prev, z_prev = t, z
nL = len(zeros_L)
c0L = mp.fsum((k + 1 - mp.mpf(1) / 2) - nsmooth(zeros_L[k], 4, 0) for k in range(20, nL)) / (nL - 20)
partL = mp.fsum(fz(t) for t in zeros_L)
tlL = tail(fz, mp.mpf(TL), nL, q=4, c0=c0L)
rec('L(s, chi_-4): zeros below 400, first three, empirical constant of N(T)', [nL, [s(x, 10) for x in zeros_L[:3]], s(c0L, 5)])
rec('L part from the zeros below 400 plus smooth tail, minus the closed form', [s(partL, 12), s(tlL, 8), s(partL + tlL - Lpart, 3)])

# (4) Lemma Gamma_sq and the Polya form


def D(z):
    return mp.zeta(z) / mp.zeta(2 * z)


def Gsq(z):
    return (2 * z - 1) / (z - 1) * 2 ** z * mp.pi ** (-(z + 1) / 2) * mp.gamma((z + 1) / 2)


lem = {}
for z in (mp.mpf(3), mp.mpc('0.3', '2.1'), mp.mpf('0.75'), mp.mpf('-1.5')):
    lem[str(z)] = s(abs(D(z) / (Gsq(z) * xi(z) / xi(2 * z)) - 1), 3)
rec('|D(s)/(Gamma_sq(s) xi(s)/xi(2s)) - 1|', lem)
rec('s = 3/4: zeta(3/4)/zeta(3/2), Gamma_sq(3/4), xi(3/4)/xi(3/2)', [s(D(mp.mpf('0.75')), 12), s(Gsq(mp.mpf('0.75')), 12), s(xi(mp.mpf('0.75')) / xi(mp.mpf('1.5')), 12)])


def Phi(u):
    q = mp.exp(-mp.pi * mp.exp(2 * u))
    a9, a5 = mp.exp(9 * u / 2), mp.exp(5 * u / 2)
    tot, n = mp.mpf(0), 1
    while True:
        term = (2 * mp.pi ** 2 * n ** 4 * a9 - 3 * mp.pi * n ** 2 * a5) * q ** (n * n)
        tot += term
        if n > 2 and abs(term) < mp.eps * abs(tot) * mp.mpf(10) ** -10:
            break
        n += 1
    return 2 * tot


h = mp.mpf(1) / 400
nodes = [(j * h, Phi(j * h)) for j in range(0, 1601)]


def polya(w):
    # int_R Phi(u) e^{wu} du = 2 int_0^oo Phi(u) cosh(wu) du, trapezoid rule
    return 2 * h * (mp.fsum(p * mp.cosh(w * u) for u, p in nodes) - nodes[0][1] / 2)


xih = xi(mp.mpf(1) / 2)
pol = {}
for w in (mp.mpf(0), mp.mpc('0.2', '1'), mp.mpf('0.25'), mp.mpf(1)):
    pol[str(w)] = s(abs(polya(w) / xi(mp.mpf(1) / 2 + w) - 1), 3)
rec('|int Phi(u) e^{wu} du / xi(1/2 + w) - 1|', pol)
EX1, EX2 = polya(mp.mpf('0.25')) / xih, polya(mp.mpf(1)) / xih
rec('s = 3/4: E e^{X/4}, E e^{X}, Gamma_sq(3/4) E e^{X/4}/E e^{X} - zeta(3/4)/zeta(3/2)',
    [s(EX1, 15), s(EX2, 15), s(Gsq(mp.mpf('0.75')) * EX1 / EX2 - D(mp.mpf('0.75')), 3)])
sg = []
for a in ('-0.5', '-1.5', '-2.5', '-3.5', '-4.5', '-5.5', '-6.5'):
    a = mp.mpf(a)
    sg.append({'alpha': s(a, 3), 'zeta(a)/zeta(2a)': s(D(a), 8), '(2a-1)/(a-1)': s((2 * a - 1) / (a - 1), 6),
               'Gamma((a+1)/2)': s(mp.gamma((a + 1) / 2), 6), 'xi(a)/xi(2a)': s(xi(a) / xi(2 * a), 6)})
rec('signs at negative alpha', sg)

# (5) residues at the zeros and the trivial pole at s = -1
mp.mp.dps = 30
res = []
for k in range(10):
    rho = mp.mpc(mp.mpf(1) / 2, GAM[k])
    zh, dz = mp.zeta(rho / 2), mp.zeta(rho, derivative=1)
    R = zh / (rho * dz)
    A = mp.gamma(rho / 2) * zh / (2 * dz)
    B = (1 + rho / 2) * mp.gamma(rho / 2) * zh * 2 ** (-rho / 2) / (2 * dz)
    R1 = zh / ((rho - 2) * dz)
    res.append({'k': k + 1, 'gamma': s(GAM[k], 12), 'R': cs(R, 6), '|R|': s(abs(R), 6), '|A|': s(abs(A), 4), '|B|': s(abs(B), 4),
                '|R1|': s(abs(R1), 4), 'zeta(rho/2)': cs(zh, 6)})
rec('residues at the first ten zeros', res)
c1 = -mp.zeta(-1) / (2 * mp.zeta(-2, derivative=1))
c2 = c1 * (mp.digamma(2) + mp.zeta(-1, derivative=1) / mp.zeta(-1) - mp.zeta(-2, derivative=2) / mp.zeta(-2, derivative=1))
rec('c_1 = -zeta(-1)/(2 zeta\'(-2)), c_2 (double pole of Gamma(s) D(s) x^s at s = -1)', [s(c1, 15), s(c2, 15)])
rec('residue coefficient of (1+s) Gamma(s) D(s) (X/2)^s at s = -1: -zeta(-1)/zeta\'(-2) = 2 c_1', s(2 * c1, 15))
out['seconds'] = time.time() - t0
json.dump(out, open('sqfree_constants.json', 'w'), indent=1)
print('done', round(time.time() - t0, 1))
