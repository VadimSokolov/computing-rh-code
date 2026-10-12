# maxexp_checks.py: the numbers of Section sec:ch16:max (ch/ch16.tex, Section 8.5.4) that code/maxlaw.py does not compute, and
# Table tab:maxexp:checks. Added for the book in October 2026, from the first author's note on maxima of exponentials
# (misc/reviews/nick-maxexp/, second version, Sections 1 to 5 and its table of checks), recomputed independently.
# Notation of the book: L(s, Delta) = sum tau(n) n^{-s-11/2}, xi(s, Delta) = (2 pi)^{-s-11/2} Gamma(s+11/2) L(s, Delta);
# F(y) = prod (1 - e^{-2 pi n y})^24 = P(M* <= y), F_1(y) = prod (1 - e^{-2 pi n y}) = P(M_1 <= y), b_n = pi (n^2 - 1)/12.
#  1. the density of M* in Ramanujan's form 2 pi F(y)(1 - E_2(iy)) and its functional equation (quasi modularity of E_2);
#  2. the mixed Poisson probabilities P_k = P(N(2 pi M*) <= k - 1) = L(k - 11/2, Delta), k = 1, ..., 11, by quadrature against
#     the density of M*, the five identities (2 pi)^{-k} (k-1)! P_k = (2 pi)^{k-12} (11-k)! P_{12-k}, Hecke's values, the gamma
#     mixture form E F(G_a/2 pi) at a = 11/2 (s = 0), and L(17/2, Delta) = E Q(14, 2 pi M*) and L(29/2, Delta) = E Q(20, 2 pi M*)
#     against Hecke's formula and the Dirichlet series;
#  3. Riemann's incomplete gamma formula for pi^{-s/2} Gamma(s/2) zeta(s) at s = 2 and 1/2 + 5i;
#  4. eta: L(s, chi_12) = E Q(s/2, pi M_1/12) at s = 0.6 + 5i and 1/2, the symmetric form xi(s, chi_12) = A(s) + A(1 - s), and
#     E exp(-pi M_1/12) = L(2, chi_12) = pi^2/(6 sqrt 3);
#  5. the conductor 11 curve: L(1, E) = E exp(-2 pi M_11) against the series 2 sum a_n e^{-2 pi n/sqrt 11}/n;
#  6. theta_4: the triple product against the theta series, E exp(-3 M_4), E K and the moment formula at s = 3, 1/2 + 5i;
#  7. M_1: the weights p_k of the Renyi type identity, -psi_1'/psi_1 at 2 three ways, the Laplace exponent at 2 as a Levy
#     integral, the two term expansion of the Levy density at 0 and its tail;
#  8. M*: -log E exp(-theta M*) against 2 sqrt(2 pi theta) - (25/4) log theta + (21/4) log 2 pi + (1/2) log 2, and the
#     cumulants kappa_1, ..., kappa_8 of M* tilted by exp(-theta M*) at theta = 0.5, 2, 8, 30 (all positive if M* is ID);
#  9. finite maxima: for random families of 2 to 6 rates (rational, exact grouping of equal subset sums) the poles and zeros
#     of the Laplace transform and the canonical function k(x) = sum_poles e^{-px} - sum_zeros e^{-zx} on a grid (the law is
#     infinitely divisible exactly when k >= 0), and the same for the truncations of the eta law (rates 2 pi n, n <= N, once)
#     and of the Delta law (rates 2 pi n, n <= N, 24 times each), which tend to M_1 and M* as N grows.
# mpmath only. Run from any folder: python3 maxexp_checks.py; writes maxexp_checks.json there.
import json, time, random
from fractions import Fraction
import mpmath as mp

mp.mp.dps = 40
t0 = time.time()
out = {}
def put(key, val):
    out[key] = val
    print(key, val, round(time.time() - t0, 1), flush=True)
S = lambda x, n=25: mp.nstr(x, n)

NT = 400
sig = [0] * (NT + 1)
for d in range(1, NT + 1):
    for m in range(d, NT + 1, d):
        sig[m] += d
tau = [0] * (NT + 1); tau[1] = 1
for n in range(2, NT + 1):
    tau[n] = -24 * sum(sig[m] * tau[n - m] for m in range(1, n)) // (n - 1)
assert tau[2] == -24 and tau[3] == 252 and tau[11] == 534612
H = mp.mpf(11) / 2
TINY = mp.mpf(10) ** (-mp.mp.dps - 8)

def F(y):
    y = mp.mpf(y)
    if y < 1:
        return y ** (-12) * mp.exp(2 * mp.pi * y - 2 * mp.pi / y) * F(1 / y)
    p = mp.mpf(1); n = 1
    while True:
        q = mp.exp(-2 * mp.pi * n * y)
        if q < TINY: break
        p *= (1 - q) ** 24; n += 1
    return p
def dlogF(y):
    if y < 1:
        return -12 / y + 2 * mp.pi + 2 * mp.pi / y ** 2 - dlogF(1 / y) / y ** 2
    return 48 * mp.pi * mp.fsum(sig[m] * mp.exp(-2 * mp.pi * m * y) for m in range(1, 80))
def dens(y):
    return F(y) * dlogF(y)
BRK = [0, 0.15, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24, 32]
def A_hecke(s, N=30):
    a = s + H
    return mp.fsum(tau[n] * (2 * mp.pi * n) ** (-a) * mp.gammainc(a, 2 * mp.pi * n) for n in range(1, N + 1))
def xi_hecke(s):
    return A_hecke(s) + A_hecke(1 - s)

# 1. the density and E_2
def E2(y):   # E_2(iy) = 1 - 24 sum sigma(m) e^{-2 pi m y}, y >= 0.4
    return 1 - 24 * mp.fsum(sig[m] * mp.exp(-2 * mp.pi * m * y) for m in range(1, NT + 1))
y = mp.mpf('0.5')
lhs = dens(1 / y)
rhs = 2 * mp.pi * y ** 12 * mp.exp(2 * mp.pi / y - 2 * mp.pi * y) * F(y) * (1 + y ** 2 * E2(y) - 6 * y / mp.pi)
put('density', {'f(2) = F(2) (log F)\'(2)': S(lhs), '2 pi F(2)(1 - E_2(2i))': S(2 * mp.pi * F(2) * (1 - E2(2))),
                'functional equation at y = 1/2, right side': S(rhs), 'integral of the density': S(mp.quad(dens, BRK), 30)})

# 2. mixed Poisson probabilities
P = {}
for k in range(1, 12):
    P[k] = mp.quad(lambda y: mp.gammainc(k, 2 * mp.pi * y, mp.inf, regularized=True) * dens(y), BRK)
pois = {}
for k in range(1, 12):
    lam_k = (2 * mp.pi) ** (-k) * mp.factorial(k - 1) * P[k]
    pois[k] = {'P(N(2 pi M*) <= k-1) = L(k - 11/2, Delta)': S(P[k]), '(2 pi)^-k (k-1)! P_k': S(lam_k),
               'Hecke xi(k - 11/2)': S(xi_hecke(mp.mpf(k) - H))}
put('mixed_poisson', pois)
put('functional_equation_identities', {'%d vs %d' % (k, 12 - k): S((2 * mp.pi) ** (-k) * mp.factorial(k - 1) * P[k] - (2 * mp.pi) ** (k - 12) * mp.factorial(11 - k) * P[12 - k], 5) for k in range(1, 6)})
put('E exp(-2 pi M*)', {'psi*(2 pi) by quadrature 2 pi int e^{-2 pi y} F(y) dy': S(2 * mp.pi * mp.quad(lambda y: mp.exp(-2 * mp.pi * y) * F(y), BRK + [48, 64])),
                         'P_1': S(P[1]), 'L(-9/2, Delta) from Hecke, 2 pi xi(-9/2)': S(2 * mp.pi * xi_hecke(-mp.mpf(9) / 2))})
# gamma mixture form at a = 11/2: E F(G_a/2 pi), G_a of density x^{a-1} e^{-x}/Gamma(a)
a = H
gm = mp.quad(lambda x: x ** (a - 1) * mp.exp(-x) / mp.gamma(a) * F(x / (2 * mp.pi)), [0, 1, 2, 4, 8, 16, 32, 64, 128])
put('gamma_mixture_s0', {'E F(G_{11/2}/2 pi)': S(gm), 'L(0, Delta) from Hecke': S(xi_hecke(0) * (2 * mp.pi) ** H / mp.gamma(H)),
                         'E Q(11/2, 2 pi M*)': S(mp.quad(lambda y: mp.gammainc(a, 2 * mp.pi * y, mp.inf, regularized=True) * dens(y), BRK))})
put('L(17/2, Delta)', {'E Q(14, 2 pi M*)': S(mp.quad(lambda y: mp.gammainc(14, 2 * mp.pi * y, mp.inf, regularized=True) * dens(y), BRK), 30),
                       'Hecke, xi(17/2) (2 pi)^14/13!': S(xi_hecke(mp.mpf(17) / 2) * (2 * mp.pi) ** 14 / mp.factorial(13), 30),
                       'sum tau(n) n^-14, n <= 400 (tail below 1e-19)': S(mp.fsum(mp.mpf(tau[n]) / mp.mpf(n) ** 14 for n in range(1, NT + 1)), 30)})
put('L(29/2, Delta)', {'E Q(20, 2 pi M*)': S(mp.quad(lambda y: mp.gammainc(20, 2 * mp.pi * y, mp.inf, regularized=True) * dens(y), BRK)),
                       'sum tau(n) n^-20, n <= 400': S(mp.fsum(mp.mpf(tau[n]) / mp.mpf(n) ** 20 for n in range(1, NT + 1)))})

# 3. Riemann's formula
def riemann(s, N=12):
    t = -1 / s - 1 / (1 - s)
    t += mp.fsum(mp.pi ** (-s / 2) * n ** (-s) * mp.gammainc(s / 2, mp.pi * n * n) + mp.pi ** (-(1 - s) / 2) * n ** (-(1 - s)) * mp.gammainc((1 - s) / 2, mp.pi * n * n) for n in range(1, N + 1))
    return t
rie = {}
for s in [mp.mpf(2), mp.mpc(0.5, 5)]:
    rie[S(s, 4)] = {'formula': S(riemann(s)), 'pi^{-s/2} Gamma(s/2) zeta(s)': S(mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s))}
put('riemann_formula', rie)

# 4. eta
def F1(y):
    y = mp.mpf(y)
    if y < 1:
        return y ** (-mp.mpf(1) / 2) * mp.exp(mp.pi * y / 12 - mp.pi / (12 * y)) * F1(1 / y)
    p = mp.mpf(1); n = 1
    while True:
        q = mp.exp(-2 * mp.pi * n * y)
        if q < TINY: break
        p *= 1 - q; n += 1
    return p
def dlogF1(y):
    if y < 1:
        return -1 / (2 * y) + mp.pi / 12 + mp.pi / (12 * y ** 2) - dlogF1(1 / y) / y ** 2
    return 2 * mp.pi * mp.fsum(sig[m] * mp.exp(-2 * mp.pi * m * y) for m in range(1, 80))
BRK1 = [0, 0.15, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 100]
chi12 = [0, 1, 0, 0, 0, -1, 0, -1, 0, 0, 0, 1]
s = mp.mpc('0.6', 5)
Lmax = mp.quad(lambda y: mp.gammainc(s / 2, mp.pi * y / 12, mp.inf, regularized=True) * F1(y) * dlogF1(y), BRK1)
Lref = mp.dirichlet(s, chi12)
def eta_ix(x):
    return mp.exp(-mp.pi * x / 12) * F1(x)
def A_eta(s):
    return mp.quad(lambda x: eta_ix(x) * x ** (s / 2 - 1), [1, 2, 4, 8, 16, 32, 64, 128, 256])
def A_eta_max(s):   # (pi/12)^{-s/2} E Gamma(s/2, (pi/12) max(M_1, 1))
    g = lambda y: mp.gammainc(s / 2, mp.pi * max(y, 1) / 12) * F1(y) * dlogF1(y)
    return (mp.pi / 12) ** (-s / 2) * mp.quad(g, BRK1)
xi12 = (12 / mp.pi) ** (s / 2) * mp.gamma(s / 2) * Lref
put('eta', {'L(0.6+5i, chi_12) = E Q(s/2, pi M_1/12)': S(Lmax), 'mpmath dirichlet': S(Lref),
            'xi(0.6+5i, chi_12)': S(xi12), 'A(s) + A(1-s), A by quadrature of eta': S(A_eta(s) + A_eta(1 - s)),
            'A(s) + A(1-s), A from the law of max(M_1, 1)': S(A_eta_max(s) + A_eta_max(1 - s)),
            'L(1/2, chi_12) = P(12 G_{1/4} > pi M_1) = E Q(1/4, pi M_1/12)': S(mp.quad(lambda y: mp.gammainc(mp.mpf(1) / 4, mp.pi * y / 12, mp.inf, regularized=True) * F1(y) * dlogF1(y), BRK1)),
            'L(1/2, chi_12), mpmath': S(mp.dirichlet(mp.mpf(1) / 2, chi12)),
            'E exp(-pi M_1/12) = L(2, chi_12), quadrature': S(mp.pi / 12 * mp.quad(lambda y: mp.exp(-mp.pi * y / 12) * F1(y), BRK1 + [200, 400, 800, 1600])),
            'pi^2/(6 sqrt 3)': S(mp.pi ** 2 / (6 * mp.sqrt(3))), 'L(2, chi_12), mpmath': S(mp.dirichlet(2, chi12))})

# 5. the conductor 11 curve, f = eta(tau)^2 eta(11 tau)^2 = q prod (1 - q^n)^2 (1 - q^{11n})^2 = sum a_n q^n
NA = 1200
poly = [0] * NA; poly[0] = 1
def mulfac(poly, step, power):
    for _ in range(power):
        for i in range(NA - 1, step - 1, -1):
            poly[i] -= poly[i - step]
for n in range(1, NA):
    mulfac(poly, n, 2)
    if 11 * n < NA:
        mulfac(poly, 11 * n, 2)
an = [0] + poly[:NA - 1]   # a_n = coefficient of q^n, n >= 1
assert an[1:12] == [1, -2, -1, 2, 1, 2, -2, 0, -2, -2, 1]
L1_series = 2 * mp.fsum(mp.mpf(an[n]) / n * mp.exp(-2 * mp.pi * n / mp.sqrt(11)) for n in range(1, NA))
def f11(y):   # f(iy) by the product for y >= 1/sqrt(11), by f(i/(11y)) = 11 y^2 f(iy) below
    y = mp.mpf(y)
    if y < 1 / mp.sqrt(11):
        return f11(1 / (11 * y)) / (11 * y * y)
    p = mp.exp(-2 * mp.pi * y); n = 1
    while True:
        q = mp.exp(-2 * mp.pi * n * y)
        if q < TINY: break
        p *= (1 - q) ** 2 * (1 - q ** 11) ** 2; n += 1
    return p
def f11_product(y):
    y = mp.mpf(y); p = mp.exp(-2 * mp.pi * y)
    for n in range(1, 4000):
        q = mp.exp(-2 * mp.pi * n * y)
        if q < TINY: break
        p *= (1 - q) ** 2 * (1 - q ** 11) ** 2
    return p
F11 = lambda y: mp.exp(2 * mp.pi * y) * f11(y)
L1_max = 2 * mp.pi * mp.quad(lambda y: mp.exp(-2 * mp.pi * y) * F11(y), [0, 0.05, 0.1, 0.2, 1 / mp.sqrt(11), 0.5, 1, 2, 4, 8, 16, 32])
put('E11', {'a_2..a_13': an[2:14], 'Fricke check at y = 0.2: product / (f(i/2.2)/(11 * 0.04)) - 1': S(f11_product('0.2') / (f11_product(1 / (11 * mp.mpf('0.2'))) / (11 * mp.mpf('0.04'))) - 1, 5),
            'L(1, E) = E exp(-2 pi M_11)': S(L1_max), 'L(1, E), series 2 sum a_n e^{-2 pi n/sqrt 11}/n': S(L1_series)})

# 6. theta_4, M_4 and the Kolmogorov variable
def F4_product(y):
    y = mp.mpf(y); p = mp.mpf(1); m = 1
    while True:
        q1 = mp.exp(-2 * mp.pi * m * y); q2 = mp.exp(-mp.pi * (2 * m - 1) * y)
        if q2 < TINY: break
        p *= (1 - q1) * (1 - q2) ** 2; m += 1
    return p
def F4(y):   # theta_4(y): product for y >= 0.25, theta_4(y) = y^{-1/2} theta_2(1/y) below
    y = mp.mpf(y)
    if y >= mp.mpf('0.25'):
        return F4_product(y)
    return y ** (-mp.mpf(1) / 2) * 2 * mp.fsum(mp.exp(-mp.pi * (m + mp.mpf(1) / 2) ** 2 / y) for m in range(0, 30))
th4series = lambda y: 1 + 2 * mp.fsum((-1) ** n * mp.exp(-mp.pi * n * n * mp.mpf(y)) for n in range(1, 60))
B4 = [0, 0.05, 0.1, 0.25, 0.5, 1, 2, 4, 8, 16, 32]
EK = mp.sqrt(mp.pi / 2) * mp.quad(lambda y: (1 - F4(y)) / (2 * mp.sqrt(y)), B4 + [64, 128])
mom = {}
for s in [mp.mpf(3), mp.mpc(0.5, 5)]:
    eps = mp.mpf('0.05')   # E M_4^{s/2} = eps^{s/2} - (s/2) int_0^eps y^{s/2-1} F4 + (s/2) int_eps^oo y^{s/2-1} (1 - F4)
    EM = eps ** (s / 2) - (s / 2) * mp.quad(lambda y: y ** (s / 2 - 1) * F4(y), [0, mp.mpf('0.01'), mp.mpf('0.02'), eps]) \
        + (s / 2) * mp.quad(lambda y: y ** (s / 2 - 1) * (1 - F4(y)), [eps] + B4[3:] + [64, 128])
    mom[S(s, 4)] = {'E[(2K^2)^{s/2}]/(2 Gamma(1+s/2))': S(mp.pi ** (s / 2) * EM / (2 * mp.gamma(1 + s / 2))), '(1 - 2^{1-s}) zeta(s)': S((1 - 2 ** (1 - s)) * mp.zeta(s))}
put('theta4', {'product / theta series - 1 at y = 0.3, 0.7': [S(F4_product('0.3') / th4series('0.3') - 1, 5), S(F4_product('0.7') / th4series('0.7') - 1, 5)],
               'E exp(-3 M_4), quadrature': S(3 * mp.quad(lambda y: mp.exp(-3 * y) * F4(y), B4)),
               'sqrt(3 pi)/sinh sqrt(3 pi)': S(mp.sqrt(3 * mp.pi) / mp.sinh(mp.sqrt(3 * mp.pi))),
               'E K, quadrature of 1 - theta_4': S(EK), 'sqrt(pi/2) log 2': S(mp.sqrt(mp.pi / 2) * mp.log(2)), 'moments': mom})

# 7. M_1
b = lambda n: mp.pi * (n * n - 1) / 12
pk = {}
for k in [1, 2, 3, 10]:
    lo, c, hi = b(6 * k - 1), b(6 * k), b(6 * k + 1)
    pk[k] = S((1 - lo / c) * hi / (hi - lo), 15)
def lt1(th):
    w = mp.sqrt(1 - 12 * mp.mpc(th) / mp.pi)
    return mp.pi * (1 - w ** 2) * mp.sin(mp.pi * w / 3) / (2 * mp.sqrt(3) * w * mp.cos(mp.pi * w / 2))
th = mp.mpf(2)
d_closed = -mp.diff(lambda x: mp.log(mp.re(lt1(x))), th)
d_sum = mp.nsum(lambda k: 1 / (th + b(6 * k - 1)) + 1 / (th + b(6 * k + 1)) - 1 / (th + b(6 * k)), [1, mp.inf])
# int_0^oo N_1(z) (theta + z)^{-2} dz, exactly over the first 60000 atoms, then N_1 by its asymptote sqrt(3z/pi)/3
KB = 20000
atoms = sorted([(b(n), 1) for n in range(5, 6 * KB + 2) if n % 6 in (1, 5)] + [(b(6 * k), -1) for k in range(1, KB + 1)])
Nv = 0; acc = mp.mpf(0)
for j in range(len(atoms) - 1):
    Nv += atoms[j][1]
    acc += Nv * (1 / (th + atoms[j][0]) - 1 / (th + atoms[j + 1][0]))
Nv += atoms[-1][1]; Z0 = atoms[-1][0]
acc += mp.quad(lambda z: (mp.sqrt(12 * z / mp.pi + 1) / 6 - mp.mpf(1) / 2) / (th + z) ** 2, [Z0, mp.inf])   # N_1 averages X/6 - 1/2, X = sqrt(12z/pi + 1)
def xnu1_series(x):   # x nu_1(x) = k_1(x) = sum_k e^{-b_{6k-1} x} + e^{-b_{6k+1} x} - e^{-b_{6k} x}
    return mp.fsum(mp.exp(-b(6 * k - 1) * x) + mp.exp(-b(6 * k + 1) * x) - mp.exp(-b(6 * k) * x) for k in range(1, 400))
def xnu1_theta(x):    # the theta form of k_1, from Poisson summation
    tht = lambda al: mp.fsum(mp.exp(-12 * mp.pi * (m + al) ** 2 / x) for m in range(-40, 41))
    return mp.exp(mp.pi * x / 12) / 2 * (1 + mp.sqrt(12 / x) * (tht(0) / 6 - tht(mp.mpf(1) / 2) / 2 - 2 * tht(mp.mpf(1) / 3) / 3)) - 1
def xnu1(x):
    x = mp.mpf(x)
    return xnu1_series(x) if x > 1 else xnu1_theta(x)
lev = mp.quad(lambda x: (1 - mp.exp(-2 * x)) * xnu1(x) / x, [0, mp.mpf('1e-6'), mp.mpf('1e-3'), 0.05, 0.25, 1, 2, 4, 8, 16, 32])
small = {S(x, 3): S(2 * x * (xnu1(x) / x - x ** (-mp.mpf(3) / 2) / (2 * mp.sqrt(3))), 12) for x in [mp.mpf('1e-2'), mp.mpf('1e-4'), mp.mpf('1e-6')]}
large = {S(x, 3): S(xnu1(x) * mp.exp(2 * mp.pi * x), 12) for x in [mp.mpf(2), mp.mpf(5), mp.mpf(10)]}
both = {S(x, 3): [S(xnu1_theta(x), 20), S(xnu1_series(x), 20)] for x in [mp.mpf('0.05'), mp.mpf('0.3'), mp.mpf(1), mp.mpf(2)]}
put('M1', {'p_k': pk, '-psi_1\'/psi_1 at 2, closed form': S(d_closed), 'sum over the Thorin atoms': S(d_sum),
           'int N_1(z)(2+z)^-2 dz': S(acc, 20), 'Laplace exponent at 2, closed form': S(-mp.log(mp.re(lt1(2)))),
           'Laplace exponent at 2, int (1 - e^{-2x}) nu_1(x) dx': S(lev),
           '2x(nu_1(x) - x^{-3/2}/(2 sqrt 3)) at small x (tends to -1)': small, 'x nu_1(x) e^{2 pi x} at large x (tends to 1)': large,
           'k_1(x): theta form and series': both})

# 8. M*: the Laplace transform at large theta, and the cumulants of the tilted laws
def psi_quad(th):
    th = mp.mpf(th); ys = mp.sqrt(2 * mp.pi / th)
    pts = sorted(set([mp.mpf(0)] + [ys * f for f in [0.25, 0.5, 0.75, 1, 1.5, 2, 3, 5, 8]] + [mp.mpf(x) for x in [1, 2, 4, 8, 16, 32]]))
    return th * mp.quad(lambda y: mp.exp(-th * y) * F(y), pts)
asy = {}
for th in [10, 100, 1000, 10000, 100000]:
    v = psi_quad(th)
    approx = mp.mpf(th) ** (mp.mpf(25) / 4) * (2 * mp.pi) ** (-mp.mpf(21) / 4) * 2 ** (-mp.mpf(1) / 2) * mp.exp(-2 * mp.sqrt(2 * mp.pi * th))
    asy[str(th)] = {'-log psi*': S(-mp.log(v), 15), '2 sqrt(2 pi theta)': S(2 * mp.sqrt(2 * mp.pi * th), 15), 'psi*/asymptote': S(v / approx, 10)}
put('psi_star_large_theta', asy)
mp.mp.dps = 50
cum = {}
for th in [mp.mpf('0.5'), mp.mpf(2), mp.mpf(8), mp.mpf(30)]:
    ys = mp.sqrt(2 * mp.pi / th)
    pts = sorted(set(BRK + [ys * f for f in [0.5, 1, 2]]))
    mu = [mp.quad(lambda y: y ** j * mp.exp(-th * y) * dens(y), pts) for j in range(9)]
    c = [(-1) ** j * mu[j] / mp.factorial(j) for j in range(9)]   # psi*(theta + h) = sum c_j h^j
    l = [mp.log(c[0])] + [mp.mpf(0)] * 8
    for j in range(1, 9):
        l[j] = (j * c[j] - mp.fsum(i * l[i] * c[j - i] for i in range(1, j))) / (j * c[0])
    kap = [(-1) ** j * mp.factorial(j) * l[j] for j in range(1, 9)]   # cumulants of the tilted law
    cum[S(th, 3)] = [S(k, 8) for k in kap]
put('tilted_cumulants_M_star', cum)
mp.mp.dps = 40

# 9. finite maxima of independent exponential variables
def grouped(rates):
    # P(max <= y) = prod_i (1 - e^{-lambda_i y}) = 1 + sum_p c_p e^{-p y}, the exponents p grouped exactly (rates are Fractions)
    coef = {Fraction(0): 1}
    for lam in rates:
        new = dict(coef)
        for e, c in coef.items():
            new[e + lam] = new.get(e + lam, 0) - c
        coef = new
    return {p: c for p, c in coef.items() if p != 0 and c != 0}
def max_transform(rates):
    # psi(theta) = 1 + sum_p c_p theta/(theta + p) = P(theta)/Q(theta), Q = prod_p (theta + p); coefficients in increasing powers
    coef = grouped(rates)
    poles = sorted(coef)
    Q = [Fraction(1)]
    for p in poles:
        Q = [Fraction(0)] + Q
        for i in range(len(Q) - 1):
            Q[i] += p * Q[i + 1]
    Pn = list(Q)
    for p in poles:
        # Q/(theta + p) by synthetic division from the top, then times c_p theta
        d = len(Q) - 1
        qq = [Fraction(0)] * d
        qq[d - 1] = Q[d]
        for i in range(d - 1, 0, -1):
            qq[i - 1] = Q[i] - p * qq[i]
        assert Q[0] == p * qq[0]
        for i, x in enumerate(qq):
            Pn[i + 1] += coef[p] * x
    while Pn and Pn[-1] == 0:
        Pn.pop()
    return poles, coef, Pn
def analyse(rates):
    poles, coef, Pn = max_transform(rates)
    deg = len(Pn) - 1
    cf = [mp.mpf(x.numerator) / x.denominator for x in reversed(Pn)]
    roots = []
    if deg > 0:
        roots = None
        for steps, extra in [(400, 400), (3000, 1500), (12000, 4000)]:
            try:
                roots = mp.polyroots(cf, maxsteps=steps, extraprec=extra); break
            except mp.NoConvergence:
                pass
        if roots is None:
            return {'n': len(rates), 'poles': len(poles), 'zeros': deg, 'failed': True}
    resid = max([abs(mp.polyval(cf, r)) / mp.polyval([abs(c) for c in cf], abs(r)) for r in roots] + [mp.mpf(0)])
    zs = [-r for r in roots]                      # psi vanishes at theta = -z
    pm = poles[0]; pmf = mp.mpf(pm.numerator) / pm.denominator
    pf = [mp.mpf(p.numerator) / p.denominator for p in poles]
    def k(x):   # the canonical function: -log psi(theta) = int_0^oo (1 - e^{-theta x}) k(x) dx/x; poles are simple
        return mp.fsum(mp.exp(-p * x) for p in pf) - mp.re(mp.fsum(mp.exp(-z * x) for z in zs))
    xs = [mp.mpf(10) ** (-3 + mp.mpf(4) * j / 60) for j in range(61)] + [mp.mpf(j) / 4 for j in range(41, 161)]
    vals = [k(x / pmf) * mp.exp(x) for x in xs]
    nonreal = sum(1 for z in zs if abs(mp.im(z)) > mp.mpf(10) ** -20)
    minre = min([mp.re(z) for z in zs]) / pmf if zs else mp.inf
    return {'n': len(rates), 'poles': len(poles), 'zeros': deg, 'nonreal zeros': nonreal, 'min Re z / min pole': float(minre),
            'min of k(x) e^{x min pole} on the grid': float(min(vals)), 'k(0+) = poles - zeros': len(poles) - deg,
            'largest relative residual of the roots': float(resid)}
rng = random.Random(20261009)
fam = []
for n, cnt in [(2, 300), (3, 300), (4, 300), (5, 300), (6, 40)]:
    for _ in range(cnt):
        rates = [Fraction(rng.randint(100, 5000), 1000) for _ in range(n)]
        fam.append(analyse(rates))
    print('families with', n, 'rates done', round(time.time() - t0, 1), flush=True)
summ = {}
for n in [2, 3, 4, 5, 6]:
    fs = [f for f in fam if f['n'] == n and not f.get('failed')]
    summ[n] = {'families': len([f for f in fam if f['n'] == n]), 'root finding failed': len([f for f in fam if f['n'] == n and f.get('failed')]),
               'with nonreal zeros': sum(1 for f in fs if f['nonreal zeros'] > 0), 'max nonreal zeros': max(f['nonreal zeros'] for f in fs),
               'min over families of min Re z / min pole': min(f['min Re z / min pole'] for f in fs),
               'min over families of min k(x) e^{x min pole}': min(f['min of k(x) e^{x min pole} on the grid'] for f in fs),
               'families with k < 0 somewhere on the grid': sum(1 for f in fs if f['min of k(x) e^{x min pole} on the grid'] < 0),
               'largest relative residual': max(f['largest relative residual of the roots'] for f in fs)}
put('finite_maxima_summary', summ)
special = {}
for name, rates in [('1,2', [1, 2]), ('1,2,3', [1, 2, 3]), ('1,2,3,4', [1, 2, 3, 4]), ('1,2,3,4,5', [1, 2, 3, 4, 5]), ('1,1,1', [1, 1, 1]),
                    ('1,10,100', [1, 10, 100]), ('1,1.01,1.02,5', [1, Fraction(101, 100), Fraction(102, 100), 5])]:
    special[name] = analyse([Fraction(r) for r in rates])
put('finite_maxima_special', special)
# truncations of the eta and Delta laws: rates 1, ..., N once (eta, in units of 2 pi), and 24 times each (Delta)
trunc = {}
for N in range(4, 11):
    trunc['eta, N = %d' % N] = analyse([Fraction(r) for r in range(1, N + 1)])
    print('eta truncation', N, round(time.time() - t0, 1), flush=True)
for N in (2, 3):
    trunc['Delta, N = %d' % N] = analyse([Fraction(r) for r in range(1, N + 1) for _ in range(24)])
    print('Delta truncation', N, round(time.time() - t0, 1), flush=True)
put('truncations', trunc)
out['seconds'] = time.time() - t0
json.dump(out, open('maxexp_checks.json', 'w'), indent=1)
print('done', round(time.time() - t0, 1))
