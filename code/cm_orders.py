# cm_orders.py: complete monotonicity of the heat trace at every order on an interval, for Section ch:heat (ch/ch12.tex).
# Added for the book in September 2026, from the checks of the authors' note on complete monotonicity at all orders
# (misc/reviews/nick-cm/), and sharpened at large orders in October 2026 (misc/reviews/nick-pathset/). With
# M_k(t) = Gamma(k+1/2)/(2 pi t^(k+1/2)) the explicit formula gives, for every k >= 1,
#   2 c_k(t)/M_k(t) >= E_k(t) - S(t) - P_k(t),
# where E_k(t) is the average of m(r) = Re psi(1/4 + i r/2) - log pi under the weight r^(2k) e^(-r^2 t) on [0, oo),
# S(t) = 2 sum_n Lambda(n) n^(-1/2) A(x_n) e^(-x_n^2/2) with x_n = log n/(2 sqrt t) bounds the prime terms through
# |1F1(k+1/2; 1/2; -x^2)| <= A(x) e^(-x^2/2), A(x) = max(2.62, (pi (2e x^2 + 1/2))^(1/4)), uniformly in k (the second
# branch from DLMF (18.14.9), e^(-x^2/2)|H_n(x)| <= (2^n n!)^(1/2); the note used Cramer's constant 1.0865 there), and
# P_k(t) = 4 pi 4^(-k) e^(t/4) t^(k+1/2)/Gamma(k+1/2) is the polar term. For orders k > K = 10^17, part (b) of the lemma
# replaces A(x_n) by a constant C < 1 + 6e-7 for the prime powers n <= 2^24 once t >= T0; S_C(t) denotes S(t) so changed.
# E_k increases in k and decreases in t, S, S_C and P_k increase in t, and P_k decreases in k for t <= 2, so one value of
# t settles a whole interval:
#  (a) T0 = 0.019 with k = 1: every order k >= 1 on (0, T0], with no information on the zeros (the note's theorem, which
#      has 0.018 with Cramer's constant);
#  (b) T1 with k = K + 1 and S_C: every order k > K on [T0, T1], which with (a) and Theorem thm:ch12:cm (orders k <= K at
#      every t > 0) makes W completely monotone on (0, T1].
# Every inequality of (a) and (b) is evaluated in ball arithmetic (python-flint, Arb) and counts only if it holds for the
# whole ball. Outside the certificate the script computes the margins in floating point (mpmath, numpy), and recomputes
# the numbers of the note, with its constant 1.0865, for comparison.
# Run from a folder containing zeros_1700.txt: python3 cm_orders.py [T1]; writes cm_orders.json there.
import json
import math
import sys
import time

import mpmath as mp
import numpy as np
from flint import arb, acb, ctx

ctx.prec = 128
pi, half, E = arb.pi(), arb(1) / 2, arb(1).exp()
K = 10**17
T0 = arb(19) / 1000
out, holds = {'certificate': {}, 'float': {}}, []
T00 = time.time()


def s(x, n=12):
    return x.str(n, radius=False) if isinstance(x, arb) else x


def check(name, cond, **vals):
    holds.append(bool(cond))
    out['certificate'][name] = dict(holds=bool(cond), **{k: s(v) for k, v in vals.items()})
    print(name, bool(cond), round(time.time() - T00, 1), flush=True)


def m(r):  # m(r) = Re psi(1/4 + i r/2) - log pi, for a real ball r
    return acb(arb(1) / 4, r / 2).digamma().real - pi.log()


C1, C2 = arb('2.62'), arb(1)
C2T = arb('2.08')  # the tail constant printed in the book, an upper bound for c2 = (pi (2e + 1/2))^(1/4) = 2.078


def A_up(x):  # an upper bound for A(x) on the ball x
    u, v = C1.upper(), (C2 * (pi * (2 * E * x * x + half)).sqrt().sqrt()).upper()
    return u if u > v else v


def prime_powers(N):  # (n, log p) for the prime powers n = p^j <= N
    sieve = np.ones(N + 1, bool)
    sieve[:2] = False
    for p in range(2, int(N**0.5) + 1):
        if sieve[p]:
            sieve[p * p::p] = False
    res = []
    for p in np.nonzero(sieve)[0].tolist():
        q = p
        while q <= N:
            res.append((q, p))
            q *= p
    return sorted(res)


def S_up(t, N):  # an upper bound for S(t): the prime powers n <= N in Arb, and the tail n > N bounded as below
    lt, total, logs = 2 * t.sqrt(), arb(0), {}
    for n, p in prime_powers(N):
        if p not in logs:
            logs[p] = arb(p).log()
        x = arb(n).log() / lt
        total += (2 * logs[p] / arb(n).sqrt() * A_up(x) * (-x * x / 2).exp()).upper()
    # For n > N >= e^10 and x_n >= 1: e^(-x_n^2/2) <= n^(-q) with q = log N/(8t), A(x) <= 2.62 + 2.08 sqrt x (check
    # tail_constant), and log n, (log n)^(3/2) <= n^(1/2), so the tail is at most
    # 2 (2.62 + 2.08 (2 sqrt t)^(-1/2)) sum_{n > N} n^(-q) <= 2 (2.62 + 2.08 (2 sqrt t)^(-1/2)) N^(1-q)/(q - 1).
    q = arb(N).log() / (8 * t)
    assert N >= math.exp(10) and q > 1 and arb(N).log() / lt >= 1
    tail = (2 * (C1 + C2T / lt.sqrt()) * arb(N) ** (1 - q) / (q - 1)).upper()
    return (total + tail).upper(), tail


def P_up(k, t):  # the polar bound P_k(t), upper endpoint
    return (4 * pi * arb(4) ** (-k) * (t / 4).exp() * t ** (k + half) / (arb(k) + half).gamma()).upper()


# (a) E_1(T0) from below by lower sums, since m is nondecreasing on [0, oo): over [0, 64] on dyadic pieces of width 1/1024
# up to 32 and 1/128 beyond, with the exact weights w = int r^2 e^(-r^2 t) dr from the antiderivative
# G(r) = sqrt(pi) erf(r sqrt t)/(4 t^(3/2)) - r e^(-r^2 t)/(2t); beyond 64 the integrand is positive and is dropped.
def E1_low(t):
    D, st = pi.sqrt() / (4 * t ** (arb(3) / 2)), t.sqrt()
    G = lambda r: pi.sqrt() * (r * st).erf() / (4 * t ** (arb(3) / 2)) - r * (-r * r * t).exp() / (2 * t)
    grid = [arb(i) / 1024 for i in range(32 * 1024 + 1)] + [32 + arb(j) / 128 for j in range(1, 32 * 128 + 1)]
    Gs = [G(r) for r in grid]
    N = arb(0)
    for i in range(len(grid) - 1):
        N += m(grid[i]).lower() * (Gs[i + 1] - Gs[i])
    return (N / D).lower(), m(arb(64))


e1, m64 = E1_low(T0)
check('a_m64_positive', m64 > 0, m64=m64)
Sa, taila = S_up(T0, 30000)
Pa = P_up(1, T0)
check('a_margin', e1 - Sa - Pa > 0, E1_lower=e1, S_upper=Sa, S_tail=taila, P_upper=Pa, margin=e1 - Sa - Pa)
check('a_tail_below_1e-299', taila < arb(10) ** -299, S_tail=taila)
c2 = C2 * (pi * (2 * E + half)).sqrt().sqrt()
check('tail_constant', c2 < C2T, c2=c2)
# Kershaw: Gamma(y+1)/Gamma(y+1/2) < (y - 1/2 + (3/4)^(1/2))^(1/2), and (3/4)^(1/2) - 1/2 = (sqrt 3 - 1)/2 < 0.37.
check('kershaw_constant', (arb(3).sqrt() - 1) / 2 < arb('0.37'), kershaw=(arb(3).sqrt() - 1) / 2)
# At y = 0, outside Kershaw's range y > 0, the ratio is 1/sqrt(pi), which lies in [sqrt(1/4), sqrt(0.37)].
check('kershaw_at_0', (half < 1 / pi.sqrt()) and (1 / pi.sqrt() < arb('0.37').sqrt()), ratio=1 / pi.sqrt())
# A(x) e^(-x^2/2) decreases on (0, oo): its first branch does, its second does for x >= 1/sqrt 2 (logarithmic derivative
# below 1/(2x) - x), and below 1/sqrt 2 the second branch is under 2.62, being increasing with this value at 1/sqrt 2.
b2 = C2 * (pi * (E + half)).sqrt().sqrt()
check('A_branch_at_1_over_sqrt2', b2 < C1, second_branch=b2)
check('lemma_constants', (arb(2).sqrt() + arb('1.2') < C1) and (arb('1.37').sqrt() < arb('1.2')),
      sqrt2_plus_1_2=arb(2).sqrt() + arb('1.2'), k1_tail=arb('1.37').sqrt())

# (b) For k > K the binomial sum in the proof of Lemma lem:ch12:hermite is split at j = eps k instead of k/2 (part (b) of
# the lemma): |F_k(x)| <= C_eps(k) e^(-x^2/2) whenever e^2 x^2 <= 4 eps k, where
#   C_eps(k) = sqrt((k + 0.37)/((1 - eps) k + 1/4)) + 2 sqrt(k + 0.37) e^(-eps k).
# The first term decreases in k because 1/4 < 0.37 (1 - eps), and the second for k > 1/(2 eps), so for k >= K + 1 the
# constant is at most C = C_eps(K + 1), and the condition holds for every x <= X = 2 sqrt(eps (K + 1))/e. For t >= T0
# every prime power n <= 2^24 has x_n <= log(2^24)/(2 sqrt T0) < X, and the prime powers above 2^24 keep A(x_n).
EPSJ, NB = arb(1) / 10**6, 2**24
K1 = arb(K + 1)
C_large = (((K1 + arb('0.37')) / ((1 - EPSJ) * K1 + arb(1) / 4)).sqrt()
           + 2 * (K1 + arb('0.37')).sqrt() * (-EPSJ * K1).exp()).upper()
X_large = 2 * (EPSJ * K1).sqrt() / E
x_cut = arb(NB).log() / (2 * T0.sqrt())
check('b_lemma_monotone', (arb(1) / 4 < arb('0.37') * (1 - EPSJ)) and (K1 > 1 / (2 * EPSJ)))
check('b_lemma_constant', C_large < 1 + 6 * arb(10) ** -7, C=C_large)
check('b_lemma_range', x_cut < X_large, x_at_cutoff=x_cut, X=X_large)


def SC_up(t, N):  # an upper bound for S_C(t): the constant C for the prime powers n <= N, the tail n > N as in S_up
    lt, total, logs = 2 * t.sqrt(), arb(0), {}
    for n, p in prime_powers(N):
        if p not in logs:
            logs[p] = arb(p).log()
        ln = arb(n).log()
        total += 2 * logs[p] / arb(n).sqrt() * (-(ln * ln) / (8 * t)).exp()
    q = arb(N).log() / (8 * t)
    assert N >= math.exp(10) and q > 1 and arb(N).log() / lt >= 1
    tail = (2 * (C1 + C2T / lt.sqrt()) * arb(N) ** (1 - q) / (q - 1)).upper()
    return (C_large * total + tail).upper(), tail


# E_{K+1}(T1) from below: the weight makes v = r^2 t a Gamma(alpha) variable, alpha = K + 3/2, and with r1 = sqrt(a/T1),
# a = alpha (1 - eps)^2, P(r < r1) = P(v < a) <= (a/alpha)^alpha e^(alpha - a) (Chernoff). Since m is nondecreasing,
# E_{K+1}(T1) >= m(r1) - (m(r1) - m(0)) P(r < r1).
mp.mp.dps = 30


def m_f(r):
    return float(mp.re(mp.digamma(mp.mpf(1) / 4 + 0.5j * r)) - mp.log(mp.pi))


PP = prime_powers(10**6)
PN = np.array([n for n, p in PP], float)
PL = np.array([math.log(p) for n, p in PP])


def A_f(x, c=1.0):  # c = 1 for the book (DLMF (18.14.9)), c = 1.0865 for the note (Cramer)
    return np.maximum(2.62, c * (math.pi * (2 * math.e * x * x + 0.5)) ** 0.25)


def S_f(t, env=False, c=1.0):
    x = np.log(PN) / (2 * math.sqrt(t))
    a = 0.9 * x ** (1 / 3) if env else A_f(x, c)
    return float(np.sum(2 * PL / np.sqrt(PN) * a * np.exp(-x * x / 2)))


def S_plain_f(t):  # S(t) with A replaced by one, twice the sum S(t) of Section ch:lessons
    x = np.log(PN) / (2 * math.sqrt(t))
    return float(np.sum(2 * PL / np.sqrt(PN) * np.exp(-x * x / 2)))


def P_f(k, t):
    return float(4 * mp.pi * mp.mpf(4) ** (-k) * mp.exp(t / 4) * mp.mpf(t) ** (k + 0.5) / mp.gamma(k + 0.5))


def E_f(k, t):
    w = lambda r: r ** (2 * k) * mp.exp(-r * r * t)
    num = mp.quad(lambda r: w(r) * (mp.re(mp.digamma(mp.mpf(1) / 4 + 0.5j * r)) - mp.log(mp.pi)),
                  [0, 2, 5, 10, 20, 30, 40, 60, 100, 200, mp.inf])
    return float(num / (mp.gamma(k + 0.5) / (2 * mp.mpf(t) ** (k + 0.5))))


def EK_f(t):  # E_{K+1}(t) = (psi(K + 3/2) - log t)/2 - log 2 pi, up to a relative correction of order t/K
    return float((mp.digamma(K + 1.5) - mp.log(t)) / 2 - mp.log(2 * mp.pi))


def bisect(f, lo, hi, it=60):
    flo = f(lo)
    for _ in range(it):
        mid = (lo + hi) / 2
        fm = f(mid)
        if (fm > 0) == (flo > 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return (lo + hi) / 2


# The phase free limit: with C = 1 the bound E_{K+1} - S_C is E_{K+1} - S_plain, and its root bounds every T1 that
# ignores the phases of the primes; with A(x_n) for every n (the lemma's part (a) alone) the root is far smaller.
t_root = bisect(lambda t: EK_f(t) - S_plain_f(t), 0.5, 2.0, 50)  # P_{K+1} is below e^(-1000) there
out['float']['combined_root'] = t_root
out['float']['combined_root_A'] = bisect(lambda t: EK_f(t) - S_f(t), 0.05, 1.0, 50)
alpha, eps = arb(K) + arb(3) / 2, arb(1) / 8192
a = alpha * (1 - eps) ** 2
logp = alpha * (2 * (1 - eps).log() + 1 - (1 - eps) ** 2)
check('b_chernoff', logp < -1000, log_bound=logp)
m0 = m(arb(0))
Pb = arb(-1000).exp().upper()


def b_bounds(T1):  # a lower bound for E_{K+1}(T1) and an upper bound for S_C(T1)
    mr1 = m((a / T1).sqrt())
    return ((mr1 - (mr1 - m0) * arb(-1000).exp()).lower(),) + SC_up(T1, NB)


# T1 is the value given, or else the largest multiple of 1/100 below the float root at which the margin holds.
if len(sys.argv) > 1:
    T1 = arb(sys.argv[1])
    eK, Sb, tailb = b_bounds(T1)
else:
    j0 = math.floor(t_root * 100)
    for j in range(j0, max(j0 - 10, 0), -1):
        T1 = arb(j) / 100
        eK, Sb, tailb = b_bounds(T1)
        if eK - Sb - Pb > 0:
            break
check('b_T1_between_T0_and_2', (T0 < T1) and (T1 <= 2), T1=T1)
logP = (4 * pi).log() - (K + 1) * arb(4).log() + T1 / 4 + (K + arb(3) / 2) * T1.log() - alpha.lgamma()
check('b_polar_small', logP < -1000, log_P=logP)
check('b_tail_below_1e-3', tailb < arb(10) ** -3, S_tail=tailb)
check('b_margin', eK - Sb - Pb > 0, T1=T1, E_K1_lower=eK, S_C_upper=Sb, S_tail=tailb, P_upper=Pb, margin=eK - Sb - Pb,
      E_K1_float=EK_f(float(T1.mid())))
out['T1'] = s(T1)
out['all_hold'] = all(holds)
print('certificate', all(holds), 'T1', s(T1), 'float root', t_root, flush=True)

# Floating point margins, and the numbers of the note with its constant 1.0865 (outside the certificate).
F = out['float']
note = {0.005: (0.828955, 2.92e-5, 0.00126, 0.8277), 0.010: (0.481685, 0.00996, 0.00355, 0.4682),
        0.015: (0.278178, 0.0669, 0.00654, 0.2047), 0.018: (0.186522, 0.1255, 0.00860, 0.0524),
        0.019: (0.159319, 0.1481, 0.00933, 0.0019), 0.020: (0.133501, 0.1719, 0.01008, -0.0484)}
rows = []
for t, vals in note.items():
    e, sv, sn, p = E_f(1, t), S_f(t), S_f(t, c=1.0865), P_f(1, t)
    rows.append(dict(t=t, E1=e, S=sv, P=p, margin=e - sv - p, S_note=sn, margin_note=e - sn - p, note=vals))
F['table'] = rows
F['root_margin'] = bisect(lambda t: E_f(1, t) - S_f(t) - P_f(1, t), 0.019, 0.021, 40)
F['root_margin_note'] = bisect(lambda t: E_f(1, t) - S_f(t, c=1.0865) - P_f(1, t), 0.018, 0.020, 40)
F['combined_root_note'] = bisect(lambda t: EK_f(t) - S_f(t, c=1.0865), 0.05, 1.0, 50)
F['combined_margin_A'] = {t: EK_f(t) - S_f(t) for t in (0.1, 0.2, 0.3, 0.4, 0.5)}
F['combined_margin'] = {t: EK_f(t) - float(C_large.mid()) * S_plain_f(t) for t in (0.1, 0.5, 1.0, 1.2, 1.25, 1.29)}
F['S_plain'] = {t: S_plain_f(t) for t in (0.31, 0.45, 0.51, 0.57, 1.0, 1.18, 1.29)}
F['root_E1'] = bisect(lambda t: E_f(1, t), 0.02, 0.03, 40)
F['root_margin_envelope'] = bisect(lambda t: E_f(1, t) - S_f(t, True) - P_f(1, t), 0.018, 0.024, 40)
print('table and roots', round(time.time() - T00, 1), flush=True)

# The explicit formula at t = 0.018 for k = 1, 3, 10: 2 c_k/M_k from the first 1700 zeros against
# E_k - 2 sum_{n <= 10^4} Lambda(n) n^(-1/2) F_k(x_n) + 2 (-1/4)^k e^(t/4)/M_k, with F_k(x) = e^(-x^2) 1F1(-k; 1/2; x^2).
mp.mp.dps = 50
g = [mp.mpf(l.split()[-1]) for l in open('zeros_1700.txt') if l.strip() and not l.startswith('#')]
t = mp.mpf(18) / 1000
ef = []
for k in (1, 3, 10):
    Mk = mp.gamma(k + 0.5) / (2 * mp.pi * t ** (k + 0.5))
    lhs = 2 * sum(y ** (2 * k) * mp.exp(-y * y * t) for y in g) / Mk
    prime = sum(2 * mp.log(p) / mp.sqrt(n) * mp.exp(-(mp.log(n) / (2 * mp.sqrt(t))) ** 2)
                * mp.hyp1f1(-k, 0.5, (mp.log(n) / (2 * mp.sqrt(t))) ** 2) for n, p in PP if n <= 10**4)
    rhs = mp.mpf(E_f(k, 0.018)) - prime + 2 * (-mp.mpf(1) / 4) ** k * mp.exp(t / 4) / Mk
    ef.append(dict(k=k, lhs=mp.nstr(lhs, 15), rhs_float_E=mp.nstr(rhs, 12)))
F['explicit_formula'] = ef
F['explicit_formula_note'] = ['0.193230228622', '0.746236943326', '1.35591436993']
print('explicit formula', round(time.time() - T00, 1), flush=True)


# |F_k(x)| e^(x^2/2) = pi^(1/4) (4^k (k!)^2/(2k)!)^(1/2) |psi_{2k}(x)|, with psi_n the normalised Hermite functions
# computed by their three term recurrence; checked against mpmath at a few points.
def hermite_sup(xs, kmax):
    p0, p1 = np.pi ** -0.25 * np.exp(-xs * xs / 2), np.pi ** -0.25 * np.exp(-xs * xs / 2) * np.sqrt(2) * xs
    vals = [p0]
    for n in range(1, 2 * kmax):
        p0, p1 = p1, np.sqrt(2 / (n + 1)) * xs * p1 - np.sqrt(n / (n + 1)) * p0
        if n % 2 == 1:
            vals.append(p1)
    ks = np.arange(kmax + 1)
    fac = np.pi ** 0.25 * np.exp(0.5 * (ks * math.log(4) + 2 * np.array([math.lgamma(k + 1) for k in ks])
                                        - np.array([math.lgamma(2 * k + 1) for k in ks])))
    return np.abs(np.array(vals)) * fac[:, None]  # row k: |F_k(x)| e^(x^2/2)


mp.mp.dps = 60
spot = []
for k, x in [(1, 0.5), (5, 3.0), (50, 10.0), (200, 7.0)]:
    ref = float(abs(mp.hyp1f1(k + 0.5, 0.5, -mp.mpf(x) ** 2)) * mp.exp(mp.mpf(x) ** 2 / 2))
    spot.append(dict(k=k, x=x, recurrence=float(hermite_sup(np.array([x]), k)[k, 0]), mpmath=ref))
F['hermite_spot_checks'] = spot
xs = np.round(np.arange(0.5, 20.0001, 0.005), 3)
H = hermite_sup(xs, 1200)
R = H / A_f(xs)[None, :]
F['ratio_max_k_ge_1'] = float(R[1:].max())
F['ratio_max_k_ge_0'] = float(R.max())
i, j = np.unravel_index(np.argmax(R[1:]), R[1:].shape)
F['ratio_argmax'] = dict(k=int(i) + 1, x=float(xs[j]))
F['ratio_max_note'] = float((H / A_f(xs, 1.0865)[None, :]).max())
F['sup_k_x_le_2'] = float(H[:, xs <= 2].max())
xs2 = np.round(np.arange(5, 30.0001, 0.01), 3)
H2 = hermite_sup(xs2, 4000)
sup, arg = H2.max(axis=0), H2.argmax(axis=0)
F['envelope_over_x13'] = dict(min=float((sup / xs2 ** (1 / 3)).min()), max=float((sup / xs2 ** (1 / 3)).max()),
                              k_over_x2_min=float((arg / xs2 ** 2).min()), k_over_x2_max=float((arg / xs2 ** 2).max()))
F['seconds'] = round(time.time() - T00, 1)
print('hermite', F['seconds'], flush=True)
json.dump(out, open('cm_orders.json', 'w'), indent=1)
