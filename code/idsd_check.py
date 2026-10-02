# idsd_check.py: the numerical inequalities in the proof of Theorem thm:ch12:idsd (Section ch:heat, ch/ch12.tex), that
# W > 0 and W' < 0 on (0, oo) when the zeros below height 1000 lie on the line. Added for the book in September 2026. The
# proof follows [P26, Theorem A] with the elementary constants of an independent check (misc/gpt-review/item130/).
# Each inequality of the proof is evaluated in ball arithmetic (python-flint, Arb) and counts only if it holds for the
# whole ball. The script also evaluates the constants of the cruder version of the proof (tau = 1e-6, height 10^4) and,
# as a cross check that is not part of the proof, the terms of the explicit formula at t = tau by mpmath quadrature.
# Run from any folder: python3 idsd_check.py; writes idsd_check.json there.
import json
import mpmath as mp
from flint import arb, acb, ctx
ctx.prec = 256
pi, gE, log2 = arb.pi(), arb.const_euler(), arb(2).log()
half = arb(1) / 2
out, holds = {}, []

def s(x, n=12):
    return x.str(n, radius=False) if isinstance(x, arb) else x

def check(name, cond, **vals):
    holds.append(bool(cond))
    out[name] = dict(holds=bool(cond), **{k: s(v) for k, v in vals.items()})

def m(r):  # m(r) = Re psi(1/4 + i r/2) - log pi, the function of the archimedean term
    return acb(arb(1) / 4, arb(r) / 2).digamma().real - pi.log()

def N_bounds(T):  # Trudgian's bound: |N(T) - (T/2pi) log(T/(2 pi e)) - 7/8| <= 0.112 log T + 0.278 log log T + 2.510 + 0.2/T, T >= e
    T = arb(T)
    main = T / (2 * pi) * (T / (2 * pi)).log() - T / (2 * pi) + arb(7) / 8
    err = arb('0.112') * T.log() + arb('0.278') * T.log().log() + arb('2.510') + arb('0.2') / T
    return main - err, main + err

# Small t: the archimedean term
tau = arb(1) / 20000
m0 = -gE - pi / 2 - 3 * log2 - pi.log()
check('m0', m0 > -6 and (m0 - m(0)).contains(0), m0=m0, psi_quarter_minus_log_pi=m(0))
binet = lambda r: arb(9) / (8 * arb(r) ** 2) + arb('1e-10')  # bound for |m(r) - log(r/(2 pi))| from Binet's formula, r >= 20
y = arb(10)
check('binet_pieces', (-(pi * y)).exp() * (1 / pi + 1 / (pi ** 2 * y)) / (1 - (-(pi * y)).exp()) < arb('1e-10')
      and 1 / (12 * (3 * y ** 2 / 4 - arb(1) / 16)) < 1 / (8 * y ** 2), note='tail of the Binet integral and its middle piece, at y = r/2 = 10')
check('binet_bound_r20', binet(20) < arb('0.003'), bound=binet(20))
spots = {}
for r in [20, 25, 30, 50, 100, 1000, 10000]:
    d = abs(m(r) - (arb(r) / (2 * pi)).log())
    spots[str(r)] = dict(diff=s(d), bound=s(binet(r)), holds=bool(d < binet(r)))
check('binet_spot_values', all(v['holds'] for v in spots.values()), values=spots)
check('m20', m(20) < 6 and ((arb(10) / pi).log() + arb('0.003')) < 6, m20=m(20))
I20 = 120 + 4 * pi + 20 * (arb(10) / pi).log() - 20
check('integral_0_20', I20 < 136, bound=I20)
alpha0 = -gE / 2 - (4 * pi).log()
alpha1 = alpha0 + 1
a0_psi = acb(half).digamma().real / 2 - (2 * pi).log()
a1_psi = acb(3 * half).digamma().real / 2 - (2 * pi).log()
check('alpha', (alpha0 - a0_psi).contains(0) and (alpha1 - a1_psi).contains(0), alpha0=alpha0, alpha1=alpha1)
E = [272 * arb(20) ** (2 * k) * tau ** (k + half) / (k + half).gamma() + arb('0.003') for k in (0, 1)]
check('E_at_tau', E[0] < arb('1.089') and E[1] < arb('0.047'), E0=E[0], E1=E[1])
ell = arb(20000).log() / 2
check('ell_tau', ell > arb('4.95'), ell=ell)

# Small t: the prime and polar terms. With y_n = (log n)^2/4t and q = log 2/(8t): e^{-y_n} and |y_n - 1/2| e^{-y_n} are at
# most e^{-y_n/2} (the maximum of (y + 1/2) e^{-y/2} is 2 e^{-3/4} < 1), e^{-y_n/2} <= n^{-q}, Lambda(n) <= log n < sqrt(n), and
# sum_{n>=2} n^{-q} <= 2^{-q} (1 + 2/(q - 1)) < 2 * 2^{-q} for q > 3, so |P^(k)(t)| < t^{-k-1/2} 2^{-q} = t^{-k-1/2} e^{-kappa/t}.
kappa = log2 ** 2 / 8
q = log2 / (8 * tau)
check('prime_constants', 2 * (arb(-3) / 4).exp() < 1 and q > 3 and 1 + 2 / (q - 1) < 2 and kappa / tau > 1201,
      q_tau=q, kappa_over_tau=kappa / tau)
R = [4 * pi / (k + half).gamma() * (-(kappa / tau)).exp() for k in (0, 1)]
check('prime_ratio', R[0] < arb('1e-500') and R[1] < arb('1e-500'), R0=R[0], R1=R[1])
polar = 2 * pi.sqrt() * tau ** (3 * half) * (tau / 4).exp()
check('polar_ratio', polar < arb('1e-5'), polar=polar)
K0 = -alpha0 + E[0] + R[0]
K1 = -alpha1 + E[1] + R[1] + polar
check('small_t_margins', K0 < arb('3.92') and K1 < arb('1.87') and ell - K0 > 0 and ell - K1 > 0,
      K0=K0, K1=K1, note='c0 >= s0 (ell - K0) and c1 >= s1 (ell - K1) on (0, tau]; ell >= ell(tau) there')

# Large t: the zeros
T = arb(1000)
H0 = (-(T ** 2 - arb(1) / 4) * tau).exp() * (T ** 2 + 1 / tau)
H1 = 2 * (-(T ** 2 - arb(1) / 4) * tau).exp() * (T ** 4 + 2 * T ** 2 / tau + 2 / tau ** 2)
L0 = (-2500 * tau).exp()
check('tail', H0 < arb('2e-16') and H1 < arb('4.1e-10') and L0 > arb('0.88') and 400 * L0 > 350
      and H0 < L0 and H1 < 400 * L0, H0=H0, H1=H1, L0=L0)
lo20, hi20 = N_bounds(20); lo50, hi50 = N_bounds(50)
lo100, hi100 = N_bounds(100); lo200, hi200 = N_bounds(200); lo1000, hi1000 = N_bounds(1000)
check('zero_counts', hi20 < 5 and lo50 > 6 and hi100 < lo200 and hi1000 < T ** 2,
      N20_upper=hi20, N50_lower=lo50, N100_upper=hi100, N200_lower=lo200, N1000_upper=hi1000,
      note='a zero in (20,50] and one in (100,200]; N(r) <= r^2 for r >= 1000 since the upper bound grows like r log r')

# The cruder version: tau' = 1e-6, height 10^4, only m(0) > -6 and m(60) > 1 needed
tc, Tc = arb('1e-6'), arb(10000)
m60_low = -6 + arb(57600) / 14401 + arb(577).log() / 2
Ac = 1 / (4 * (pi * tc).sqrt()) - 210 / pi
A1c = 1 / (8 * pi.sqrt() * tc ** (3 * half)) - 252000 / pi
Pc = 1 / (pi * tc).sqrt() * (-(log2 ** 2) / (4 * tc)).exp()
H0c = (-(Tc ** 2 - arb(1) / 4) * tc).exp() * (Tc ** 2 + 1 / tc)
H1c = 2 * (-(Tc ** 2 - arb(1) / 4) * tc).exp() * (Tc ** 4 + 2 * Tc ** 2 / tc + 2 / tc ** 2)
check('crude_version', m0 > -6 and m60_low > 1 and m(60) > m60_low and Ac > 55 and A1c > arb('1e7')
      and Pc < arb('0.01') and H0c < arb('2e-21') and H1c < arb('3e-13') and (-2500 * tc).exp() > arb('0.99'),
      m60_lower=m60_low, m60=m(60), A0_lower=Ac, A1_lower=A1c, P_upper=Pc, H0=H0c, H1=H1c)

# Cross-check, not part of the proof: the terms of the explicit formula at t = tau, 30 digits
mp.mp.dps = 30
t = mp.mpf(1) / 20000
mm = lambda r: mp.re(mp.digamma(mp.mpf(1) / 4 + 0.5j * r)) - mp.log(mp.pi)
cuts = [0, 1, 2, 5, 10, 20, 50, 100, 200, 400, 800, mp.inf]
A = mp.quad(lambda r: mp.exp(-t * r * r) * mm(r), cuts) / (2 * mp.pi)
A1 = mp.quad(lambda r: r * r * mp.exp(-t * r * r) * mm(r), cuts) / (2 * mp.pi)
lam = {}  # von Mangoldt function below 60; the terms with n >= 60 are below e^{-83000}
for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59]:
    k = p
    while k < 60:
        lam[k] = mp.log(p); k *= p
P = sum(lam[n] / mp.sqrt(n) * mp.exp(-mp.log(n) ** 2 / (4 * t)) for n in lam) / (2 * mp.sqrt(mp.pi * t))
s0, s1 = 1 / (4 * mp.sqrt(mp.pi * t)), mp.gamma(1.5) * t ** -1.5 / (4 * mp.pi)
Wt = mp.exp(t / 4) + A - P
out['crosscheck_at_tau'] = dict(A=mp.nstr(A, 15), A_over_s0=mp.nstr(A / s0, 10), minusAprime_over_s1=mp.nstr(A1 / s1, 10),
                                P=mp.nstr(P, 5), W=mp.nstr(Wt, 15), W_over_s0=mp.nstr(Wt / s0, 10),
                                minusWprime_over_s1=mp.nstr((A1 - mp.exp(t / 4) / 4) / s1, 10),
                                ell_plus_alpha0=mp.nstr(-mp.log(t) / 2 - mp.euler / 2 - mp.log(4 * mp.pi), 10),
                                note='not part of the proof: quadrature in floating point')
out['all_hold'] = all(holds)
json.dump(out, open('idsd_check.json', 'w'), indent=1)
print('all inequalities hold:', all(holds), '| checks:', len(holds))
for k, v in out.items():
    if isinstance(v, dict) and 'holds' in v and not v['holds']:
        print('FAILED', k, v)
print(json.dumps(out['crosscheck_at_tau'], indent=1))
