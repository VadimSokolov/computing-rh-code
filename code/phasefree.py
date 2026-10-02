# code/phasefree.py (September 2026, from misc/reviews/nick-audit/C-calibration): the reach of phase free bounds on the
# prime side of the explicit formula eq:ch12:Garith for the Gaussian smoothing G(u,t) of the zeros, \mathcal G(u,t) in
# the book (Section ch:lessons).
# Bounding |cos(u log n)| by one bounds the prime term by S(t)/sqrt(2 pi t), with
#   S(t) = sum_{n>=2} Lambda(n) n^{-1/2} e^{-(log n)^2/8t},
# so that the bound proves G(u,t) + G(-u,t) > 0 exactly when S(t) < sqrt(2 pi t) (pol(u,t) + A(u,t)), with the polar term
# pol(u,t) = 2 e^{t/2-2tu^2} cos 2tu and the archimedean term A(u,t) = (1/2pi) int_R e^{-2t(r-u)^2} m(r) dr. Since S increases
# with t and sqrt(2 pi t) A(u,t) is close to log(u/2pi)/2, this holds for t below a threshold t0(u), the root of
# S(t0) = sqrt(2 pi t0) (pol + A)(u,t0); also reported are the root of the approximation 2S(t) = log(u/2pi) and the ratio
# sqrt(2 pi t0) A(u,t0) / (log(u/2pi)/2). S(t) is the prime sum of eq:ch2:Warith at time 2t, so the explicit formula
# gives it from the archimedean integral:
#   S(t) = 2 sqrt(2 pi t) (e^{t/2} + (1/4pi) int_R e^{-2t r^2} m(r) dr - W(2t)),  m(r) = Re psi(1/4 + ir/2) - log pi,
# with W(2t) from the first 1700 zeros, below 1e-8 for t >= 0.05 (the zeros above 2197.26 add less than 1e-4000 there).
# The sieve of Lambda(n) up to 2e7, with the tail bounded by psi(x) < 1.03883 x (Rosser and Schoenfeld), checks S(t) for
# t <= 1.5. Reads zeros_1700.txt; writes phasefree.json.
# About a minute on one core.
import json

import mpmath as mp
import numpy as np

mp.mp.dps = 30
G2 = [mp.mpf(x) ** 2 for x in open('zeros_1700.txt').read().split()]
m = lambda r: mp.re(mp.digamma(mp.mpf(1) / 4 + 0.5j * r)) - mp.log(mp.pi)


def S_explicit(t):
    t = mp.mpf(t)
    R = 30 / mp.sqrt(t)
    pts = [mp.mpf(0)] + [mp.mpf(2) ** j for j in range(-2, 12) if mp.mpf(2) ** j < R] + [R]
    arch = 2 * mp.quad(lambda r: mp.exp(-2 * t * r * r) * m(r), pts) / (4 * mp.pi)
    W2t = mp.fsum(mp.exp(-2 * t * x) for x in G2)
    return 2 * mp.sqrt(2 * mp.pi * t) * (mp.exp(t / 2) + arch - W2t), W2t


def sieve_lambda(N):   # Lambda(n) for n <= N as a float array
    lam = np.zeros(N + 1)
    isp = np.ones(N + 1, bool); isp[:2] = False
    for p in range(2, int(N ** 0.5) + 1):
        if isp[p]:
            isp[p * p::p] = False
    for p in np.nonzero(isp)[0]:
        q = int(p)
        while q <= N:
            lam[q] = np.log(p); q *= int(p)
    return lam


N = 2 * 10 ** 7
lam = sieve_lambda(N)
n = np.nonzero(lam)[0]
ln, ll = lam[n], np.log(n.astype(float))


def S_sieve(t):
    s = float(np.sum(ln * np.exp(-0.5 * ll - ll ** 2 / (8 * t))))
    f = lambda x: x ** -0.5 * mp.exp(-mp.log(x) ** 2 / (8 * t))   # decreasing for x > 1
    L = mp.log(N)
    integral = mp.exp(t / 2) * mp.sqrt(2 * mp.pi * t) * mp.erfc((L - 2 * t) / mp.sqrt(8 * t))   # int_N^inf f
    psiN = float(np.sum(lam))
    tail = 1.03883 * (N * f(N) + integral) - psiN * f(N)
    return s, float(tail)


out = dict(S=[], t0=[])
for t in ['0.5', '1', '1.5', '1.72', '2', '3', '4', '6']:
    S, W2t = S_explicit(t)
    row = dict(t=float(t), S_explicit=mp.nstr(S, 15), W_2t=mp.nstr(W2t, 5), polar_model=mp.nstr(mp.sqrt(8 * mp.pi * mp.mpf(t)) * mp.exp(mp.mpf(t) / 2), 10))
    if float(t) <= 1.5:
        s, tail = S_sieve(float(t))
        row.update(S_sieve_2e7=s, sieve_tail_bound=tail)
    out['S'].append(row)
    print(row, flush=True)
def PA(us, t):   # pol(u,t) + A(u,t) and A(u,t), with the working precision raised by the digits of u
    with mp.workdps(mp.mp.dps + int(mp.log10(mp.mpf(us))) + 5):
        u, t = mp.mpf(us), mp.mpf(t)
        L = 12 / mp.sqrt(t)   # |r - u| > L adds less than 1e-120
        A = mp.quad(lambda r: mp.exp(-2 * t * (r - u) ** 2) * m(r), [u - L + 2 * L * i / 24 for i in range(25)]) / (2 * mp.pi)
        return 2 * mp.exp(t / 2 - 2 * t * u * u) * mp.cos(2 * t * u) + A, A


for u in ['185.24', '1977.22', '1e3', '1e6', '1e9', '3e12', '1e20', '1e50']:
    half_log = mp.log(mp.mpf(u) / (2 * mp.pi)) / 2
    t0 = mp.findroot(lambda t: S_explicit(t)[0] - mp.sqrt(2 * mp.pi * t) * PA(u, t)[0], (mp.mpf('0.2'), mp.mpf(5)), solver='anderson')
    t0_log = mp.findroot(lambda t: S_explicit(t)[0] - half_log, (mp.mpf('0.2'), mp.mpf(5)), solver='anderson')
    ratio = mp.sqrt(2 * mp.pi * t0) * PA(u, t0)[1] / half_log
    out['t0'].append(dict(u=u, t0=mp.nstr(t0, 8), t0_log_approximation=mp.nstr(t0_log, 8), A_ratio_at_t0=mp.nstr(ratio, 8),
                          two_loglog_u=mp.nstr(2 * mp.log(mp.log(mp.mpf(u))), 6)))
    print(out['t0'][-1], flush=True)
json.dump(out, open('phasefree.json', 'w'), indent=1)
