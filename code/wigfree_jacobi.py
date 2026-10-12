# wigfree_jacobi.py: the Jacobi parameters of the Stieltjes moment sequence and the local arcsine (Abel) prediction, for the
# subsection on the semicircle as a local model in the free probability section (misc/reviews/nick-wigner/merge/wigfree.tex,
# label sec:wigfree:jacobi). Added for the book in October 2026, from the first author's note on the semicircle law and the
# Thorin framework (misc/reviews/nick-wigner/source-v7/, version 7, Section 9 and its table), recomputed independently.
# The Stieltjes moments m_n = sum_j a_j^{-n-1} (Table tab:ch3:diag) are the power sums p_{n+1} of the a_j = gamma_j^2, and
# log G(u) = sum_n (-1)^{n+1} p_n u^n / n gives them from the Taylor coefficients M_k / ((2k)! M_0) of G at 0, with
# M_k = int_R v^{2k} Phi(v) dv by the trapezoid rule (step 1/600 on [0, 3.6], 1000 digits), as in jacobi_zeros.py. No zero
# enters. J_50 = L^{-1} H' L^{-T} with H = [m_{i+j}] = L L^T and H' = [m_{i+j+1}], i, j < 50; alpha_n = J[n-1, n-1] and
# beta_n = |J[n-1, n]| (n counted from 1). The local arcsine model: a Jacobi matrix whose parameters vary slowly has, near
# index n, the arcsine density of states on [alpha_n - 2 beta_n, alpha_n + 2 beta_n]; if every such interval touches 0, its
# right end L(n) = alpha_n + 2 beta_n is about 2 alpha_n and about 4 beta_n, and matching the count of eigenvalues above x with
# N(x^{-1/2}) is an Abel equation whose solution is n(L) = int_{T < L^{-1/2}} dN(T) / sqrt(1 - L T^2). With the smooth density
# dN = (2 pi)^{-1} log(T / 2 pi) dT on T > 2 pi this script solves n(L) = n for L = L_n (the note's convention) and
# n(L) = n - 1/2 (the midpoint convention, which the sum over integer indices suggests), and the cruder rule of version 3 of the
# note, alpha_n = (1/2) gamma_{2n/pi}^{-2}, with gamma_m the inverse of the smooth count (T / 2 pi) log(T / 2 pi e) + 7/8. It
# records the parameters, the three predictions, the ratios, their mean and standard deviation over 17 <= n <= 49, and the
# ratios beta_n^2 / (alpha_n alpha_{n+1}), which are 1/4 for constant parameters whose interval touches 0.
# The precision is checked by repeating the Cholesky step at 600 digits. Run from any folder: python3 wigfree_jacobi.py; writes
# wigfree_jacobi.json there (a few minutes on one core).
import json
import time

import mpmath as mp

t0 = time.time()
out = {}


def rec(k, v):
    out[k] = v
    print(k, v, round(time.time() - t0, 1), flush=True)


def s(x, n=12):
    return mp.nstr(x, n)


mp.mp.dps = 1000


def phi_nodes(h, V):
    # (v, Phi(v)) at v = j h, Phi(v) = 2 sum_n (2 pi^2 n^4 e^{9v/2} - 3 pi n^2 e^{5v/2}) exp(-pi n^2 e^{2v}) (eq:ch2:Phi)
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


def moments(nodes, h, K):
    # M_k = int_R v^{2k} Phi(v) dv = 2 int_0^oo v^{2k} Phi(v) dv by the trapezoid rule, k = 0, ..., K
    m = [mp.mpf(0)] * (K + 1)
    for v, ph in nodes:
        w = 2 * h * ph
        if v == 0:
            m[0] += w / 2
            continue
        vv = v * v
        p = mp.mpf(1)
        for k in range(K + 1):
            m[k] += w * p
            p *= vv
    return m


def power_sums(a):
    # log(sum_k a_k u^k), a_0 = 1: n psi_n = n a_n - sum_{0<j<n} j psi_j a_{n-j}; p_n = (-1)^{n+1} n psi_n
    K = len(a) - 1
    psi = [mp.mpf(0)] * (K + 1)
    for n in range(1, K + 1):
        psi[n] = (n * a[n] - mp.fsum(j * psi[j] * a[n - j] for j in range(1, n))) / n
    return [None] + [(-1) ** (n + 1) * n * psi[n] for n in range(1, K + 1)]


h = mp.mpf(1) / 600
nodes = phi_nodes(h, mp.mpf('3.6'))
K = 101
M = moments(nodes, h, K)
xi_half = mp.mpf(-1) / 8 * mp.pi ** mp.mpf('-0.25') * mp.gamma(mp.mpf('0.25')) * mp.zeta(mp.mpf('0.5'))
rec('M_0 against xi(1/2), relative difference', s(abs(M[0] / xi_half - 1), 3))
a = [M[k] / (mp.factorial(2 * k) * M[0]) for k in range(K + 1)]
p = power_sums(a)
rec('p_1, ..., p_4', [s(p[j], 15) for j in range(1, 5)])
N = 50
H = mp.matrix(N, N)
H1 = mp.matrix(N, N)
for i in range(N):
    for j in range(N):
        H[i, j] = p[i + j + 1]
        H1[i, j] = p[i + j + 2]
L = mp.cholesky(H)
Li = mp.inverse(L)
J = Li * H1 * Li.T
alpha = [J[i, i] for i in range(N)]
beta = [abs(J[i, i + 1]) for i in range(N - 1)]
offtri = max(abs(J[i, j]) for i in range(N) for j in range(N) if abs(i - j) > 1)
rec('largest entry of J_50 off the three diagonals', s(offtri, 3))
rec('alpha_1..alpha_50', [s(x, 10) for x in alpha])
rec('beta_1..beta_49', [s(x, 10) for x in beta])
rec('alpha_n / (2 beta_n), n = 1..49', [s(alpha[i] / (2 * beta[i]), 4) for i in range(N - 1)])
rec('beta_n^2 / (alpha_n alpha_{n+1}), n = 1..49', [s(beta[i] ** 2 / (alpha[i] * alpha[i + 1]), 4) for i in range(N - 1)])

# precision: the power sums rounded to 600 digits and the Cholesky step repeated at 600 digits
mp.mp.dps = 600
H6 = mp.matrix(N, N)
H16 = mp.matrix(N, N)
for i in range(N):
    for j in range(N):
        H6[i, j] = +p[i + j + 1]
        H16[i, j] = +p[i + j + 2]
L6 = mp.cholesky(H6)
Li6 = mp.inverse(L6)
J6 = Li6 * H16 * Li6.T
mp.mp.dps = 1000
dal = max(abs(J6[i, i] / alpha[i] - 1) for i in range(N))
dbe = max(abs(abs(J6[i, i + 1]) / beta[i] - 1) for i in range(N - 1))
rec('largest relative change of alpha_n, beta_n when the Cholesky step runs at 600 instead of 1000 digits', [s(dal, 3), s(dbe, 3)])
rec('log10 det H_50 and least diagonal entry of L', [s(mp.log10(mp.fprod(L[i, i] for i in range(N)) ** 2), 8), s(min(L[i, i] for i in range(N)), 5)])

# the predictions, at 30 digits
mp.mp.dps = 30
al = [mp.mpf(x) for x in alpha]
be = [mp.mpf(x) for x in beta]


def n_of_L(L):
    # n(L) = int_{2 pi}^{T_s} (2 pi)^{-1} log(T / 2 pi) dT / sqrt(1 - L T^2), T_s = L^{-1/2}; T = T_s sin(phi)
    Ts = 1 / mp.sqrt(L)
    if Ts <= 2 * mp.pi:
        return mp.mpf(0)
    lo = mp.asin(2 * mp.pi / Ts)
    return mp.quad(lambda ph: mp.log(Ts * mp.sin(ph) / (2 * mp.pi)) * Ts / (2 * mp.pi), [lo, (lo + mp.pi / 2) / 2, mp.pi / 2])


def L_of_n(n):
    # solve n(L) = n in the variable log L by the Anderson Bjorck bracketing method: n(L) = 0 at L = 1/(4 pi^2), and n(L) > 49
    # at L = 10^{-7} (T_s = 3162)
    lo = mp.log(mp.mpf(10) ** -7)
    hi = -2 * mp.log(2 * mp.pi) - mp.mpf(10) ** -12
    return mp.e ** mp.findroot(lambda lg: n_of_L(mp.e ** lg) - n, (lo, hi), solver='anderson')


def gam_smooth(m):
    # the inverse of the smooth count (T / 2 pi) log(T / 2 pi e) + 7/8
    return mp.findroot(lambda T: T / (2 * mp.pi) * mp.log(T / (2 * mp.pi * mp.e)) + mp.mpf(7) / 8 - m, 2 * mp.pi * m / mp.log(m + 2) + 10)


pred = {}
for n in range(1, N):
    Ln = L_of_n(mp.mpf(n))
    Lh = L_of_n(mp.mpf(n) - mp.mpf(1) / 2) if n > 1 else None
    crude = 1 / (2 * gam_smooth(2 * mp.mpf(n) / mp.pi) ** 2)
    pred[n] = (Ln / 2, (Lh / 2) if Lh is not None else None, crude)
    # check of the Abel solution: n(L) at the solution
    if n in (5, 25, 49):
        rec('check n(L_n) - n at n = %d' % n, s(n_of_L(Ln) - n, 3))
rows = {}
for n in (5, 9, 13, 17, 21, 25, 29, 33, 37, 41, 45, 49):
    P, Ph, C = pred[n]
    rows[n] = {'alpha_n': s(al[n - 1], 4), '2 beta_n': s(2 * be[n - 1], 4), 'L_n/2 (Abel, n(L) = n)': s(P, 4),
               'alpha_n/(L_n/2)': s(al[n - 1] / P, 3), '2beta_n/(L_n/2)': s(2 * be[n - 1] / P, 3),
               'L/2 (Abel, n(L) = n - 1/2)': s(Ph, 4) if Ph is not None else None,
               'alpha_n/(L/2), midpoint': s(al[n - 1] / Ph, 3) if Ph is not None else None,
               '2beta_n/(L/2), midpoint': s(2 * be[n - 1] / Ph, 3) if Ph is not None else None,
               'crude (1/2) gamma_{2n/pi}^{-2}': s(C, 3), 'alpha_n/crude': s(al[n - 1] / C, 3), '2beta_n/crude': s(2 * be[n - 1] / C, 3)}
rec('table rows', rows)


def stats(xs):
    # mean, sample standard deviation, population standard deviation (the note's np.std), least and largest value
    m = mp.fsum(xs) / len(xs)
    ss = mp.fsum((x - m) ** 2 for x in xs)
    return [s(m, 4), s(mp.sqrt(ss / (len(xs) - 1)), 3), s(mp.sqrt(ss / len(xs)), 3), s(min(xs), 4), s(max(xs), 4)]


rng = range(17, 50)
rec('17 <= n <= 49, mean, sd, population sd, min, max of alpha_n/(L_n/2), note convention', stats([al[n - 1] / pred[n][0] for n in rng]))
rec('17 <= n <= 49, mean, sd, population sd, min, max of 2beta_n/(L_n/2), note convention', stats([2 * be[n - 1] / pred[n][0] for n in rng]))
rec('17 <= n <= 49, mean, sd, population sd, min, max of alpha_n/(L/2), midpoint convention', stats([al[n - 1] / pred[n][1] for n in rng]))
rec('17 <= n <= 49, mean, sd, population sd, min, max of 2beta_n/(L/2), midpoint convention', stats([2 * be[n - 1] / pred[n][1] for n in rng]))
rec('17 <= n <= 49, mean, sd, population sd, min, max of alpha_n/crude', stats([al[n - 1] / pred[n][2] for n in rng]))
rec('17 <= n <= 49, mean, sd, population sd, min, max of 2beta_n/crude', stats([2 * be[n - 1] / pred[n][2] for n in rng]))
rec('17 <= n <= 49, mean, sd, population sd, min, max of alpha_n/(2 beta_n)', stats([al[n - 1] / (2 * be[n - 1]) for n in rng]))
rec('17 <= n <= 49, mean, sd, population sd, min, max of beta_n^2/(alpha_n alpha_{n+1})', stats([be[n - 1] ** 2 / (al[n - 1] * al[n]) for n in rng]))
rec('1 <= n <= 49, mean, sd of beta_n^2/(alpha_n alpha_{n+1})', stats([be[n - 1] ** 2 / (al[n - 1] * al[n]) for n in range(1, 50)])[:2])
rec('crude over Abel (note convention) at n = 17, 33, 49', [s(pred[n][2] / pred[n][0], 4) for n in (17, 33, 49)])
out['series for the figure'] = {'n': list(range(1, N)), 'alpha': [s(al[n - 1], 8) for n in range(1, N)],
                                '2beta': [s(2 * be[n - 1], 8) for n in range(1, N)],
                                'abel': [s(pred[n][0], 8) for n in range(1, N)],
                                'abel_mid': [s(pred[n][1], 8) if pred[n][1] is not None else None for n in range(1, N)],
                                'crude': [s(pred[n][2], 8) for n in range(1, N)]}
out['seconds'] = round(time.time() - t0, 1)
json.dump(out, open('wigfree_jacobi.json', 'w'), indent=1)
print('done', out['seconds'])
