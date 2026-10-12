# jacobi_zeros.py: the ordinates from the Taylor coefficients at the centre, the paragraphs on the Stieltjes moments and
# Table tab:ch12:jacobi in Section ch:heat (ch/ch12.tex). Added for the book in October 2026, from the first author's note on
# the semicircle law and the Thorin framework (misc/reviews/nick-wigner/, Section 6 and Table 2 there), recomputed independently.
# The Stieltjes moments of Table tab:ch3:diag, m_n = sum_j a_j^{-n-1} = int_0^oo t^n W(t) dt/n!, are the power sums p_{n+1},
# and log G(u) = sum_j log(1 + u/a_j) = sum_{n>=1} (-1)^{n+1} p_n u^n/n gives them from the Taylor coefficients of
# G(u) = xi(1/2 + sqrt u)/xi(1/2) = Xi(i sqrt u)/Xi(0) = sum_k M_k u^k/((2k)! M_0), M_k = int_R v^{2k} Phi(v) dv, with Phi the
# kernel of eq:ch2:Phi. M_k by the trapezoid rule with step 1/600 on [0, 3.6] at 1000 digits: Phi is analytic in the strip
# |Im v| < pi/4, so the error is below 10^{-1200}, and the part beyond 3.6 is below 10^{-1700}. No zero enters the computation.
# With H = [m_{i+j}] and H' = [m_{i+j+1}], i, j < N, the Cholesky factorisation H = L L^T exists exactly when H is positive
# definite, and then J_N = L^{-1} H' L^{-T} is the Jacobi matrix of order N of the moment sequence. Under RH its eigenvalues are
# the nodes of the N point Gauss rule of the measure sum_j a_j^{-1} delta_{1/a_j}, and 1/sqrt(eigenvalue) approximates the
# ordinates from the bottom up. For N = 10, 20, ..., 50 the script records how many ordinates are recovered to within 10^{-6},
# the errors of some of them, log10 det H, the least diagonal entry of L, and the Jacobi parameters of order 50. It then finds
# the precision needed: the Taylor coefficients rounded to D digits and the whole computation repeated at D digits, the least
# D, in steps of 10, from which on (three steps in a row) the same ordinates are recovered to within 10^{-6}.
# Reference ordinates from Arb (acb.zeta_zeros). Run from any folder: python3 jacobi_zeros.py; writes jacobi_zeros.json there.
import json
import time

import mpmath as mp
from flint import acb, ctx

t0 = time.time()
out = {}


def rec(k, v):
    out[k] = v
    print(k, v, round(time.time() - t0, 1), flush=True)


def s(x, n=12):
    return mp.nstr(x, n)


mp.mp.dps = 1000
ctx.prec = 900
GAM = [mp.mpf(z.imag.mid().str(260, radius=False)) for z in acb.zeta_zeros(1, 30)]
ctx.prec = 53


def phi_nodes(h, V):
    # (v, Phi(v)) at v = j h, j = 0, ..., V/h, with Phi(v) = 2 sum_n (2 pi^2 n^4 e^{9v/2} - 3 pi n^2 e^{5v/2}) exp(-pi n^2 e^{2v})
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
        while True:   # qn = q^(n^2)
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
    # M_k = int_R v^{2k} Phi(v) dv = 2 int_0^oo v^{2k} Phi(v) dv, trapezoid rule, k = 0, ..., K
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
    # log(sum_k a_k u^k) = sum_n psi_n u^n with a_0 = 1: n psi_n = n a_n - sum_{0<j<n} j psi_j a_{n-j}; p_n = (-1)^{n+1} n psi_n
    K = len(a) - 1
    psi = [mp.mpf(0)] * (K + 1)
    for n in range(1, K + 1):
        psi[n] = (n * a[n] - mp.fsum(j * psi[j] * a[n - j] for j in range(1, n))) / n
    return [None] + [(-1) ** (n + 1) * n * psi[n] for n in range(1, K + 1)]


def jacobi(p, N):
    # H = [m_{i+j}], H' = [m_{i+j+1}] with m_i = p_{i+1}; None if the Cholesky factorisation fails
    H = mp.matrix(N, N)
    H1 = mp.matrix(N, N)
    for i in range(N):
        for j in range(N):
            H[i, j] = p[i + j + 1]
            H1[i, j] = p[i + j + 2]
    try:
        L = mp.cholesky(H)
    except ValueError:
        return None
    Li = mp.inverse(L)
    return Li * H1 * Li.T, L


def ordinates(J):
    return sorted(1 / mp.sqrt(x) for x in mp.eigsy(J, eigvals_only=True) if x > 0)


def count_close(zs, tol=mp.mpf('1e-6')):
    n = 0
    for k, z in enumerate(zs):
        if k < len(GAM) and abs(z - GAM[k]) < tol:
            n += 1
        else:
            break
    return n


h = mp.mpf(1) / 600
nodes = phi_nodes(h, mp.mpf('3.6'))
K = 101
M = moments(nodes, h, K)
xi_half = mp.mpf(-1) / 8 * mp.pi ** mp.mpf('-0.25') * mp.gamma(mp.mpf('0.25')) * mp.zeta(mp.mpf('0.5'))
rec('M_0 = Xi(0) against xi(1/2) = -(1/8) pi^(-1/4) Gamma(1/4) zeta(1/2)', [s(M[0], 30), s(xi_half, 30), s(abs(M[0] / xi_half - 1), 3)])
a = [M[k] / (mp.factorial(2 * k) * M[0]) for k in range(K + 1)]
p = power_sums(a)
rec('m_0, ..., m_3 (= p_1, ..., p_4)', [s(p[j], 15) for j in range(1, 5)])
full = {}
for N in (10, 20, 30, 40, 50):
    J, L = jacobi(p, N)
    zs = ordinates(J)
    cnt = count_close(zs)
    full[N] = cnt
    errs = {'gamma_%d' % (j + 1): s(abs(zs[j] - GAM[j]), 2) for j in sorted({0, 1, 4, 9, 14, 19, cnt - 1, cnt}) if j < min(len(zs), len(GAM))}
    d = {'ordinates within 1e-6': cnt, 'errors': errs,
         'log10 det H_N': s(mp.fsum(2 * mp.log10(L[i, i]) for i in range(N)), 6),
         'least diagonal entry of L': s(min(L[i, i] for i in range(N)), 4),
         'last diagonal entry of L, squared (det H_N / det H_{N-1})': s(L[N - 1, N - 1] ** 2, 4),
         'recovered ordinates beyond the count, first three': [s(z, 10) for z in zs[cnt:cnt + 3]]}
    if N == 50:
        d['alpha_1, ..., alpha_5 (Jacobi diagonal)'] = [s(J[i, i], 4) for i in range(5)]
        d['beta_1, ..., beta_5 (Jacobi off diagonal)'] = [s(J[i, i + 1], 4) for i in range(5)]
        d['max |J_ij| off the three diagonals'] = s(max(abs(J[i, j]) for i in range(N) for j in range(N) if abs(i - j) > 1), 3)
    rec('order %d' % N, d)


def trial(N, D):
    with mp.workdps(D):
        aD = [+x for x in a]   # the Taylor coefficients rounded to D digits
        r = jacobi(power_sums(aD), N)
        if r is None:
            return 0, None
        zs = ordinates(r[0])
        return count_close(zs), abs(zs[0] - GAM[0])


prec = {}
for N in (10, 20, 30, 40, 50):
    D = 20
    hist = []
    while D <= 1000:
        cnt, e1 = trial(N, D)
        hist.append((D, cnt, e1))
        if len(hist) >= 3 and all(c == full[N] for _, c, _ in hist[-3:]):
            break
        D += 10
    Dmin = hist[-3][0] if len(hist) >= 3 and all(c == full[N] for _, c, _ in hist[-3:]) else None
    prec[N] = {'least digits D': Dmin, 'error of gamma_1 at D': s(hist[-3][2], 3) if Dmin else None,
               'counts at D - 30, ..., D + 20': [c for d_, c, _ in hist if Dmin and Dmin - 30 <= d_ <= Dmin + 20]}
    rec('precision for order %d' % N, prec[N])
out['seconds'] = time.time() - t0
json.dump(out, open('jacobi_zeros.json', 'w'), indent=1)
print('done', round(time.time() - t0, 1))
