# jacobi_dbn.py: the Hankel test of the Stieltjes moments along the de Bruijn Newman flow, the paragraph at the end of
# Section sec:ch21:collisions (ch/ch21.tex). Added for the book in October 2026, from the first author's note on the semicircle
# law and the Thorin framework (misc/reviews/nick-wigner/, Section 7 there), recomputed independently. For H_t of Section
# sec:ch21:flow, H_t(i sqrt s) = sum_k s^k/(2k)! int_0^oo u^{2k} e^{tu^2} Phi_dBN(u) du, and with v = 2u the integral is
# 2^{-2k-3} M_k(t), M_k(t) = int_R v^{2k} e^{t v^2/4} Phi(v) dv, Phi the kernel of eq:ch2:Phi. So the power sums of the zeros z_j
# of H_t (z_j = 2 gamma_j at t = 0), their Hankel matrices and Jacobi matrices follow from M_k(t) as in jacobi_zeros.py
# (trapezoid rule with step 1/600 on [0, 3.6]; the node values of Phi are computed once and reweighted for each t).
# (1) For N = 10, 20, 30, 40, 50, at 100 + 16 N digits: positive definiteness of the order N Hankel matrix on the grid
#     t = 0, -0.05, ..., -4, the time t*_N at which it first fails going back from 0, by bisection to 10^{-7}, and the closest
#     adjacent pair among the ordinates that the order N Jacobi matrix returns just above t*_N; N = 30 again with step 1/550 at
#     650 digits.
# (2) The collisions of the real zeros. At 110 digits with step 1/100, the zeros of H_t that start at 2 gamma_j, j <= 64, are
#     followed back from t = 0 to t = -2.05 (dz/dt = H_zz/H_z, steps that keep each zero well inside its gaps, Newton's method
#     at each step). When two neighbours come within 0.05 of each other, their collision is located as the double zero
#     H_t = dH_t/dz = 0 by Newton's method in (z, t) (H_t = -H_zz holds exactly for the trapezoid sum), recorded and removed.
#     The top three zeros only follow; the collisions of the zeros that start below index 59 are reported and used in (3)
#     and (6) (dbn_collisions.py finds the same ones).
# (3) The pair behind each failure. At t = t*_N - 10^{-5} and at t = t*_N - delta (delta = 0.01, or half the time to the next
#     collision below t*_N if that is less), every pair that has collided by then is followed into the complex plane (from
#     z_c + i sqrt(2 (t_c - t)) near the collision, then by continuation with complex Newton steps), refined with the nodes and
#     the digits of the Hankel test, and its contribution 2 Re z^{-2n} removed from the power sums; the script records
#     whether the order N matrix is then positive definite, for each pair, for each pair alone (the contributions of all the
#     other pairs removed) and for all of them together.
# (4) The order 30 Hankel matrix at six times between the collisions of the pairs that start at (2 gamma_19, 2 gamma_20)
#     and (2 gamma_13, 2 gamma_14).
# (5) The zeros of H_t from the order 30 Jacobi matrix at t = 0.2, 0, -0.4 and -1.0, against Newton's method on H_t, and at
#     t = 0 against gamma_n.
# (6) The failure curve: the first order n_fail(t) <= 50 at which a pivot of the Hankel matrix is nonpositive (one elimination
#     of the order 50 matrix at 900 digits gives it), on the grid t = 0, -0.01, ..., -1.30, so that t*_N = sup{t : n_fail(t) <= N};
#     at one time on each plateau of n_fail, the tests of (3): n_fail again after the contribution of each complex pair is
#     removed, with each pair alone, and with all of them removed.
# Reference ordinates from Arb (acb.zeta_zeros). Run from any folder: python3 jacobi_dbn.py; writes jacobi_dbn.json there.
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


ctx.prec = 300
with mp.workdps(80):
    GAM = [mp.mpf(z.imag.mid().str(80, radius=False)) for z in acb.zeta_zeros(1, 60)]
ctx.prec = 53


def phi_nodes(h, V, dps):
    # (v, Phi(v)) at v = j h, j = 0, ..., V/h, with Phi(v) = 2 sum_n (2 pi^2 n^4 e^{9v/2} - 3 pi n^2 e^{5v/2}) exp(-pi n^2 e^{2v})
    nodes = []
    with mp.workdps(dps):
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


def power_sums(a):
    # log(sum_k a_k u^k) = sum_n psi_n u^n with a_0 = 1: n psi_n = n a_n - sum_{0<j<n} j psi_j a_{n-j}; p_n = (-1)^{n+1} n psi_n
    K = len(a) - 1
    psi = [mp.mpf(0)] * (K + 1)
    for n in range(1, K + 1):
        psi[n] = (n * a[n] - mp.fsum(j * psi[j] * a[n - j] for j in range(1, n))) / n
    return [None] + [(-1) ** (n + 1) * n * psi[n] for n in range(1, K + 1)]


def jacobi(p, N):
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


def flow_power_sums(nodes, h, t, K, dps):
    # power sums of the zeros of H_t from M_k(t), k = 0, ..., K
    with mp.workdps(dps):
        t = mp.mpf(t)
        m = [mp.mpf(0)] * (K + 1)
        for v, ph in nodes:
            w = 2 * h * ph * mp.exp(t * v * v / 4)
            if v == 0:
                m[0] += w / 2
                continue
            vv = v * v
            pw = mp.mpf(1)
            for k in range(K + 1):
                m[k] += w * pw
                pw *= vv
        a = [m[k] / (mp.factorial(2 * k) * m[0]) / mp.mpf(4) ** k for k in range(K + 1)]
        return power_sums(a)


def hankel_ok(nodes, h, t, N, dps):
    with mp.workdps(dps):
        return jacobi(flow_power_sums(nodes, h, t, 2 * N + 1, dps), N) is not None


def flow_ordinates(nodes, h, t, N, dps):
    with mp.workdps(dps):
        r = jacobi(flow_power_sums(nodes, h, t, 2 * N + 1, dps), N)
        return None if r is None else ordinates(r[0])


def bisect(nodes, h, N, dps, hi, lo, tol=mp.mpf('1e-7')):
    # hi: positive definite, lo: not
    while hi - lo > tol:
        mid = (hi + lo) / 2
        if hankel_ok(nodes, h, mid, N, dps):
            hi = mid
        else:
            lo = mid
    return hi, lo


rec('2 gamma_j for j = 12, 13, 14, 15 and 18, 19, 20, 21', [s(2 * GAM[j - 1], 9) for j in (12, 13, 14, 15, 18, 19, 20, 21)])
# (1) the order N Hankel test along the flow
h = mp.mpf(1) / 600
nodes = phi_nodes(h, mp.mpf('3.6'), 920)
tstar = {}
for N in (10, 20, 30, 40, 50):
    dps = 100 + 16 * N
    prev = mp.mpf(0)
    if not hankel_ok(nodes, h, prev, N, dps):
        tstar[N] = 'not positive definite at t = 0'
        rec('order %d' % N, tstar[N])
        continue
    fail = None
    for i in range(1, 81):
        t = -mp.mpf(i) / 20
        if not hankel_ok(nodes, h, t, N, dps):
            fail = t
            break
        prev = t
    if fail is None:
        tstar[N] = {'digits': dps, 'positive definite on the grid down to': s(prev, 4)}
        rec('order %d' % N, tstar[N])
        continue
    hi, lo = bisect(nodes, h, N, dps, prev, fail)
    zs = flow_ordinates(nodes, h, hi, N, dps)[:N // 2]
    gaps = [(zs[j + 1] - zs[j], j) for j in range(len(zs) - 1)]
    g, j = min(gaps)
    tstar[N] = {'digits': dps, 'first failure on the grid': s(fail, 4), 't*_N bracket (positive definite, not)': [s(hi, 10), s(lo, 10)],
                'closest adjacent ordinates just above t*_N (indices from 1, values, gap)': [j + 1, j + 2, s(zs[j], 10), s(zs[j + 1], 10), s(g, 3)]}
    rec('order %d' % N, tstar[N])
# N = 30 again, step 1/550 and 650 digits
h6 = mp.mpf(1) / 550
nodes6 = phi_nodes(h6, mp.mpf('3.6'), 670)
if isinstance(tstar.get(30), dict) and 't*_N bracket (positive definite, not)' in tstar[30]:
    hi30, lo30 = [mp.mpf(x) for x in tstar[30]['t*_N bracket (positive definite, not)']]
    hi6, lo6 = bisect(nodes6, h6, 30, 650, hi30 + mp.mpf('1e-4'), lo30 - mp.mpf('1e-4'))
    rec('order 30, step 1/550 and 650 digits: t*_30 bracket', [s(hi6, 10), s(lo6, 10)])

# (2) the collisions: the real zeros of H_t that start at 2 gamma_j, j <= 64, tracked back to t = -2.05 at 110 digits
mp.mp.dps = 110
KT, KR = 64, 58
ctx.prec = 400
GAMT = [mp.mpf(z.imag.mid().str(110, radius=False)) for z in acb.zeta_zeros(1, KT)]
ctx.prec = 53
hD = mp.mpf(1) / 100
nodesD = phi_nodes(hD, mp.mpf('3.6'), 130)
UD = [v / 2 for v, _ in nodesD]
UD2 = [u * u for u in UD]
UD3 = [u * u * u for u in UD]


def weights(t):
    # trapezoid weights of H_t(z) = (1/4) int_0^oo e^{t v^2/4} Phi(v) cos(z v/2) dv, step 1/100
    t = mp.mpf(t)
    return [hD * ph * mp.exp(t * v * v / 4) / 4 * (mp.mpf(1) / 2 if v == 0 else 1) for v, ph in nodesD]


def parts(z, W):
    # H, H_z, H_zz, H_zzz at real z for the weights W of one t; H_t = -H_zz and H_zt = -H_zzz (backward heat equation)
    x = z * hD / 2
    c1, s1 = mp.cos(x), mp.sin(x)
    c, sn = mp.mpf(1), mp.mpf(0)
    H = Hz = Hzz = Hzzz = mp.mpf(0)
    for j in range(len(W)):
        wc = W[j] * c
        ws = W[j] * sn
        H += wc
        Hz -= UD[j] * ws
        Hzz -= UD2[j] * wc
        Hzzz += UD3[j] * ws
        c, sn = c * c1 - sn * s1, sn * c1 + c * s1
    return H, Hz, Hzz, Hzzz


def newton(z, W, tol=mp.mpf(10) ** -30, maxit=40):
    P = None
    for _ in range(maxit):
        P = parts(z, W)
        dz = P[0] / P[1]
        z -= dz
        if abs(dz) < tol * (1 + abs(z)):
            return z, True, P
    return z, False, P


def double_zero(z, t, maxit=40):
    # Newton's method for H_t(z) = dH_t/dz(z) = 0 in (z, t)
    for _ in range(maxit):
        H, Hz, Hzz, Hzzz = parts(z, weights(t))
        a, b, c, d = Hz, -Hzz, Hzz, -Hzzz
        det = a * d - b * c
        dz = (d * H - b * Hz) / det
        dt = (-c * H + a * Hz) / det
        z -= dz
        t -= dt
        if abs(dz) < mp.mpf(10) ** -25 and abs(dt) < mp.mpf(10) ** -25:
            return z, t, parts(z, weights(t)), True
    return z, t, None, False


t = mp.mpf(0)
W = weights(t)
Z = [newton(2 * g, W)[0] for g in GAMT]
ids = list(range(1, KT + 1))
rec('(2) largest |z_j(0) - 2 gamma_j|, j <= %d' % KT, s(max(abs(Z[k] - 2 * GAMT[k]) for k in range(KT)), 3))
TEND = mp.mpf('-2.05')
GTH = mp.mpf('0.05')
cols = []
lost = []
steps = 0
P = [parts(z, W) for z in Z]
while t > TEND and len(Z) > 6:
    # one step back in time: dz/dt = H_zz/H_z; the top three zeros do not control the step and are dropped when they fail
    steps += 1
    V = [p_[2] / p_[1] for p_ in P]
    ctrl = len(Z) - 3
    close = [(Z[k + 1] - Z[k], k) for k in range(ctrl) if Z[k + 1] - Z[k] < GTH and V[k + 1] - V[k] > 0]
    tn = None
    if close:
        g, k = min(close)
        zc, tc, Pc, ok = double_zero((Z[k] + Z[k + 1]) / 2, t - g * g / 8)
        if ok and Z[k] - g < zc < Z[k + 1] + g and t - mp.mpf('0.01') < tc < t:
            i1, i2 = ids[k], ids[k + 1]
            cols.append({'ids': [i1, i2], 'zc': zc, 'tc': tc})
            rec('(2) collision %d' % len(cols), {'pair (indices at t = 0)': [i1, i2], '2 gamma at t = 0': [s(2 * GAMT[i1 - 1], 9), s(2 * GAMT[i2 - 1], 9)],
                                                  'double zero (z, t)': [s(zc, 12), s(tc, 12)], 'theta = z/2': s(zc / 2, 9), 'H, H_z there': [s(Pc[0], 3), s(Pc[1], 3)],
                                                  'isolated pair estimate -D_0^2/8': s(-(2 * GAMT[i2 - 1] - 2 * GAMT[i1 - 1]) ** 2 / 8, 5)})
            del Z[k:k + 2]
            del ids[k:k + 2]
            del P[k:k + 2]
            V = [p_[2] / p_[1] for p_ in P]
            tn = tc - mp.mpf('1e-9')
            ctrl = len(Z) - 3
        else:
            tn = t - g * g / 32
    if tn is None:
        dt = mp.mpf('0.05')
        for k in range(ctrl):
            r = V[k + 1] - V[k]
            if r > 0:
                dt = min(dt, mp.mpf('0.3') * (Z[k + 1] - Z[k]) / r)
        for k in range(ctrl):
            gl = Z[k] - Z[k - 1] if k > 0 else Z[0]
            if V[k] != 0:
                dt = min(dt, mp.mpf('0.2') * min(gl, Z[k + 1] - Z[k]) / abs(V[k]))
        tn = t - dt
    for attempt in range(40):
        Wn = weights(tn)
        newZ, newP, bad = [], [], None
        for k, z in enumerate(Z):
            zz, ok, Pz = newton(z + V[k] * (tn - t), Wn)
            if not ok:
                bad = k
                break
            newZ.append(zz)
            newP.append(Pz)
        if bad is None:
            for k in range(len(newZ) - 1):
                if not (newZ[k + 1] - newZ[k] > mp.mpf('1e-12')) or (newP[k][1] > 0) == (newP[k + 1][1] > 0):
                    bad = k + 1
                    break
        if bad is None:
            break
        if bad >= ctrl:
            lost.append({'ids': ids[bad:], 't': s(tn, 8)})
            Z, ids, P, V = Z[:bad], ids[:bad], P[:bad], V[:bad]
            ctrl = len(Z) - 3
            continue
        tn = t + (tn - t) / 2
    else:
        rec('(2) tracking failed at t', s(t, 10))
        break
    t, Z, P = tn, newZ, newP
low = [c_ for c_ in cols if max(c_['ids']) <= KR]
rec('(2) summary', {'steps': steps, 'final t': s(t, 6), 'tracked zeros left': len(Z), 'dropped at the top': lost,
                    'collisions with both indices <= %d (pair, z, theta, t)' % KR: [[c_['ids'], s(c_['zc'], 9), s(c_['zc'] / 2, 8), s(c_['tc'], 9)] for c_ in low]})


# (3) the pair behind each failure: remove the contribution 2 Re z^{-2n} of one complex pair from the power sums
def cparts(z, Wl, us):
    # H, H_z, H_zz at complex z for trapezoid weights Wl at the nodes v = 2 us
    E1 = mp.exp(1j * z * (us[1] - us[0]))
    F1 = 1 / E1
    E = F = mp.mpc(1)
    H = Hz = Hzz = mp.mpc(0)
    for j in range(len(Wl)):
        c = (E + F) / 2
        sn = (E - F) / 2j
        H += Wl[j] * c
        Hz -= Wl[j] * us[j] * sn
        Hzz -= Wl[j] * us[j] ** 2 * c
        E *= E1
        F *= F1
    return H, Hz, Hzz


def cnewton(z, Wl, us, tol, maxit=60):
    for _ in range(maxit):
        H, Hz, Hzz = cparts(z, Wl, us)
        dz = H / Hz
        z -= dz
        if abs(dz) < tol * (1 + abs(z)):
            return z, True
    return z, False


def follow(zc, tc, targets):
    # the zero z_c + i sqrt(2 (t_c - t)) + ... of the pair born at (z_c, t_c), continued to the times in targets (descending)
    res = {}
    d0 = mp.mpf('1e-4')
    tt = tc - d0
    z, ok = cnewton(zc + 1j * mp.sqrt(2 * d0), weights(tt), UD, mp.mpf(10) ** -40)
    if not ok:
        return None
    for T in targets:
        while tt > T:
            dt = min(mp.mpf('0.5') * (tc - tt), mp.mpf('0.01'), tt - T)
            H, Hz, Hzz = cparts(z, weights(tt), UD)
            tt -= dt
            z, ok = cnewton(z - (Hzz / Hz) * dt, weights(tt), UD, mp.mpf(10) ** -40)
            if not ok:
                return None
        res[T] = z
    return res


def contrib(z, K):
    # the contribution 2 Re z^{-2n}, n = 1, ..., K, of the conjugate pair z, conj(z) to the power sums
    w = 1 / (z * z)
    c, r = [None], w
    for n in range(1, K + 1):
        c.append(2 * mp.re(r))
        r *= w
    return c


def chol_ok(p, N):
    Hm = mp.matrix(N, N)
    for i in range(N):
        for j in range(N):
            Hm[i, j] = p[i + j + 1]
    try:
        mp.cholesky(Hm)
        return True
    except ValueError:
        return False


U6 = [v / 2 for v, _ in nodes]
tests = {}
for N, d in tstar.items():
    if not isinstance(d, dict) or 't*_N bracket (positive definite, not)' not in d:
        continue
    ts = mp.mpf(d['t*_N bracket (positive definite, not)'][0])
    below = [c_['tc'] for c_ in low if c_['tc'] < ts]
    delta = min(mp.mpf('0.01'), (ts - max(below)) / 2) if below else mp.mpf('0.01')
    tests[N] = [ts - mp.mpf('1e-5'), ts - delta]
paths = {}
for c_ in low:
    tg = sorted({T for Ts in tests.values() for T in Ts if T < c_['tc']}, reverse=True)
    if tg:
        paths[tuple(c_['ids'])] = follow(c_['zc'], c_['tc'], tg)
for N in sorted(tests):
    dps = 100 + 16 * N
    for T in tests[N]:
        with mp.workdps(dps + 20):
            p = flow_power_sums(nodes, h, T, 2 * N + 1, dps + 20)
            W6 = [2 * h * ph * mp.exp(T * v * v / 4) * (mp.mpf(1) / 2 if v == 0 else 1) for v, ph in nodes]
            items = []
            for c_ in low:
                key = tuple(c_['ids'])
                if c_['tc'] <= T or paths.get(key) is None or T not in paths[key]:
                    continue
                z, ok = cnewton(mp.mpc(paths[key][T]), W6, U6, mp.mpf(10) ** -(dps - 30))
                items.append((key, c_, z, ok, contrib(z, 2 * N + 1)))
            tot = [None] + [mp.fsum(it[4][n] for it in items) for n in range(1, 2 * N + 2)]
            rows = []
            for key, c_, z, ok, cz in items:
                q = [None] + [p[n] - cz[n] for n in range(1, 2 * N + 2)]
                qo = [None] + [p[n] - tot[n] + cz[n] for n in range(1, 2 * N + 2)]
                rows.append({'pair': list(key), 'collision t': s(c_['tc'], 8), 'theta = z/2 at t': s(z / 2, 9), 'refined': ok,
                             'positive definite without it': chol_ok(q, N), 'positive definite with it alone': chol_ok(qo, N)})
            qa = [None] + [p[n] - tot[n] for n in range(1, 2 * N + 2)]
            rec('(3) order %d at t = %s' % (N, s(T, 10)), {'t*_N minus t': s(mp.mpf(tstar[N]['t*_N bracket (positive definite, not)'][0]) - T, 3),
                                                        'positive definite': chol_ok(p, N), 'complex pairs among the zeros that start below index %d' % (KR + 1): rows,
                                                        'positive definite without all of them': chol_ok(qa, N)})

# (4) the order 30 Hankel matrix between the collisions of (2 gamma_19, 2 gamma_20) and (2 gamma_13, 2 gamma_14)
between = {}
for tb in ('-1.276', '-1.278', '-1.280', '-1.282', '-1.284', '-1.2852'):
    between[tb] = hankel_ok(nodes, h, mp.mpf(tb), 30, 580)
rec('(4) order 30 Hankel matrix positive definite between the two collisions', between)

# (5) the zeros of H_t from J_30 at four times
zt = {}
for T in ('0.2', '0', '-0.4', '-1.0'):
    T_ = mp.mpf(T)
    zs = flow_ordinates(nodes, h, T_, 30, 580)
    Wt = weights(T_)
    direct = [newton(z, Wt)[0] for z in zs[:13]]
    zt[T] = {'theta_n = z_n/2 for n = 1, 2, 3, 5, 10 (J_30)': [s(zs[n - 1] / 2, 10) for n in (1, 2, 3, 5, 10)],
             'z_n for n = 1, 2, 3, 5, 10 (J_30)': [s(zs[n - 1], 10) for n in (1, 2, 3, 5, 10)],
             '|J_30 - Newton| in theta, n = 1, ..., 13': [s(abs(zs[n] - direct[n]) / 2, 2) for n in range(13)]}
    if T == '0':
        zt[T]['|theta_n - gamma_n|, n = 1, ..., 13'] = [s(abs(zs[n] / 2 - GAM[n]), 2) for n in range(13)]
rec('(5) zeros of H_t from J_30', zt)


# (6) the failure curve n_fail(t) and the pairs behind it
def first_bad_pivot(p, N):
    # index (from 1) of the first nonpositive pivot of Gaussian elimination on [p_{i+j+1}], i, j < N, or None
    A = [[p[i + j + 1] for j in range(N)] for i in range(N)]
    for k in range(N):
        if A[k][k] <= 0:
            return k + 1
        for i in range(k + 1, N):
            f = A[i][k] / A[k][k]
            for j in range(k + 1, N):
                A[i][j] -= f * A[k][j]
    return None


NS = 50
dpsS = 100 + 16 * NS
grid = [-mp.mpf(i) / 100 for i in range(131)]
nf, ps = {}, {}
for T in grid:
    with mp.workdps(dpsS + 20):
        ps[T] = flow_power_sums(nodes, h, T, 2 * NS + 1, dpsS + 20)
        nf[T] = first_bad_pivot(ps[T], NS)
plat = []
for T in grid:
    if plat and plat[-1][2] == nf[T]:
        plat[-1][1] = T
    else:
        plat.append([T, T, nf[T]])
rec('(6) n_fail(t) on t = 0, -0.01, ..., -1.30: plateaus (from, to, first failing order, None if none up to 50)', [[s(a, 3), s(b, 3), c] for a, b, c in plat])
reps = [grid[(grid.index(a) + grid.index(b)) // 2] for a, b, c in plat if c is not None]
paths6 = {}
for c_ in low:
    tg = sorted([T for T in reps if T < c_['tc']], reverse=True)
    if tg:
        paths6[tuple(c_['ids'])] = follow(c_['zc'], c_['tc'], tg)
for T in reps:
    with mp.workdps(dpsS + 20):
        p = ps[T]
        W6 = [2 * h * ph * mp.exp(T * v * v / 4) * (mp.mpf(1) / 2 if v == 0 else 1) for v, ph in nodes]
        items = []
        for c_ in low:
            key = tuple(c_['ids'])
            if c_['tc'] <= T or paths6.get(key) is None or T not in paths6[key]:
                continue
            z, ok = cnewton(mp.mpc(paths6[key][T]), W6, U6, mp.mpf(10) ** -(dpsS - 30))
            items.append((key, c_, z, ok, contrib(z, 2 * NS + 1)))
        tot = [None] + [mp.fsum(it[4][n] for it in items) for n in range(1, 2 * NS + 2)]
        rows = []
        for key, c_, z, ok, cz in items:
            q = [None] + [p[n] - cz[n] for n in range(1, 2 * NS + 2)]
            qo = [None] + [p[n] - tot[n] + cz[n] for n in range(1, 2 * NS + 2)]
            rows.append({'pair': list(key), 'collision t': s(c_['tc'], 8), 'theta = z/2 at t': s(z / 2, 9), 'refined': ok,
                         'n_fail without it': first_bad_pivot(q, NS), 'n_fail with it alone': first_bad_pivot(qo, NS)})
        qa = [None] + [p[n] - tot[n] for n in range(1, 2 * NS + 2)]
        rec('(6) at t = %s' % s(T, 3), {'n_fail': nf[T], 'complex pairs among the zeros that start below index %d' % (KR + 1): rows,
                                         'n_fail without all of them': first_bad_pivot(qa, NS)})
out['seconds'] = time.time() - t0
json.dump(out, open('jacobi_dbn.json', 'w'), indent=1)
print('done', round(time.time() - t0, 1))
