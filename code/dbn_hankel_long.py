# dbn_hankel_long.py: the Hankel test of the Stieltjes moments further back along the de Bruijn Newman flow, to t = -20: Table
# tab:ch21:nfaillong and the paragraph before it in Section sec:ch21:hankel (ch/ch21.tex), which illustrate Theorem
# thm:ch12:inertia (ch/ch12.tex). Added for the book in October 2026, from the triage of the first author's chapter on Pakes
# (misc/reviews/nick-pakes/triage-A.md, items A36 and A37 and Recommendation 6, whose check scripts A_flow.py, A_flow_n30.py and
# A_flow_removal.py it replaces); none of the first author's code is used.
# Conventions of Section sec:ch21:flow: 8 H_t(2 theta) = int_R Phi(v) e^{t v^2/4} cos(theta v) dv, Phi the kernel of eq:ch2:Phi,
# so that the zeros theta_j(t) of theta -> H_t(2 theta) are the ordinates gamma_j at t = 0; G_t(u) = H_t(2i sqrt u)/H_t(0) =
# sum_k M_k(t) u^k/((2k)! M_0(t)) with M_k(t) = int_R v^{2k} e^{t v^2/4} Phi(v) dv, and the logarithm of G_t gives the moments
# m_n(t) = sum_j a_j(t)^{-n-1}, a_j(t) = theta_j(t)^2 (Newton's identities); H_N(t) = [m_{i+j}(t)] and H'_N(t) = [m_{i+j+1}(t)],
# i, j < N. A conjugate pair theta, conj(theta) of zeros contributes 2 Re theta^{-2n-2} to m_n(t).
# (1) The real zeros below theta = 88. The trapezoid sum for 8 H_t(2 theta) with step 1/100 on |v| <= 3.3 is evaluated in ball
#     arithmetic (python-flint, 230 bits) on the grid theta = 10, 10.02, ..., 88 for t = 0, -0.05, ..., -20, and its sign at
#     each point is read from the ball (the number of points where the ball contains zero is recorded). The sign changes are
#     followed from one time to the next and labelled by their index k at t = 0, where they lie within 0.01 of gamma_k; between
#     two times a pair of neighbours may disappear (a collision) and zeros may leave at the top, since the zeros rise as t
#     decreases. Each collision is then located as a double zero, 8 H_t(2 theta) = d/dtheta 8 H_t(2 theta) = 0, by Newton's
#     method in (theta, t) at 100 digits on the same trapezoid sum (node values at 400 bits), for which the backward heat
#     equation d/dt = -(1/4) d^2/dtheta^2 holds exactly.
# (2) Each collided pair is followed into the complex plane from theta_c + i sqrt((t_c - t)/2) just below its collision time,
#     by Euler steps in t, each corrected by Newton's method (ball arithmetic at 256 bits on the midpoints, the same trapezoid
#     sum), to the times of the table.
# (3) At t = 0, -1, -2, -3, -4, -5, -6, -8, -10, -12, -15, -20: the numbers iota_N(t) and iota'_N(t) of negative eigenvalues of
#     H_N(t) and H'_N(t) for N <= 30, from the signs of the pivots of Gaussian elimination (Sylvester's law of inertia), in two
#     settings, the moments from the trapezoid rule on [0, 3.6] with step 1/300 at 600 digits and with step 1/450 at 900 digits
#     (mpmath); in the second setting also the numbers of negative eigenvalues of the diagonally scaled matrices of orders 10
#     and 30.
# (4) The removal tests at the same times, t < 0: each pair of (2) is refined by Newton's method on the trapezoid sum of the
#     setting at its working precision, and its contribution is removed from the m_n(t); the first failing orders of H_N(t) and
#     H'_N(t), N <= 30, with one pair removed, with all of them removed, and with all but one removed (each pair alone).
# (5) The rows of Table tab:ch21:nfaillong: the lowest conjugate pair, its position, the real zeros below it, the number of
#     pairs among the zeros that start below theta = 88 that have collided, the first failing orders, (iota_10, iota'_10) and
#     (iota_30, iota'_30), and the orders N at which the first, second, ... negative eigenvalue appears.
# Reference ordinates from Arb (acb.zeta_zeros). Run from any folder: python3 dbn_hankel_long.py; writes dbn_hankel_long.json
# there. About six minutes on one core; python3 dbn_hankel_long.py quick runs a coarse version of every step.
import json
import sys
import time

import mpmath as mp
from flint import acb, arb, arb_mat, ctx

QUICK = len(sys.argv) > 1 and sys.argv[1] == 'quick'
T0 = time.time()
out = {'quick run': QUICK}


def rec(k, v):
    out[k] = v
    txt = json.dumps(v)
    print(k, txt if len(txt) < 4000 else txt[:4000] + ' ...', round(time.time() - T0, 1), flush=True)
    json.dump(out, open('dbn_hankel_long.json', 'w'), indent=1)


def s(x, n=12):
    return mp.nstr(x, n)


TTAB = [-1, -2, -3, -4, -5, -6, -8, -10, -12, -15, -20]  # the times of (2) to (5)
NMAX = 30
SETTINGS = ((300, 600), (450, 900))  # (1/step, digits) of the moments
DTHI, DTI, TEND = 50, 20, -20  # theta step 1/50, t step 1/20
if QUICK:
    TTAB, NMAX, SETTINGS, DTHI, DTI, TEND = [-1, -2, -4], 12, ((150, 250), (200, 400)), 25, 10, -4

with mp.workdps(70):
    ctx.prec = 300
    GAM = [mp.mpf(z.imag.mid().str(70, radius=False)) for z in acb.zeta_zeros(1, 30)]
rec('gamma_1, ..., gamma_30 (Arb)', [s(g, 10) for g in GAM])

# Step (1): the real zeros below theta = 88, t = 0 to TEND
ctx.prec = 230
HV, JV = 100, 330  # nodes v_j = j/100, j <= 330
pi = arb.pi()


def phi_arb(v):
    # Phi(v) = 2 sum_n (2 pi^2 n^4 e^{9v/2} - 3 pi n^2 e^{5v/2}) exp(-pi n^2 e^{2v})
    e2, a9, a5 = (2 * v).exp(), (9 * v / 2).exp(), (5 * v / 2).exp()
    tot, n = arb(0), 1
    while True:
        a = pi * n * n * e2
        tot += (2 * pi * pi * n ** 4 * a9 - 3 * pi * n * n * a5) * (-a).exp()
        if n >= 3 and a > arb(0.75 * ctx.prec + 40):
            break
        n += 1
    return 2 * tot


ctx.prec = 400
BHI = [(1 if j == 0 else 2) * phi_arb(arb(j) / HV) / HV for j in range(JV + 1)]  # trapezoid weights of int_R = 2 int_0^oo
ctx.prec = 230
V = [arb(j) / HV for j in range(JV + 1)]
B = [arb(b.mid()) for b in BHI]
NTH = 78 * DTHI
TH = [10 + arb(i) / DTHI for i in range(NTH + 1)]
THF = [10 + i / DTHI for i in range(NTH + 1)]
CM = arb_mat(NTH + 1, JV + 1)
for i in range(NTH + 1):
    for j in range(JV + 1):
        CM[i, j] = (TH[i] * V[j]).cos()


def signs(t):
    # signs of the trapezoid sum for 8 H_t(2 theta) on the grid (0 where the ball contains zero)
    Wc = arb_mat(JV + 1, 1)
    for j in range(JV + 1):
        Wc[j, 0] = B[j] * (t * V[j] * V[j] / 4).exp()
    Fm = CM * Wc
    return [1 if Fm[i, 0] > 0 else (-1 if Fm[i, 0] < 0 else 0) for i in range(NTH + 1)]


def zeros_of(sg):
    # midpoints of the cells where the sign changes, between consecutive points with a determined sign
    idx = [i for i in range(len(sg)) if sg[i] != 0]
    return [(THF[a] + THF[b]) / 2 for a, b in zip(idx[:-1], idx[1:]) if sg[a] * sg[b] < 0]


sg = signs(arb(0))
undet = sg.count(0)
prev = zeros_of(sg)
labels = list(range(1, len(prev) + 1))
rec('(1) grid', {'theta step': 1 / DTHI, 't step': 1 / DTI, 'nodes': JV + 1, 'grid points': NTH + 1})
rec('(1) t = 0: zeros on the grid, largest |theta - gamma_k|', [len(prev), round(max(abs(prev[k] - float(GAM[k])) for k in range(len(prev))), 4)])
snap = {0: (list(prev), list(labels))}
events, warn = [], []
NT = -TEND * DTI
for m in range(1, NT + 1):
    t = -m / DTI
    sg = signs(arb(-m) / DTI)
    undet += sg.count(0)
    cur = zeros_of(sg)
    d = len(prev) - len(cur)
    if d < 0:
        warn.append(['more zeros than before', t, len(prev), len(cur)])
        prev, labels = cur, [None] * len(cur)
        continue
    ntop = sum(1 for z in prev if z > 87.7)
    best = None
    for e in range(min(d, ntop) + 1):
        if (d - e) % 2:
            continue
        c = (d - e) // 2
        base = prev[:len(prev) - e]
        nb = len(base)
        combos = [()] if c == 0 else ([(a,) for a in range(nb - 1)] if c == 1 else
                                      ([(a, b) for a in range(nb - 1) for b in range(a + 2, nb - 1)] if c == 2 else []))
        for cb in combos:
            gone = {i for a in cb for i in (a, a + 1)}
            rem = [z for i, z in enumerate(base) if i not in gone]
            cost = max([abs(x - y) for x, y in zip(rem, cur)] + [0])
            if best is None or cost < best[0]:
                best = (cost, e, cb)
    if best is None or best[0] > 0.3:
        warn.append(['matching', t, len(prev), len(cur), best])
        if best is None:
            break
    cost, e, cb = best
    base, lab = prev[:len(prev) - e], labels[:len(labels) - e]
    gone = {i for a in cb for i in (a, a + 1)}
    for a in cb:
        events.append({'labels': [lab[a], lab[a + 1]], 'theta at t_hi': [base[a], base[a + 1]], 't_hi': t + 1 / DTI, 't_lo': t})
    labels = [l_ for i, l_ in enumerate(lab) if i not in gone]
    prev = cur
    if m % DTI == 0:
        snap[-m // DTI] = (list(cur), list(labels))
rec('(1) tracking: warnings, undetermined signs', [warn, undet])

# the collisions as double zeros, Newton's method in (theta, t) at 100 digits
mp.mp.dps = 100
N1 = [(mp.mpf(j) / HV, mp.mpf(BHI[j].mid().str(115, radius=False))) for j in range(JV + 1)]
H1 = mp.mpf(1) / HV


def parts_real(th, t):
    # F, F', F'', F''' at real theta, F the trapezoid sum for 8 H_t(2 theta)
    c1, s1 = mp.cos(th * H1), mp.sin(th * H1)
    c, sn = mp.mpf(1), mp.mpf(0)
    F = F1 = F2 = F3 = mp.mpf(0)
    for v, b in N1:
        w = b * mp.exp(t * v * v / 4)
        F += w * c
        F1 -= w * v * sn
        F2 -= w * v * v * c
        F3 += w * v ** 3 * sn
        c, sn = c * c1 - sn * s1, sn * c1 + c * s1
    return F, F1, F2, F3


def newton_real(th, t):
    for _ in range(60):
        F, F1, _, _ = parts_real(th, t)
        d = F / F1
        th -= d
        if abs(d) < mp.mpf(10) ** -40:
            return th, True
    return th, False


def double_zero(th, t):
    # F = F' = 0 in (theta, t); dF/dt = -F''/4 and dF'/dt = -F'''/4
    for _ in range(80):
        F, F1, F2, F3 = parts_real(th, t)
        a, b, c, d = F1, -F2 / 4, F2, -F3 / 4
        det = a * d - b * c
        dth, dt = (d * F - b * F1) / det, (-c * F + a * F1) / det
        th, t = th - dth, t - dt
        if abs(dth) < mp.mpf(10) ** -40 and abs(dt) < mp.mpf(10) ** -40:
            return th, t, True
    return th, t, False


cols = []
for ev in events:
    th_hi = mp.mpf(ev['t_hi'])
    za, oka = newton_real(mp.mpf(ev['theta at t_hi'][0]), th_hi)
    zb, okb = newton_real(mp.mpf(ev['theta at t_hi'][1]), th_hi)
    tt, za_, zb_ = th_hi, za, zb
    # if the pair is still wide at t_hi, follow it down by Euler steps until it is narrow, then solve for the double zero
    for _ in range(400):
        if zb_ - za_ < mp.mpf('0.05'):
            break
        Pa, Pb = parts_real(za_, tt), parts_real(zb_, tt)
        dt = (zb_ - za_) ** 2 / 8
        ya, ok1 = newton_real(za_ - Pa[2] / (4 * Pa[1]) * dt, tt - dt)
        yb, ok2 = newton_real(zb_ - Pb[2] / (4 * Pb[1]) * dt, tt - dt)
        if not (ok1 and ok2) or not yb - ya > 0:
            break
        tt, za_, zb_ = tt - dt, ya, yb
    thc, tc, ok = double_zero((za_ + zb_) / 2, tt - (zb_ - za_) ** 2 / 2)
    valid = ok and oka and okb and za - mp.mpf('0.05') < thc < zb + mp.mpf('0.05') and mp.mpf(ev['t_lo']) - mp.mpf('0.01') < tc < th_hi
    F, F1, F2, _ = parts_real(thc, tc)
    cols.append({'labels': ev['labels'], 'theta_c': thc, 't_c': tc, 'valid': valid, 'grid bracket': [ev['t_lo'], ev['t_hi']],
                 'residuals': [s(F, 3), s(F1, 3)], "F''": s(F2, 3)})
cols.sort(key=lambda c_: -c_['t_c'])
rec('(1) collisions (pair labelled by the indices at t = 0, theta_c, t_c, valid, grid bracket of t)',
    [[c_['labels'], s(c_['theta_c'], 10), s(c_['t_c'], 10), c_['valid'], c_['grid bracket'], c_['residuals'], c_["F''"]] for c_ in cols])
rec('(1) real zeros at the table times (labels, theta to 0.01)', {str(T): {'labels': snap[T][1], 'theta': [round(x, 2) for x in snap[T][0]]}
                                                               for T in [0] + TTAB if T in snap})

# Step (2): the pairs in the complex plane
ctx.prec = 256
VA = [arb(j) / HV for j in range(JV + 1)]
VA2 = [v * v for v in VA]
BA = [arb(b.mid()) for b in BHI]
HA = arb(1) / HV


def wts(t):
    ta = arb(mp.nstr(t, 40))
    return [BA[j] * (ta * VA2[j] / 4).exp() for j in range(JV + 1)]


def parts_c(th, w):
    # F, F', F'' at complex theta for the weights w
    E1 = (acb(0, 1) * th * HA).exp()
    E1i = 1 / E1
    E, Ei = acb(1), acb(1)
    F = Fd = Fdd = acb(0)
    for j in range(JV + 1):
        cj, sj = (E + Ei) / 2, (E - Ei) / acb(0, 2)
        F += w[j] * cj
        Fd -= w[j] * VA[j] * sj
        Fdd -= w[j] * VA2[j] * cj
        E, Ei = E * E1, Ei * E1i
    return F, Fd, Fdd


def mid(z):
    return acb(z.real.mid(), z.imag.mid())


def cnewton(th, w, maxit=40):
    for _ in range(maxit):
        F, Fd, _ = parts_c(th, w)
        d = mid(F / Fd)
        th = mid(th) - d
        if abs(d).mid() < arb('1e-25'):
            return th, True
    return th, False


def follow(thc, tc, targets):
    # the zero theta_c + i sqrt((t_c - t)/2) + ... of the pair, continued to the times in targets (descending, below t_c)
    d0 = mp.mpf('1e-4')
    tt = tc - d0
    th, ok = cnewton(acb(arb(mp.nstr(thc, 40)), arb(mp.nstr(mp.sqrt(d0 / 2), 40))), wts(tt))
    if not ok:
        return None, 'start'
    res, steps = {}, 0
    for T in targets:
        T = mp.mpf(T)
        while tt > T:
            dt = min(mp.mpf('0.5') * (tc - tt), mp.mpf('0.05'), tt - T)
            while True:
                _, Fd, Fdd = parts_c(th, wts(tt))
                pred = mid(th - Fdd / (4 * Fd) * acb(arb(mp.nstr(dt, 40))))
                new, ok = cnewton(pred, wts(tt - dt))
                if ok and float(abs(new - pred).mid()) < 0.02 and new.imag > 0:
                    break
                dt /= 2
                if dt < mp.mpf('1e-8'):
                    return None, 'step'
            th, tt = new, tt - dt
            steps += 1
        res[float(T)] = th
    return res, steps


paths = {}
for c_ in cols:
    tg = [T for T in TTAB if T < c_['t_c']]
    if tg and c_['valid']:
        r, info = follow(c_['theta_c'], c_['t_c'], tg)
        paths[tuple(c_['labels'])] = r
        rec('(2) pair %s followed' % c_['labels'], {'steps': info} if r is None else
            {'steps': info, 'theta at the table times': {str(int(T)): [float(z.real.mid()), float(z.imag.mid())] for T, z in r.items()}})


# Step (3): the inertia of H_N(t) and H'_N(t), N <= NMAX
def phi_mp(v, dps):
    q = mp.exp(-mp.pi * mp.exp(2 * v))
    a9, a5 = mp.exp(9 * v / 2), mp.exp(5 * v / 2)
    tot, n = mp.mpf(0), 1
    while True:
        term = (2 * mp.pi ** 2 * n ** 4 * a9 - 3 * mp.pi * n ** 2 * a5) * q ** (n * n)
        tot += term
        if n > 2 and abs(term) < abs(tot) * mp.mpf(10) ** (-dps - 20):
            break
        n += 1
    return 2 * tot


def node_values(hinv, dps, vmax=mp.mpf('3.6')):
    mp.mp.dps = dps
    hh = mp.mpf(1) / hinv
    return [(j * hh, hh * phi_mp(j * hh, dps) * (1 if j == 0 else 2)) for j in range(int(vmax * hinv) + 1)]


def moments(nodes, t, K):
    # m_0(t), ..., m_{K-1}(t): Newton's identities on log G_t, G_t(u) = sum_k M_k u^k/((2k)! M_0)
    M = [mp.mpf(0)] * (K + 1)
    ta = mp.mpf(t)
    for v, w in nodes:
        ww = w * mp.exp(ta * v * v / 4)
        p, vv = mp.mpf(1), v * v
        for k in range(K + 1):
            M[k] += ww * p
            p *= vv
    g = [M[k] / (mp.factorial(2 * k) * M[0]) for k in range(K + 1)]
    lg = [mp.mpf(0)] * (K + 1)
    for k in range(1, K + 1):
        lg[k] = (k * g[k] - mp.fsum(j * lg[j] * g[k - j] for j in range(1, k))) / k
    return [(-1) ** (k - 1) * k * lg[k] for k in range(1, K + 1)]  # m_{k-1} = sum_j a_j^{-k}


def inertia(m, shift, n):
    # iota_1, ..., iota_n of [m_{i+j+shift}]: numbers of negative pivots of Gaussian elimination (None after a zero pivot)
    A = [[m[i + j + shift] for j in range(n)] for i in range(n)]
    seq, neg = [], 0
    for k in range(n):
        piv = A[k][k]
        if piv == 0:
            return seq + [None] * (n - k)
        neg += piv < 0
        seq.append(neg)
        for i in range(k + 1, n):
            f = A[i][k] / piv
            for j in range(k, n):
                A[i][j] -= f * A[k][j]
    return seq


def first_fail(seq):
    return next((i + 1 for i, x in enumerate(seq) if x is None or x > 0), None)


def eig_neg(m, shift, n):
    A = mp.matrix([[m[i + j + shift] for j in range(n)] for i in range(n)])
    D = [1 / mp.sqrt(abs(A[i, i])) for i in range(n)]
    E = mp.eigsy(mp.matrix([[D[i] * A[i, j] * D[j] for j in range(n)] for i in range(n)]), eigvals_only=True)
    return sum(1 for i in range(n) if E[i] < 0)


def make_FdF(nodes, t):
    ta = mp.mpf(t)
    ew = [(v, w * mp.exp(ta * v * v / 4)) for v, w in nodes]
    hh = nodes[1][0]

    def FdF(z):
        E = mp.exp(mp.j * z * hh)
        Ei = 1 / E
        Ej, Eij = mp.mpc(1), mp.mpc(1)
        S, D = mp.mpc(0), mp.mpc(0)
        for j, (v, w) in enumerate(ew):
            if j:
                Ej *= E
                Eij *= Ei
            S += w * (Ej + Eij)
            D += w * v * (Ej - Eij)
        return S / 2, -D / (2 * mp.j)  # F = sum w cos(z v), F' = -sum w v sin(z v)

    return FdF


def refine(FdF, z, dps):
    for it in range(1, 41):
        F, dF = FdF(z)
        d = F / dF
        z -= d
        if abs(d) < mp.mpf(10) ** (-(dps - 40)):
            return z, it
    return z, None


K = 2 * NMAX
times = [0] + TTAB
inert, removal, pairs_at = {}, {}, {}
for hinv, dps in SETTINGS:
    key = 'step 1/%d, %d digits' % (hinv, dps)
    nodes = node_values(hinv, dps)
    inert[key], removal[key] = {}, {}
    for T in times:
        m = moments(nodes, T, K)
        row = {'iota_N, N = 1..%d' % NMAX: inertia(m, 0, NMAX), "iota'_N, N = 1..%d" % NMAX: inertia(m, 1, NMAX)}
        row['first failing N'] = [first_fail(row['iota_N, N = 1..%d' % NMAX]), first_fail(row["iota'_N, N = 1..%d" % NMAX])]
        if (hinv, dps) == SETTINGS[-1]:
            row['negative eigenvalues of the scaled matrices, N = 10 and %d' % NMAX] = [[eig_neg(m, sh, n) for sh in (0, 1)] for n in (10, NMAX)]
        inert[key][str(T)] = row
        if T == 0:
            continue
        # (4) the pairs at time T, refined at this precision, and the removal tests
        FdF = make_FdF(nodes, T)
        prs = []
        for lab, r in paths.items():
            if r is None or float(T) not in r:
                continue
            z0 = r[float(T)]
            z, it = refine(FdF, mp.mpc(mp.mpf(z0.real.mid().str(40, radius=False)), mp.mpf(z0.imag.mid().str(40, radius=False))), dps)
            prs.append((lab, z, it))
        prs.sort(key=lambda x: x[1].real)
        C = {lab: [2 * mp.re(z ** (-2 * (n + 1))) for n in range(K)] for lab, z, it in prs}

        def minus(labs):
            return [m[n] - mp.fsum(C[l_][n] for l_ in labs) for n in range(K)]

        def ff(mm):
            return [first_fail(inertia(mm, 0, NMAX)), first_fail(inertia(mm, 1, NMAX))]

        allm = minus(list(C))
        removal[key][str(T)] = {
            'pairs (labels, theta, Newton iterations)': [[list(lab), s(z.real, 15) + ' +- ' + s(abs(z.imag), 15) + 'i', it] for lab, z, it in prs],
            'first failing N of H_N, H\'_N': row['first failing N'],
            'with one pair removed': {str(list(lab)): ff(minus([lab])) for lab, _, _ in prs},
            'with all of them removed (iota_N, iota\'_N)': [inertia(allm, 0, NMAX), inertia(allm, 1, NMAX)],
            'with all of them removed (first failing N)': ff(allm),
            'each pair alone': {str(list(lab)): ff(minus([o for o in C if o != lab])) for lab, _, _ in prs}}
        pairs_at.setdefault(key, {})[T] = prs
        print(key, T, 'done', round(time.time() - T0, 1), flush=True)
    rec('(3) inertia, ' + key, inert[key])
    rec('(4) removal tests, ' + key, removal[key])
k1, k2 = ['step 1/%d, %d digits' % x for x in SETTINGS]
same_in = all(inert[k1][str(T)]['iota_N, N = 1..%d' % NMAX] == inert[k2][str(T)]['iota_N, N = 1..%d' % NMAX] and
              inert[k1][str(T)]["iota'_N, N = 1..%d" % NMAX] == inert[k2][str(T)]["iota'_N, N = 1..%d" % NMAX] for T in times)
same_rm = all(json.dumps({a: b for a, b in removal[k1][str(T)].items() if not a.startswith('pairs')}) ==
              json.dumps({a: b for a, b in removal[k2][str(T)].items() if not a.startswith('pairs')}) for T in TTAB)
eig_ok = all(inert[k2][str(T)]['negative eigenvalues of the scaled matrices, N = 10 and %d' % NMAX] ==
             [[inert[k2][str(T)]['iota_N, N = 1..%d' % NMAX][n - 1], inert[k2][str(T)]["iota'_N, N = 1..%d" % NMAX][n - 1]] for n in (10, NMAX)] for T in times)
mp.mp.dps = 30
pos_dev = max([abs(z1 - z2) for T in TTAB for (l1, z1, _), (l2, z2, _) in zip(pairs_at[k1].get(T, []), pairs_at[k2].get(T, []))] + [mp.mpf(0)])
rec('(3), (4) agreement of the two settings: inertia, removal tests, eigenvalue counts = pivot counts, largest difference of the pair positions',
    [same_in, same_rm, eig_ok, s(pos_dev, 3)])

# Step (5): the rows of the table
rows = {}
for T in TTAB:
    prs = pairs_at[k2].get(T, [])
    reals = snap.get(T, ([], []))[0]
    inn = inert[k2][str(T)]
    sq, sq1 = inn['iota_N, N = 1..%d' % NMAX], inn["iota'_N, N = 1..%d" % NMAX]
    r = {'first failing N of H_N, H\'_N': inn['first failing N'],
         '(iota_10, iota\'_10)': [sq[9], sq1[9]] if NMAX >= 10 else None,
         '(iota_%d, iota\'_%d)' % (NMAX, NMAX): [sq[-1], sq1[-1]],  # iota_{N+1} - iota_N is 0 or 1 by interlacing
         'orders N at which iota_N reaches 1, 2, ...': [sq.index(j) + 1 for j in range(1, sq[-1] + 1)] if None not in sq else None,
         "orders N at which iota'_N reaches 1, 2, ...": [sq1.index(j) + 1 for j in range(1, sq1[-1] + 1)] if None not in sq1 else None,
         'collided pairs among the zeros that start below theta = 88': len(prs),
         'real zeros below theta = 88': len(reals)}
    if prs:
        lab, z, _ = prs[0]
        r.update({'lowest pair': list(lab), 'theta': s(z.real, 8) + ' +- ' + s(abs(z.imag), 6) + 'i',
                  'real zeros below it': sum(1 for x in reals if x < z.real),
                  'pairs below theta = 88': sum(1 for _, w, _ in prs if w.real < 88)})
    rows[str(T)] = r
rec('(5) table rows', rows)
out['seconds'] = round(time.time() - T0, 1)
json.dump(out, open('dbn_hankel_long.json', 'w'), indent=1)
print('done', out['seconds'])
