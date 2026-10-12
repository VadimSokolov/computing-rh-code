# biane_flow.py: the de Bruijn Newman flow of Polya's function against that of Xi, for Section sec:biane:flow of
# ch/biane.tex (Table tab:biane:flow, the collision parameters, the zero velocities and Figure fig:biane:flow). Added for the
# book in October 2026, from the first author's note on Biane's paper on Polya (misc/reviews/nick-biane/, Section 12 there),
# recomputed independently.
# The code writes lambda for the time of the flow (t in the book) and t for the frequency (theta in the book).
# For real lambda, P_lambda(t) = int_R e^{lambda u^2/4} e^{-2 pi cosh 2u} e^{itu} du and
# Xi_lambda(t) = int_R e^{lambda u^2/4} Phi(u) cos(tu) du, with Phi the kernel of eq:ch2:Phi, so P_0(t) = K_{it/2}(2 pi),
# Xi_0 = Xi, and Xi_lambda(t) = 8 H_lambda(2t) in the normalisation of Section sec:ch21:flow. Both by the trapezoid rule with
# step 1/250 on [0, 2] at 90 digits: the kernels are analytic in |Im u| < pi/4 and below 10^{-68} beyond u = 2, so for
# |t| <= 90, |Im t| <= 3 and -60 <= lambda <= 0 the error is far below the values, which are about e^{-pi t/4} at worst;
# the script checks P_0 against Arb's K_{it/2}(2 pi) and Xi_0 against mpmath's xi.
# Real zeros in (0, 90): sign changes on a grid of step 0.05, refined by the bracketing method of Anderson and Bjorck, for
# lambda = 0, -0.5, ..., -20 and -21, ..., -60. Zeros in the box 0.5 <= Re t <= 60.3, |Im t| <= 3: the argument principle
# on the boundary, step 0.04, halved wherever the argument changes by more than 0.5 between neighbours. Collisions of the
# pairs (gamma_9, gamma_10), (gamma_4, gamma_5), (gamma_2, gamma_3): bisection in lambda for the disappearance of the two
# sign changes in a window (step 0.01), then Newton's method for the double zero Xi_lambda = d/dt Xi_lambda = 0 in
# (t, lambda). Velocities of the first zeros at lambda = 0 from d t_k / d lambda = F_tt / (4 F_t) and from central
# differences with step 10^{-3}.
# Run from any folder: python3 biane_flow.py; writes biane_flow.json there (about half an hour on one core).
import json
import time

import mpmath as mp
from flint import acb, arb, ctx

t0 = time.time()
out = {}


def rec(k, v):
    out[k] = v
    print(k, v, round(time.time() - t0, 1), flush=True)


def s(x, n=12):
    return mp.nstr(x, n)


mp.mp.dps = 90
PI = mp.pi
H = mp.mpf(1) / 250
U = mp.mpf(2)
J = int(mp.nint(U / H))
UJ = [j * H for j in range(J + 1)]
WJ = [H / 2 if j in (0, J) else H for j in range(J + 1)]


def phi_book(u):
    q = mp.exp(-PI * mp.exp(2 * u))
    a9, a5 = mp.exp(mp.mpf(9) * u / 2), mp.exp(mp.mpf(5) * u / 2)
    tot = mp.mpf(0)
    for n in range(1, 40):
        term = (2 * PI ** 2 * n ** 4 * a9 - 3 * PI * n ** 2 * a5) * q ** (n * n)
        tot += term
        if n > 3 and abs(term) < mp.mpf(10) ** (-110):
            break
    return 2 * tot


BASE = {'P': [mp.exp(-2 * PI * mp.cosh(2 * u)) for u in UJ], 'X': [phi_book(u) for u in UJ]}


def weights(which, lam):
    lam = mp.mpf(lam)
    return [2 * w * b * mp.exp(lam * u * u / 4) for w, b, u in zip(WJ, BASE[which], UJ)]


def F(wv, t, d=0):
    # d-th t-derivative of int_R k(u) e^{lambda u^2/4} cos(tu) du, k even
    if d == 0:
        return mp.fdot(wv, [mp.cos(t * u) for u in UJ])
    if d == 1:
        return -mp.fdot(wv, [u * mp.sin(t * u) for u in UJ])
    if d == 2:
        return -mp.fdot(wv, [u * u * mp.cos(t * u) for u in UJ])
    return mp.fdot(wv, [u ** 3 * mp.sin(t * u) for u in UJ])


# ---------------------------------------------------------------- validation at lambda = 0
ctx.prec = 300
val = []
w0P, w0X = weights('P', 0), weights('X', 0)
for t in [10, 30, 60, 90]:
    ref = acb(2 * arb.pi()).bessel_k(acb(0, arb(t) / 2)).real
    ref = mp.mpf(ref.mid().str(60, radius=False))
    val.append(['P', t, s(F(w0P, mp.mpf(t)), 15), s(abs(F(w0P, mp.mpf(t)) / ref - 1), 3)])


def xi(sv):
    return sv * (sv - 1) / 2 * PI ** (-sv / 2) * mp.gamma(sv / 2) * mp.zeta(sv)


for t in [mp.mpf('14.134725'), mp.mpf(50), mp.mpf(90)]:
    ref = mp.re(xi(mp.mpf(1) / 2 + 1j * t))
    val.append(['Xi', s(t, 8), s(F(w0X, t), 15), s(abs(F(w0X, t) - ref) / abs(ref), 3)])
rec('validation at lambda = 0: function, t, quadrature value, relative error', val)
rec('kernels at u = 2 (truncation)', [s(BASE['P'][-1], 3), s(BASE['X'][-1], 3)])

# ---------------------------------------------------------------- real zeros in (0, 90)
TG = [mp.mpf(k) / 20 for k in range(1, 1801)]
CG = [[mp.cos(t * u) for u in UJ] for t in TG]
rec('grid prepared', len(TG))


def refine(wv, a, b):
    sc = max(abs(F(wv, a)), abs(F(wv, b)))
    return mp.findroot(lambda t: F(wv, t) / sc, (a, b), solver='anderson', tol=mp.mpf(10) ** -60)


def real_zeros(which, lam):
    wv = weights(which, lam)
    vals = [mp.fdot(wv, row) for row in CG]
    zs = []
    for k in range(len(TG) - 1):
        if vals[k] * vals[k + 1] < 0:
            zs.append(refine(wv, TG[k], TG[k + 1]))
    return zs, vals


LAMS = [mp.mpf(-k) / 2 for k in range(0, 41)] + [mp.mpf(-k) for k in range(21, 61)]
flow = {'P': {}, 'X': {}}
for which in ['P', 'X']:
    for lam in LAMS:
        zs, vals = real_zeros(which, lam)
        flow[which][s(lam, 6)] = [float(z) for z in zs]
    rec('flow computed for ' + which, len(LAMS))
out['flow_zeros'] = flow
tab = []
for lam in [0, -1, -2, -4, -6, -8, -10, -15, -20, -25, -30, -40, -50, -60]:
    key = s(mp.mpf(lam), 6)
    zp, zx = flow['P'][key], flow['X'][key]
    tab.append([lam, len(zp), round(zp[0], 4) if zp else None, len(zx), [round(z, 3) for z in zx[:4]]])
rec('lambda, real zeros of P_lambda in (0,90), first zero of P_lambda, real zeros of Xi_lambda in (0,90), first four', tab)
# monotone upward motion of the zeros of P_lambda: each zero increases as lambda decreases
mono = True
for i in range(len(LAMS) - 1):
    a, b = flow['P'][s(LAMS[i], 6)], flow['P'][s(LAMS[i + 1], 6)]
    m = min(len(a), len(b))
    if any(b[j] <= a[j] for j in range(m)) or len(b) > len(a):
        mono = False
rec('P_lambda: every zero increases as lambda decreases, and the count never increases', mono)


def min_gap(zs):
    best = None
    for a, b in zip(zs[:-1], zs[1:]):
        g = (b - a) * mp.log((a + b) / 2 / (2 * PI)) / (2 * PI)
        if best is None or g < best[0]:
            best = (g, a, b)
    return best


gaps = []
for lam in [0, -10, -30, -60]:
    key = s(mp.mpf(lam), 6)
    gp, gx = min_gap(flow['P'][key]), min_gap(flow['X'][key])
    gaps.append([lam, (s(gp[0], 4), round(gp[1], 3), round(gp[2], 3)) if gp else None,
                 (s(gx[0], 4), round(gx[1], 3), round(gx[2], 3)) if gx else None, len(flow['X'][key])])
rec('smallest gap in units of 2 pi/log(t/2pi): lambda, P (gap, pair), Xi (gap, pair), real zeros of Xi', gaps)

# ---------------------------------------------------------------- argument principle in the box
A0, B0, HB = mp.mpf('0.5'), mp.mpf('60.3'), mp.mpf(3)


def contour_pts(step):
    pts = []
    k = int(mp.nint((B0 - A0) / step))
    kk = int(mp.nint(2 * HB / step))
    pts += [A0 + (B0 - A0) * j / k - 1j * HB for j in range(k)]
    pts += [B0 + 1j * (-HB + 2 * HB * j / kk) for j in range(kk)]
    pts += [B0 - (B0 - A0) * j / k + 1j * HB for j in range(k)]
    pts += [A0 + 1j * (HB - 2 * HB * j / kk) for j in range(kk)]
    return pts


def winding(wv):
    pts = contour_pts(mp.mpf('0.04'))
    pts.append(pts[0])
    f = lambda t: mp.fdot(wv, [mp.cos(t * u) for u in UJ])
    vals = [f(t) for t in pts]
    tot, mx, nref = mp.mpf(0), mp.mpf(0), 0
    for i in range(len(pts) - 1):
        stack = [(pts[i], vals[i], pts[i + 1], vals[i + 1], 0)]
        while stack:
            p0, v0, p1, v1, dep = stack.pop()
            d = mp.arg(v1 / v0)
            if abs(d) > 0.5 and dep < 12:
                pm = (p0 + p1) / 2
                vm = f(pm)
                nref += 1
                stack.append((pm, vm, p1, v1, dep + 1))
                stack.append((p0, v0, pm, vm, dep + 1))
            else:
                tot += d
                mx = max(mx, abs(d))
    return tot / (2 * PI), mx, nref


ap = []
for which, lam in [('P', 0), ('P', -10), ('P', -30), ('P', -60), ('X', 0), ('X', -2), ('X', -4)]:
    wv = weights(which, lam)
    wnd, mx, nref = winding(wv)
    zs = flow[which][s(mp.mpf(lam), 6)]
    nreal = sum(1 for z in zs if A0 <= z <= B0)
    dist = min(abs(z - float(B0)) for z in zs) if zs else None
    ap.append([which, lam, s(wnd, 8), nreal, s(mx, 3), nref, round(dist, 3) if dist is not None else None])
rec('argument principle: function, lambda, winding number, real zeros in [0.5, 60.3], largest step, refinements, '
    'distance of the nearest real zero to Re t = 60.3', ap)

# ---------------------------------------------------------------- collisions of Xi_lambda
def nsign(wv, a, b, step=mp.mpf('0.01')):
    n = int(mp.nint((b - a) / step))
    vals = [F(wv, a + (b - a) * j / n) for j in range(n + 1)]
    idx = [j for j in range(n) if vals[j] * vals[j + 1] < 0]
    return len(idx), [a + (b - a) * (j + mp.mpf(1) / 2) / n for j in idx]


coll = []
for (lo_t, hi_t, lam_hi, lam_lo, label) in [(mp.mpf(47), mp.mpf('51.5'), mp.mpf(0), mp.mpf('-2.5'), '9,10'),
                                             (mp.mpf('29.5'), mp.mpf('35.5'), mp.mpf(-2), mp.mpf(-5), '4,5'),
                                             (mp.mpf(20), mp.mpf(31), mp.mpf(-10), mp.mpf(-15), '2,3')]:
    n_hi = nsign(weights('X', lam_hi), lo_t, hi_t)[0]
    n_lo = nsign(weights('X', lam_lo), lo_t, hi_t)[0]
    a, b = lam_lo, lam_hi
    for _ in range(14):
        m = (a + b) / 2
        if nsign(weights('X', m), lo_t, hi_t)[0] == n_hi:
            b = m
        else:
            a = m
    n_b, mids = nsign(weights('X', b), lo_t, hi_t)
    tc0 = (mids[0] + mids[-1]) / 2 if mids else (lo_t + hi_t) / 2

    def eqs(t, lam):
        wv = weights('X', lam)
        e = mp.exp(PI * t / 4)
        return [F(wv, t) * e, F(wv, t, 1) * e]

    sol = mp.findroot(eqs, (tc0, b), tol=mp.mpf(10) ** -40)
    coll.append([label, n_hi, n_lo, s(b, 8), s(sol[0], 10), s(sol[1], 10)])
rec('collisions of Xi_lambda: pair, sign changes before and after, bisection value of lambda, t_c, lambda_c', coll)

# ---------------------------------------------------------------- velocities at lambda = 0
vel = []
for which, g in [('P', mp.mpf('19.54')), ('X', mp.mpf('14.1347'))]:
    w0 = weights(which, 0)
    sc = lambda t: mp.exp(PI * t / 4)
    z0 = mp.findroot(lambda t: F(w0, t) * sc(t), g, tol=mp.mpf(10) ** -50)
    form = F(w0, z0, 2) / (4 * F(w0, z0, 1))
    ep = mp.mpf('0.001')
    wp, wm = weights(which, ep), weights(which, -ep)
    zp = mp.findroot(lambda t: F(wp, t) * sc(t), z0, tol=mp.mpf(10) ** -50)
    zm = mp.findroot(lambda t: F(wm, t) * sc(t), z0, tol=mp.mpf(10) ** -50)
    vel.append([which, s(z0, 12), s(form, 8), s((zp - zm) / (2 * ep), 8)])
rec('velocity of the first zero at lambda = 0: function, zero, F_tt/(4 F_t), central difference', vel)

out['seconds'] = round(time.time() - t0, 1)
json.dump(out, open('biane_flow.json', 'w'), indent=1)
print('done', out['seconds'])
