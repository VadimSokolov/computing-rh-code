# dbeta_minimum.py: the constants of Theorem thm:dbeta:minall and Table tab:dbeta:minall (ch/basepoint.tex, the comparison
# across vertical lines and the minimum at the origin below the basepoint one), the ratio table tab:dbeta:ratio, and the
# minima quoted in Section sec:dbeta:beta (ch/ch16.tex). Added for the book in October 2026, from the first author's note on
# the Dirichlet beta function in Thorin form (misc/reviews/nick-dbeta/, Section 10 there), recomputed independently with Arb
# (python-flint). Three functions: xi, xi(s, chi_4) = (4/pi)^{(s+1)/2} Gamma((s+1)/2) L(s, chi_4) and
# xi(s, chi_12) = (12/pi)^{s/2} Gamma(s/2) L(s, chi_12); v_sigma(t) = Re (f'/f)(sigma+it), b_sigma = v_sigma(0).
# For each: V = (log f)''(1/2) from the Taylor series of log f at 1/2; the elementary lower bound on the line 1+delta,
# v_{1+delta}(t) >= A(1+delta, t) - S(1+delta), with A(sigma, t) the gamma factor part (Re of the digamma term plus the
# conductor constant; for xi the terms 1/s + 1/(s-1), of positive real part, are dropped) and S(sigma) the sum of
# Lambda(n) n^{-sigma} over the n prime to the conductor (computed from zeta'/zeta); t_1 = the root of
# A(1+delta, t) - S(1+delta) = (1/2+delta) V, found by bisection, and certified at t_1 rounded up to three decimals in ball
# arithmetic (A increases in t, so the inequality holds for every larger t); the delta that minimises t_1; the margins of the
# convexity step, tau_1 - t_1 - (1/2+delta)/sqrt(3), with tau_1 the first ordinate (certified in dbeta_zeros.py).
# Also: b_sigma/(sigma - 1/2) against V; the ratio v_sigma(t)/(sigma - 1/2) on rows of sigma; grid minima of
# v_sigma(t) - b_sigma for sigma <= 1 and, for chi_4 and chi_12, for sigma beyond 1+delta (the open case), with
# v_sigma''(0) = -(log f)'''(sigma); and the point where the comparison bound at sigma = 0.75 meets b_{0.75} for chi_4.
# Further: the largest delta (step 0.05) keeping the convexity margin positive, certified as above; and, for chi_4 and
# chi_12, a proof in ball arithmetic that v_sigma(t) > b_sigma for every t != 0 and every sigma from some sigma_2 on, by
# bounding the prime part below by the terms with chi(n) = -1 alone (the comments at the step give the argument).
# Run on Hopper (a few minutes): python3 dbeta_minimum.py; writes dbeta_minimum.json in the current folder.
import json
import time

from flint import acb, acb_series, arb, ctx, dirichlet_char

t0 = time.time()
out = {}


def rec(k, v):
    out[k] = v
    print(k, v, round(time.time() - t0, 1), flush=True)


def S_(x, n=10):
    return x.str(n, radius=False)


ctx.prec = 128
CHI4 = dirichlet_char(4, 3)
CHI12 = dirichlet_char(12, 11)
PI = arb.pi()
LOGPI = PI.log()
LOG4PI = (arb(4) / PI).log()
LOG12PI = (arb(12) / PI).log()
LOG2 = arb(2).log()
LOG3 = arb(3).log()


def dlogL(s, chi):
    c = acb_series.dirichlet_l(acb_series([s, 1], prec=2), chi).coeffs()
    return c[1] / c[0]


def dlogzeta(s):
    c = acb_series.zeta(acb_series([s, 1], prec=2)).coeffs()
    return c[1] / c[0]


def LL_xi(s):
    if s == acb(1):
        # the pole of zeta'/zeta cancels 1/(s-1): (xi'/xi)(1) = -log(pi)/2 + psi(3/2)/2 + Euler's constant
        return acb(-LOGPI / 2 + arb(1.5).digamma() / 2 + arb.const_euler())
    return 1 / s + 1 / (s - 1) - LOGPI / 2 + (s / 2).digamma() / 2 + dlogzeta(s)


def LL_4(s):
    return LOG4PI / 2 + ((s + 1) / 2).digamma() / 2 + dlogL(s, CHI4)


def LL_12(s):
    return LOG12PI / 2 + (s / 2).digamma() / 2 + dlogL(s, CHI12)


def logser_xi(x, n):
    s = acb_series([x, 1], prec=n)
    return (s * (1 - s) / 2).log() - s / 2 * acb(LOGPI) + (s / 2).lgamma() + acb_series.zeta(s).log()


def logser_4(x, n):
    s = acb_series([x, 1], prec=n)
    return (s + 1) / 2 * acb(LOG4PI) + ((s + 1) / 2).lgamma() + acb_series.dirichlet_l(s, CHI4).log()


def logser_12(x, n):
    s = acb_series([x, 1], prec=n)
    return s / 2 * acb(LOG12PI) + (s / 2).lgamma() + acb_series.dirichlet_l(s, CHI12).log()


def A_xi(sg, t):
    return -LOGPI / 2 + (acb(sg, t) / 2).digamma().real / 2


def A_4(sg, t):
    return LOG4PI / 2 + ((acb(sg, t) + 1) / 2).digamma().real / 2


def A_12(sg, t):
    return LOG12PI / 2 + (acb(sg, t) / 2).digamma().real / 2


def S_xi(sg):
    return -dlogzeta(acb(sg)).real


def S_4(sg):
    return S_xi(sg) - LOG2 / (arb(2) ** arb(sg) - 1)


def S_12(sg):
    return S_4(sg) - LOG3 / (arb(3) ** arb(sg) - 1)


FUN = {
    'xi': dict(LL=LL_xi, logser=logser_xi, A=A_xi, S=S_xi, delta='3.75', tau1=14.134725141734693, H0='3e12'),
    'chi4': dict(LL=LL_4, logser=logser_4, A=A_4, S=S_4, delta='2.15', tau1=6.020948904697597, H0='2.5e7'),
    'chi12': dict(LL=LL_12, logser=logser_12, A=A_12, S=S_12, delta='1.30', tau1=3.804627614, H0='8.33e6'),
}


def t1_of(F, d, V, hi=60.0, tol=1e-10):
    # root in t of A(1+d, t) - S(1+d) - (1/2+d) V, by bisection (the function increases in t)
    sg = arb(1) + arb(d)
    c = F['S'](sg) + (arb(0.5) + arb(d)) * V
    g = lambda t: F['A'](sg, t) - c
    if g(arb(0)) > 0:
        return 0.0
    lo, hi_ = 0.0, hi
    while hi_ - lo > tol:
        m = (lo + hi_) / 2
        if g(arb(m)) > 0:
            hi_ = m
        else:
            lo = m
    return hi_


res = {}
for name, F in FUN.items():
    V = 2 * F['logser'](acb(0.5), 3).coeffs()[2].real
    F['V'] = V
    d = F['delta']
    t1 = t1_of(F, d, V)
    t1r = (int(t1 * 1000) + 1) / 1000
    sg = arb(1) + arb(d)
    margin = F['A'](sg, arb(t1r)) - F['S'](sg) - (arb(0.5) + arb(d)) * V
    # the delta minimising t_1, on a grid of step 0.05 and then 0.005
    ctx.prec = 64
    best = min(((t1_of(F, '%.2f' % (k * 0.05), V, tol=1e-6), k * 0.05) for k in range(6, 201)))
    best = min(((t1_of(F, '%.3f' % (best[1] + j * 0.005), V, tol=1e-8), best[1] + j * 0.005) for j in range(-10, 11)))
    ctx.prec = 128
    tau1 = F['tau1']
    rows = {}
    for a in ['0.55', '0.6', '0.75', '1', '1.5', '2', '3']:
        b = F['LL'](acb(arb(a))).real
        rows[a] = {'b': S_(b, 10), 'b/(a-1/2)': S_(b / (arb(a) - 0.5), 8), 'below V': bool(b / (arb(a) - 0.5) < V)}
    res[name] = {'V': S_(V, 12), 'delta': d, 't1': '%.6f' % t1, 't1 rounded up': t1r,
                 'certified A-S-(1/2+d)V at t1 rounded up': margin.str(5), 'positive': bool(margin > 0),
                 'optimal delta': '%.3f' % best[1], 't1 at optimal delta': '%.5f' % best[0],
                 'tau1': tau1, 'tau1 - t1 - 1/(2 sqrt3)': '%.4f' % (tau1 - t1r - 0.5 / 3 ** 0.5),
                 'tau1 - t1 - (1/2+delta)/sqrt3': '%.4f' % (tau1 - t1r - (0.5 + float(d)) / 3 ** 0.5),
                 'H0': F['H0'], 'b_alpha': rows}
    rec(name, res[name])

# the largest delta, in steps of 0.05, for which the convexity margin tau_1 - t_1 - (1/2+delta)/sqrt3 stays above 0.02,
# with t_1 and the margin then certified at 128 bits as above
dmax = {}
for name, F in FUN.items():
    V, tau1 = F['V'], F['tau1']
    ctx.prec = 64
    k = int(round(float(F['delta']) / 0.05))
    last = None
    while True:
        dd = '%.2f' % (k * 0.05)
        t1 = t1_of(F, dd, V, hi=200.0, tol=1e-7)
        m = tau1 - t1 - (0.5 + float(dd)) / 3 ** 0.5
        if m < 0.02:
            break
        last = (dd, t1)
        k += 1
    ctx.prec = 128
    dd, t1 = last
    t1r = (int(t1 * 1000) + 1) / 1000
    sg = arb(1) + arb(dd)
    margin = F['A'](sg, arb(t1r)) - F['S'](sg) - (arb(0.5) + arb(dd)) * V
    dmax[name] = {'delta_max': dd, 't1': '%.6f' % t1, 't1 rounded up': t1r,
                  'certified A-S-(1/2+d)V at t1 rounded up': margin.str(5), 'positive': bool(margin > 0),
                  'tau1 - t1 - (1/2+delta)/sqrt3': '%.4f' % (tau1 - t1r - (0.5 + float(dd)) / 3 ** 0.5),
                  'sigma_max = 1 + delta_max': '%.2f' % (1 + float(dd))}
rec('delta_max', dmax)


# far to the right, chi_4 and chi_12: v_sigma(t) - b_sigma = G_sigma(t) + P_sigma(t), with G_sigma the gamma part
# (1/2)[Re psi((sigma+1+it)/2) - psi((sigma+1)/2)] for chi_4 and (1/2)[Re psi((sigma+it)/2) - psi(sigma/2)] for chi_12, and
# P_sigma(t) >= -R_sigma(t), R_sigma(t) = sum over chi(n) = -1 of Lambda(n) n^{-sigma} (1 - cos(t log n)) <= the smaller of
# 2 S^-(sigma) and (t^2/2) M_2(sigma), where S^- = sum_{chi(n)=-1} Lambda(n) n^{-sigma} = (S + L'/L)/2 and
# M_2 = sum_{chi(n)=-1} Lambda(n) (log n)^2 n^{-sigma} = (S'' + (L'/L)'')/2, S = -(log lambda)' with
# lambda(s) = zeta(s) prod_{p | m} (1 - p^{-s}). G_sigma(t) increases in |t| and decreases in sigma, G_sigma(t)/t^2 decreases
# in |t|, and S^-, M_2 decrease in sigma. So on a cell [s0, s1] it suffices that G_{s1}(T) > max(2 S^-(s0), T^2 M_2(s0)/2)
# for one T > 0: then G > R at every t != 0 (|t| <= T by the second bound, |t| >= T by the first). Cells of width 1/32 on
# [1.5, 20], T in {0.25, 0.5, ..., 4}; beyond 20, T = 2, G_sigma(2) >= 1/((sigma+1)^2+4), and S^-, M_2 at sigma are at
# most 3^{-(sigma-20)} times their values at 20.
def gamma_part(name, sg, T):
    if name == 'chi4':
        return ((arb(sg) + 1 + acb(0, T)) / 2).digamma().real / 2 - ((arb(sg) + 1) / 2).digamma() / 2
    return ((arb(sg) + acb(0, T)) / 2).digamma().real / 2 - (arb(sg) / 2).digamma() / 2


def minus_part(name, sg):
    chi = CHI4 if name == 'chi4' else CHI12
    primes = [2] if name == 'chi4' else [2, 3]
    s = acb_series([acb(arb(sg)), 1], prec=4)
    cL = acb_series.dirichlet_l(s, chi).log().coeffs()
    lam = acb_series.zeta(s).log()
    for p in primes:
        lam += (1 - (-s * acb(arb(p).log())).exp()).log()
    cl = lam.coeffs()
    Sm = (-cl[1] + cL[1]).real / 2
    M2 = 3 * (cL[3] - cl[3]).real
    return Sm, M2


ctx.prec = 64
far = {}
Ts = [0.25 * j for j in range(1, 17)]
for name in ['chi4', 'chi12']:
    cells = []
    k = 0
    while True:
        s0 = 1.5 + k / 32
        if s0 >= 20:
            break
        s1 = s0 + 1 / 32
        Sm, M2 = minus_part(name, s0)
        ok = None
        for T in Ts:
            G = gamma_part(name, s1, T)
            if G > 2 * Sm and G > T * T * M2 / 2:
                ok = T
                break
        cells.append((s0, ok))
        k += 1
    j = len(cells)
    while j > 0 and cells[j - 1][1] is not None:
        j -= 1
    sigma2 = cells[j][0] if j < len(cells) else None
    Sm20, M220 = minus_part(name, 20)
    K = 2 * Sm20 if 2 * Sm20 > 2 * M220 else 2 * M220
    tail = K * (arb(21) ** 2 + 4)
    far[name] = {'cells': len(cells), 'all cells pass from sigma': sigma2,
                 'largest failing cell starts at': cells[j - 1][0] if j > 0 else None,
                 'T used at sigma2': cells[j][1] if j < len(cells) else None,
                 'S^-(4)': S_(minus_part(name, 4)[0], 6), 'M_2(4)': S_(minus_part(name, 4)[1], 6),
                 'G_4(1)': S_(gamma_part(name, 4, 1), 6),
                 'tail: K((20+1)^2+4) < 1': tail.str(5), 'tail ok': bool(tail < 1)}
rec('far_right', far)
ctx.prec = 128

# the ratio v_sigma(t)/(sigma - 1/2) along rows of sigma
rat = {}
for name, sgs, ts in [('chi4', ['0.55', '0.6', '0.75', '1', '1.5', '2.5', '3.15'], [0, 10, 40]),
                      ('xi', ['0.55', '0.6', '0.75', '1', '1.5', '2.5', '4.75'], [0, 20, 40]),
                      ('chi12', ['0.55', '0.6', '0.75', '1', '1.5', '2.3'], [0, 10, 40])]:
    F = FUN[name]
    rat[name] = {}
    for t in ts:
        vals = [F['LL'](acb(arb(sg), t)).real / (arb(sg) - 0.5) for sg in sgs]
        rat[name]['t=%d' % t] = {'sigma': sgs, 'ratio': [S_(v, 6) for v in vals],
                                 'decreasing': all(vals[k] > vals[k + 1] for k in range(len(vals) - 1))}
rec('ratio', rat)

# grid minima of v_sigma(t) - b_sigma, step 0.01, at 64 bits
ctx.prec = 64
gm = {}
for name, sgs, T in [('chi4', ['0.55', '0.6', '0.75', '1'], 100), ('chi12', ['0.55', '0.6', '0.75', '1'], 60),
                     ('chi4', ['3.5', '4', '5', '6', '8', '10', '12', '15', '20', '30'], 200),
                     ('chi12', ['2.5', '3', '4', '6', '10', '20'], 200)]:
    F = FUN[name]
    for sg in sgs:
        b = F['LL'](acb(arb(sg))).real
        mn, at = None, None
        mn_all = None
        for k in range(1, int(T * 100) + 1):
            t = k / 100
            d = float((F['LL'](acb(arb(sg), t)).real - b).mid())
            if mn_all is None or d < mn_all:
                mn_all = d
            if t >= 0.5 and (mn is None or d < mn):
                mn, at = d, t
        ctx.prec = 128
        v2 = -6 * F['logser'](acb(arb(sg)), 4).coeffs()[3].real
        ctx.prec = 64
        gm['%s sigma=%s' % (name, sg)] = {'T': T, 'b': S_(b, 8), 'min over [0.5,T] of v-b': '%.6g' % mn, 'at t': at,
                                          'min over grid of v-b': '%.3g' % mn_all, "v''(0)": S_(v2, 6)}
        print(name, sg, gm['%s sigma=%s' % (name, sg)], round(time.time() - t0, 1), flush=True)
rec('grid_minima', gm)
ctx.prec = 128

# where the comparison bound at sigma = 0.75 (delta = 2.15) meets b_{0.75}, for chi_4
F = FUN['chi4']
b75 = F['LL'](acb(0.75)).real
g = lambda t: arb(0.25) / arb(2.65) * (F['A'](arb(3.15), t) - F['S'](arb(3.15))) - b75
lo, hi = 0.0, 10.0
while hi - lo > 1e-9:
    m = (lo + hi) / 2
    if g(arb(m)) > 0:
        hi = m
    else:
        lo = m
rec('bound_meets_b075', {'b_0.75': S_(b75, 8), 't': '%.6f' % hi, '0.25 V': S_(0.25 * F['V'], 8)})
json.dump(out, open('dbeta_minimum.json', 'w'), indent=1)
print('done', round(time.time() - t0, 1))
