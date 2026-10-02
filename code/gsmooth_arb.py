# code/gsmooth_arb.py (October 2026): certificates, in ball arithmetic (python-flint, Arb), for the Gaussian smoothing
# statements of Section ch:finite. For each case of firstfail.json, the closest pair A < B of consecutive zeros in a
# window moved to the quadruple 1/2 +- delta +- i g, g = (A+B)/2, and with one zero of each pair {rho, 1-rho} of F and
# rho = 1/2 + i gamma~,
#   G(u, t) = sum exp(-2t (gamma~ - u)^2),
# in which the line zeros give positive terms and the quadruple gives 2 exp(2t(delta^2 - (g-u)^2)) cos(4 t delta (g-u)).
# G_0 denotes the same sum over the zeros of F below H = 3e12, all on the line except the quadruple [PT21].
# t1 is the failure time gauss_t found by firstfail.py rounded down to two decimals (1.71 for the pair near 1977.22).
#  (1) Cover: for u in [g-8, g+8], G_0(u, t1) >= m0 > 0, from the line zeros within 40 of g (the others
#      add positive terms) and the quadruple, on intervals [c-r, c+r] by f(c) - |f'(c)| r - sup|f''| r^2/2.
#  (2) Domination: for u >= g+8 the first line zero gamma_+ above g outweighs the modulus of the quadruple term, since
#      (gamma_+ - g)(2u - g - gamma_+) >= delta^2 + log 2/(2 t1) at u = g+8 and the left side increases with u;
#      likewise gamma_- below g for u <= g-8. So G_0(., t1) >= 0 on the whole line.
#  (3) Convolution: sqrt(2t/pi) G_0(., t) is sqrt(2 t1/pi) G_0(., t1) convolved with the centred normal density of
#      variance s^2 = 1/(4t) - 1/(4 t1), term by term. For u within 8 of g and 0 < t <= t1 this gives
#      G_0(u, t) >= c m0, c = min(Phi(2) - 1/2, 16 phi(2) sqrt(t1)): if s <= 4 the normal law puts at least
#      Phi(2) - 1/2 on the half of the window beyond u, and if s > 4 at least 16 sqrt t phi(2), while sqrt(t1/t) >= 1.
#  (4) Zeros above H, wherever they lie: a quadruple 1/2 +- b +- i gamma (|b| < 1/2) gives
#      2 exp(-2t((gamma-u)^2 - b^2)) cos(4 t b (gamma - u)), which is negative only if gamma - u > pi/(4t). So for
#      u <= g+8 their sum is at least -e^{t/2} sum exp(-2t(gamma-u)^2) over gamma > u + X, X = max(H - u, pi/(4t)),
#      and partial summation against N(s) <= s^2 (s >= 100, from Trudgian's bound) gives at most
#      e^{t/2} (1 + u/X)^2 (1 + 2tX^2) e^{-2tX^2}/(2t) <= e^{t1/2} (1 + u/X)^2 (4/pi^2) w0 (1 + w0) e^{-w0},
#      w0 = pi (H - u)/2, for every t <= t1. Hence G(u, t) >= c m0 - far > 0 for |u - g| <= 8 and 0 < t <= t1. This
#      bounds only how far these zeros can lower G: their positive part is unbounded as t -> 0. The record gives c m0
#      (lower_bound_zeros_below_H) and c m0 - far (lower_bound_for_all_t_le_t1), both rounded down.
#  (5) A certified negative value G(u*, t2) < 0 at the failure time rounded up to four decimals (or a little later),
#      with every zero accounted for: the window, at most
#      N(g-40) line zeros below it, the line zeros between g+40 and H (partial summation as in (4)) and those above H.
# The zeros are computed by Arb (acb.zeta_zeros: isolated, counted by Turing's method, refined in ball arithmetic).
# Reads firstfail.json (the pairs and the floating point failure times); writes gsmooth_arb.json. Seven cases in parallel,
# a few minutes.
import json
import math
from decimal import Decimal
from fractions import Fraction
from multiprocessing import Pool

from flint import acb, arb, ctx


def dec(x, n, up):
    """An n significant digit decimal string at most every point of the ball x (up=False) or at least every point
    (up=True): the exact binary endpoint, rounded outward in exact rational arithmetic."""
    m, k = (int(v) for v in (x.upper() if up else x.lower()).man_exp())
    q = Fraction(m)*Fraction(2)**k
    if q == 0:
        return '0'
    d = len(str(abs(q.numerator))) - len(str(q.denominator))
    while Fraction(10)**d > abs(q):
        d -= 1
    while Fraction(10)**(d + 1) <= abs(q):
        d += 1
    j = q/Fraction(10)**(d - n + 1)
    return str(Decimal(math.ceil(j) if up else math.floor(j)).scaleb(d - n + 1))


ctx.prec = 192
PI = arb.pi()
H = arb(3)*arb(10)**12
W = 40


def N_upper(T):   # Trudgian: N(T) <= (T/2pi) log(T/(2 pi e)) + 7/8 + 0.112 log T + 0.278 log log T + 2.510 + 0.2/T, T >= e
    return (T/(2*PI)*(T/(2*PI*arb(1).exp())).log() + arb(7)/8 + arb('0.112')*T.log()
            + arb('0.278')*T.log().log() + arb('2.510') + arb('0.2')/T)


def fl(x):
    return float(x.mid().str(20, radius=False))


def run_case(case):
    ctx.prec = 192
    num, den = case['delta'].split('/')
    DELTA = arb(num)/arb(den)
    t1s = '%.2f' % (math.floor(case['gauss_t']*100)/100)
    T1 = arb(t1s)
    nA = case['zero_numbers'][0]
    first = max(1, nA - 70)
    Z = acb.zeta_zeros(first, 141)
    half = arb(1)/2
    assert all(z.real == half for z in Z)
    GAM = {first + i: z.imag for i, z in enumerate(Z)}
    A, B = GAM[nA], GAM[nA + 1]
    assert abs(A - arb(case['A'])) < arb('1e-40') and abs(B - arb(case['B'])) < arb('1e-40')
    g = (A + B)/2
    gf = fl(g)
    near = [n for n in GAM if n not in (nA, nA + 1) and abs(fl(GAM[n]) - gf) <= W]
    n0, n1 = min(near), max(near)
    assert GAM[n0 - 1] < g - W and g - W < GAM[n0] and GAM[n1] < g + W and g + W < GAM[n1 + 1]
    assert len(near) == n1 - n0 - 1
    LINE = [GAM[n] for n in sorted(near)]
    zq = acb(g, -DELTA)

    def f012(u, t):
        """G over the window zeros and the quadruple at the real ball u: value, first and second derivative in u."""
        f0 = arb(0); f1 = arb(0); f2 = arb(0)
        for gam in LINE:
            d = gam - u
            e = (-2*t*d*d).exp()
            f0 += e; f1 += 4*t*d*e; f2 += (16*t*t*d*d - 4*t)*e
        d = zq - acb(u)
        E = (-2*t*d*d).exp()
        f0 += 2*E.real; f1 += 2*(4*t*d*E).real; f2 += 2*((16*t*t*d*d - 4*t)*E).real
        return f0, f1, f2

    def cover(t, lo, hi, w0=Fraction(1, 64), wmin=Fraction(1, 2**30)):
        cells = []
        a = lo
        while a < hi:
            cells.append((a, a + w0)); a += w0
        stack = cells[::-1]
        n = 0; m0 = None; at = None; fmin = None; ufmin = None
        while stack:
            a, b = stack.pop()
            c = (a + b)/2; r = (b - a)/2
            cb = arb(c.numerator)/c.denominator; rb = arb(r.numerator)/r.denominator
            F0, F1, _ = f012(cb, t)
            _, _, F2 = f012(cb + rb*arb(0, 1), t)
            lb = F0 - abs(F1)*rb - abs(F2)*rb*rb/2
            v = fl(F0)
            if fmin is None or v < fmin:
                fmin, ufmin = v, float(c)
            if lb > 0:
                n += 1
                x = lb.lower()
                if m0 is None or x < m0:
                    m0, at = x, (float(a), float(b))
                continue
            if F0 < 0:
                return dict(ok=False, negative_at=float(c), value=F0.str(10), cells=n)
            if b - a <= wmin:
                return dict(ok=False, undecided=(float(a), float(b)), lb=lb.str(6), cells=n)
            m = (a + b)/2
            stack += [(m, b), (a, m)]
        return dict(ok=True, cells=n, m0=m0, at=at, fmin=fmin, u_fmin=ufmin)

    out = dict(window=case['window'], delta=case['delta'], pair=[nA, nA + 1], A=A.str(25), B=B.str(25), g=g.str(25),
               t1=t1s, float_failure_time=case['gauss_t'], window_zero_numbers=[n0 + 1, n1 - 1], n_window=len(LINE))
    # (1) the cover of [g-8, g+8], extended to dyadic ends
    lo = Fraction(math.floor((gf - 8)*64) - 1, 64); hi = Fraction(math.ceil((gf + 8)*64) + 1, 64)
    assert arb(lo.numerator)/lo.denominator < g - 8 and arb(hi.numerator)/hi.denominator > g + 8
    cv = cover(T1, lo, hi)
    m0 = cv.pop('m0', None)
    if cv['ok']:
        cv['m0'] = m0.str(8, radius=False)
        cv['u_fmin_minus_g'] = cv['u_fmin'] - gf
    out['cover'] = cv
    # (2) domination outside the window
    up = min(n for n in GAM if GAM[n] > g and n not in (nA, nA + 1))
    dn = max(n for n in GAM if GAM[n] < g and n not in (nA, nA + 1))
    need = DELTA*DELTA + arb(2).log()/(2*T1)
    dp = GAM[up] - g; dm = g - GAM[dn]
    dom_up = dp*(16 - dp); dom_dn = dm*(16 - dm)
    out['domination'] = dict(zero_above=up, zero_below=dn, gap_above=dp.str(10), gap_below=dm.str(10), need=need.str(10),
                             above=dom_up.str(10), below=dom_dn.str(10),
                             holds=bool(dp < 16 and dm < 16 and dom_up > need and dom_dn > need))
    # (3) the constant of the convolution step
    phi2 = (-arb(2)).exp()/(2*PI).sqrt()
    Phi2 = (1 + (arb(2)/arb(2).sqrt()).erf())/2
    cconst = min(Phi2 - half, 16*phi2*T1.sqrt())
    # (4) zeros above H, for every t <= t1 and u <= g + 8
    umax = g + 8
    w0 = PI*(H - umax)/2
    X = H - umax
    far = (T1/2).exp()*(1 + umax/X)**2*4/PI**2*w0*(1 + w0)*(-w0).exp()
    out['convolution_constant'] = cconst.str(10)
    out['far_bound_log10'] = (far.log()/arb(10).log()).str(8)
    ok = bool(cv['ok'] and out['domination']['holds'] and cconst*m0 > far)
    out['positive_for_all_t_le_t1'] = ok
    # the bound for the zeros below H, and for the whole sum after the zeros above H have lowered it by at most far
    out['lower_bound_zeros_below_H'] = dec(cconst*m0, 6, False) if cv['ok'] else None
    out['lower_bound_for_all_t_le_t1'] = dec(cconst*m0 - far, 6, False) if cv['ok'] else None

    # (5) a negative value just above t1
    def upper_G(u, t):
        f0, _, _ = f012(u, t)
        below = N_upper(g - W)*(-2*t*(u - (g - W))**2).exp()                # u - gamma >= u - (g - W) > 0
        Xa = g + W - u
        above = (1 + u/Xa)**2*(1 + 2*t*Xa*Xa)*(-2*t*Xa*Xa).exp()/(2*t)        # zeros in (g+W, inf), partial summation
        Xh = H - u
        farh = (t/2).exp()*(1 + u/Xh)**2*(1 + 2*t*Xh*Xh)*(-2*t*Xh*Xh).exp()/(2*t)
        return f0 + below + above + farh, f0, below, above

    tg = case['gauss_t']; ug = case['gauss_u_minus_g']
    base = math.ceil(tg*1e4)/1e4
    for t2s in ['%.4f' % base, '%.4f' % (base + 1e-4), '%.4f' % (base + 5e-4), '%.4f' % (base + 2e-3), '%.4f' % (base + 1e-2)]:
        t2 = arb(t2s)
        i0 = int(round(ug*1000))
        best = None
        for i in range(i0 - 300, i0 + 301):                                    # u - g within 0.3 of the float minimiser
            v = fl(f012(g + arb(i)/1000, t2)[0])
            if best is None or v < best[0]:
                best = (v, i, 0)
        for j in range(-100, 101):                                              # refine to step 1e-5
            v = fl(f012(g + arb(best[1])/1000 + arb(j)/100000, t2)[0])
            if v < best[0]:
                best = (v, best[1], j)
        du = arb(best[1])/1000 + arb(best[2])/100000
        U, f0, below, above = upper_G(g + du, t2)
        out['negative'] = dict(t=t2s, u_minus_g=du.str(8, radius=False), upper=U.str(8),
                               below=below.str(3), above=above.str(3), certified=bool(U < 0))
        if U < 0:
            break
    out['all_certified'] = bool(ok and out['negative']['certified'])
    print(json.dumps(out), flush=True)
    return out


if __name__ == '__main__':
    cases = json.load(open('firstfail.json'))['cases']
    with Pool(len(cases)) as pool:
        res = pool.map(run_case, cases)
    json.dump(dict(prec=192, all_certified=all(r['all_certified'] for r in res), cases=res), open('gsmooth_arb.json', 'w'),
              indent=1)
    for r in res:
        print('%s delta %s: G > 0 for |u-g| <= 8 and every t <= %s (G >= %s); G(g%+.5f, %s) <= %s; far zeros < 10^%s; %s'
              % (r['window'], r['delta'], r['t1'], r['lower_bound_for_all_t_le_t1'], float(r['negative']['u_minus_g']),
                 r['negative']['t'], r['negative']['upper'], r['far_bound_log10'].split(' ')[0].lstrip('['), r['all_certified']))
    print('all certified:', all(r['all_certified'] for r in res))
