# maxexp_psizeros.py: the zeros of the Laplace transform psi*(theta) = E exp(-theta M*) of the Delta max law, for Section
# sec:ch16:max (ch/ch16.tex, Section 8.5.4) and Figure fig:maxexp:laws (right). Added for the book in October 2026; the first
# author's note on maxima of exponentials (misc/reviews/nick-maxexp/) computed psi* on |Re theta| < 15 pi only, where its one
# zero is real, and maxlaw.py confirms that. Here the region is |Re theta| <= 2 pi (MMAX + 1/2), |Im theta| <= 100.
# psi* is meromorphic with simple poles at -2 pi m, residue -2 pi m tau(m+1), and for every c > 0 (the cut form)
#   psi*(theta) = theta int_0^c e^{-theta y} F(y) dy + theta e^{-theta c} sum_{n>=1} tau(n) e^{-2 pi (n-1) c}/(theta + 2 pi (n-1)),
# F(y) = prod (1 - e^{-2 pi n y})^24 = y^-12 e^{2 pi y - 2 pi/y} F(1/y).
# 1. Scanning evaluator: c = 0.08, the integral by Gauss Legendre on [0.012, 0.08] (6 panels of 96 nodes; F < 1e-200 below
#    0.012), the series to n = 900, python-flint balls at 230 bits (midpoints used); checked against c = 0.05.
# 2. Real axis: the sign of psi* at 160 points in each gap (-2 pi (m+1), -2 pi m), m <= MMAX, against the residue signs, and
#    each sign change refined by bisection.
# 3. Argument principle: the change of arg psi* along the boundary of rectangles |Re| <= 2 pi (k + 1/2), |Im| <= H, and of the
#    cells (2 pi wide, 5 high) of a grid of the upper half plane, Im in [0.5, 30], each step halved until arg changes by < 0.3.
#    Each cell with winding number one is quartered down to side 1e-4, keeping the quarter of winding number one, and Newton's
#    method refines the zero.
# 4. Certificates in ball arithmetic, with c = 3/100 and the series to n = 1500, the tail bounded by |tau(n)| <= d(n) n^{11/2}
#    <= 2 n^6 (Deligne) and the integral bounded by F(y) <= y^-12 e^{2 pi y - 2 pi/y} on (0, c]: each real zero by opposite
#    signs at two points 1e-20 apart; each nonreal zero z0 by Rouche: with A >= |psi*(z0)|, B <= |psi*'(z0)| and C >= |psi*''|
#    on the box of half side r around z0, A + 2 C r^2 < B r/2 (the box contains the disc of radius r, on whose boundary
#    |psi* - psi*'(z0)(theta - z0)| <= A + C r^2/2 < B r <= |psi*'(z0)(theta - z0)|), so psi* has exactly one zero in that disc;
#    this is tested for r = 10^-6, ..., 10^-29: the radius 10^-6 gives uniqueness, the smallest radius the accuracy.
# Run from any folder: python3 maxexp_psizeros.py; writes maxexp_psizeros.json there.
import json, time
import mpmath as mp
from flint import acb, arb, ctx

t0 = time.time()
ctx.prec = 230
mp.mp.dps = 75
out = {}
def put(k, v):
    out[k] = v
    print(k, v, round(time.time() - t0, 1), flush=True)
MMAX = 60
NT = 1600
sig = [0] * (NT + 1)
for d in range(1, NT + 1):
    for m in range(d, NT + 1, d):
        sig[m] += d
tau = [0] * (NT + 1); tau[1] = 1
for n in range(2, NT + 1):
    tau[n] = -24 * sum(sig[m] * tau[n - m] for m in range(1, n)) // (n - 1)
assert tau[2] == -24 and tau[3] == 252 and tau[11] == 534612 and tau[13] == -577738
PI = arb.pi(); TWO_PI = 2 * PI
PIF = 3.141592653589793
def fl(z):
    return complex(float(z.real.mid()), float(z.imag.mid()))

def F_small(y):
    q = (-TWO_PI / y).exp()
    p = arb(1); qn = q
    for n in range(1, 80):
        p *= (1 - qn) ** 24
        qn *= q
        if qn < arb(10) ** -90:
            break
    return y ** -12 * (TWO_PI * y - TWO_PI / y).exp() * p
gl = mp.calculus.quadrature.GaussLegendre(mp.mp)
BASE = gl.calc_nodes(6, mp.mp.prec)

class Psi:
    # the cut form with Gauss Legendre for the integral (not certified; used for scanning and Newton)
    def __init__(self, c, y0, nterms=900, panels=6):
        self.c = arb(c); self.N = nterms
        ys, ws = [], []
        h = (mp.mpf(c) - mp.mpf(y0)) / panels
        for k in range(panels):
            lo = mp.mpf(y0) + k * h
            for x, w in BASE:
                ys.append(lo + h / 2 * (x + 1)); ws.append(h / 2 * w)
        self.y = [arb(mp.nstr(y, 70)) for y in ys]
        self.wF = [arb(mp.nstr(w, 70)) * F_small(arb(mp.nstr(y, 70))) for y, w in zip(ys, ws)]
        self.a = [arb(tau[n]) * (-TWO_PI * (n - 1) * self.c).exp() for n in range(1, nterms + 1)]
    def __call__(self, th, der=False):
        th = acb(th)
        I = acb(0); I1 = acb(0)
        for y, w in zip(self.y, self.wF):
            e = w * (-th * y).exp()
            I += e
            if der: I1 -= y * e
        S0 = acb(0); S1 = acb(0)
        for n in range(self.N):
            r = 1 / (th + TWO_PI * n)
            S0 += self.a[n] * r
            if der: S1 -= self.a[n] * r * r
        E = (-th * self.c).exp()
        f = th * (I + E * S0)
        if not der:
            return f
        return f, I + th * I1 + E * (1 - self.c * th) * S0 + th * E * S1
P1 = Psi('0.08', '0.012')
P2 = Psi('0.05', '0.008', nterms=1400)
chk = {}
for th in [acb(2), acb(-100, 3), acb(-300, 1), acb(-50, 80), acb(20, 100)]:
    a = P1(th); b = P2(th)
    chk[str(fl(th))] = {'psi*': str(fl(a)), 'relative difference of the cuts 0.08 and 0.05': float(abs(a - b).mid() / abs(a).mid())}
put('cut_consistency', chk)
put('psi*(2)', P1(acb(2)).real.str(25, radius=False))

# 2. the real axis
gaps = []; realzeros = []
for m in range(0, MMAX + 1):
    # m = 0 is the gap (-2 pi, 0): psi*(0) = 1 and psi* > 0 just right of -2 pi (residue -2 pi tau(2) > 0), so nothing is forced
    lo = -2 * mp.pi * (m + 1); hi = -2 * mp.pi * m
    xs = [lo + (hi - lo) * (j + mp.mpf(1) / 2) / 160 for j in range(160)]
    vals = [float(P1(acb(arb(mp.nstr(x, 50)))).real.mid()) for x in xs]
    ch = [j for j in range(159) if vals[j] * vals[j + 1] < 0]
    gaps.append({'m': m, 'tau(m+1) tau(m+2) > 0': m >= 1 and tau[m + 1] * tau[m + 2] > 0, 'sign changes': len(ch)})
    for j in ch:
        a, b = xs[j], xs[j + 1]; fa = vals[j]
        for _ in range(110):
            mid = (a + b) / 2
            fm = float(P1(acb(arb(mp.nstr(mid, 50)))).real.mid())
            if fm * fa < 0: b = mid
            else: a, fa = mid, fm
        realzeros.append((a + b) / 2)
put('gaps_summary', {'gaps': MMAX + 1, 'gaps with tau(m+1) tau(m+2) > 0': sum(1 for g in gaps if g['tau(m+1) tau(m+2) > 0']),
                     'real zeros': len(realzeros),
                     'every same sign gap has exactly one sign change, every other gap none': all((g['sign changes'] == 1) == g['tau(m+1) tau(m+2) > 0'] and g['sign changes'] <= 1 for g in gaps)})
out['gaps'] = gaps
put('real_zeros', [mp.nstr(z, 20) for z in realzeros])

# 3. argument principle
def darg(f, z0, z1, f0, f1, depth=0):
    dv = float((f1 / f0).arg().mid())
    if abs(dv) < 0.3 or depth > 25:
        return dv, depth > 25
    zm = (z0 + z1) / 2; fm = f(zm)
    d0, e0 = darg(f, z0, zm, f0, fm, depth + 1); d1, e1 = darg(f, zm, z1, fm, f1, depth + 1)
    return d0 + d1, e0 or e1
CACHE = {}
WORST = [0.0]   # the largest relative radius of the balls met on the contours (a check that 230 bits suffice)
def fval(z):
    key = (round(float(z.real.mid()), 9), round(float(z.imag.mid()), 9))
    if key not in CACHE:
        v = P1(z)
        WORST[0] = max(WORST[0], float(abs(v).rad()) / float(abs(v).mid()))
        CACHE[key] = v
    return CACHE[key]
def seg_arg(p, q, h0=0.25):
    n = max(1, int(abs(fl(q) - fl(p)) / h0 + 0.999))
    zs = [p + (q - p) * arb(j) / n for j in range(n + 1)]; fs = [fval(z) for z in zs]
    tot = 0.0; bad = False
    for j in range(n):
        d, e = darg(P1, zs[j], zs[j + 1], fs[j], fs[j + 1]); tot += d; bad = bad or e
    return tot, bad
def winding(corners):
    tot = 0.0; bad = False
    for p, q in zip(corners, corners[1:] + corners[:1]):
        d, e = seg_arg(p, q); tot += d; bad = bad or e
    return tot / (2 * PIF), bad
rect = []
for k, Hh in [(7, 60), (15, 30), (25, 30), (40, 30), (MMAX, 30), (MMAX, 100)]:
    R = TWO_PI * (k + arb(1) / 2)
    w, bad = winding([acb(-R, -Hh), acb(R, -Hh), acb(R, Hh), acb(-R, Hh)])
    nreal = sum(1 for z in realzeros if z > -2 * mp.pi * (k + mp.mpf(1) / 2))
    rect.append({'k': k, 'H': Hh, 'winding': round(w, 6), 'zeros (winding + k poles)': int(round(w)) + k, 'real zeros inside': nreal, 'depth limit': bad})
    print(rect[-1], round(time.time() - t0, 1), flush=True)
put('rectangles', rect)
KC = MMAX
xs = [-TWO_PI * (KC + arb(1) / 2) + TWO_PI * j for j in range(KC + 7)]   # up to 2 pi (5 + 1/2)
ys = [arb('0.5')] + [arb(5 * j) for j in range(1, 7)]
hor = {}; ver = {}
for i in range(len(xs) - 1):
    for j in range(len(ys)):
        hor[(i, j)] = seg_arg(acb(xs[i], ys[j]), acb(xs[i + 1], ys[j]))[0]
for i in range(len(xs)):
    for j in range(len(ys) - 1):
        ver[(i, j)] = seg_arg(acb(xs[i], ys[j]), acb(xs[i], ys[j + 1]))[0]
cells = []
for i in range(len(xs) - 1):
    for j in range(len(ys) - 1):
        w = (hor[(i, j)] + ver[(i + 1, j)] - hor[(i, j + 1)] - ver[(i, j)]) / (2 * PIF)
        if abs(w) > 0.5:
            cells.append((i, j, round(w, 4)))
put('cells_with_zeros', [{'Re/2pi': [float(xs[i].mid()) / (2 * PIF), float(xs[i + 1].mid()) / (2 * PIF)], 'Im': [float(ys[j].mid()), float(ys[j + 1].mid())], 'winding': w} for i, j, w in cells])
print('cells done', round(time.time() - t0, 1), flush=True)
# each cell is quartered, keeping the quarter of winding number one, down to side 1e-4; Newton's method then refines the zero
def wind_box(x0, x1, y0, y1):
    corners = [acb(x0, y0), acb(x1, y0), acb(x1, y1), acb(x0, y1)]
    tot = 0.0
    for a in range(4):
        p, q = corners[a], corners[(a + 1) % 4]
        n = max(4, int(abs(fl(q) - fl(p)) / 0.25 + 0.999))
        zs = [p + (q - p) * arb(j) / n for j in range(n + 1)]; fs = [P1(z) for z in zs]
        for j in range(n):
            tot += darg(P1, zs[j], zs[j + 1], fs[j], fs[j + 1])[0]
    return tot / (2 * PIF)
def locate(x0, x1, y0, y1):
    while float((x1 - x0).mid()) > 1e-4 or float((y1 - y0).mid()) > 1e-4:
        found = False
        for shift in (arb(1) / 2, arb(37) / 73, arb(29) / 61):   # a zero on a dividing line moves the line
            xm = x0 + (x1 - x0) * shift; ym = y0 + (y1 - y0) * shift
            for bx in ((x0, xm), (xm, x1)):
                for by in ((y0, ym), (ym, y1)):
                    if abs(wind_box(bx[0], bx[1], by[0], by[1]) - 1) < 0.01:
                        x0, x1 = bx; y0, y1 = by; found = True; break
                if found: break
            if found: break
        if not found:
            return None
    return x0, x1, y0, y1
def newton(z):
    # midpoints only: a ball iterate would inflate its radius at every step through the cancellation in the cut form
    for _ in range(60):
        f, d = P1(z, der=True)
        dz = f / d
        z = acb(z.real.mid(), z.imag.mid()) - acb(dz.real.mid(), dz.imag.mid())
        if float(abs(dz).mid()) < 1e-55: break
    return z
nonreal = []
for i, j, w in cells:
    box = locate(xs[i], xs[i + 1], ys[j], ys[j + 1])
    if box is None:
        nonreal.append((None, False)); continue
    z = newton(acb((box[0] + box[1]) / 2, (box[2] + box[3]) / 2))
    ok = float(box[0].mid()) - 1e-4 < fl(z).real < float(box[1].mid()) + 1e-4 and float(box[2].mid()) - 1e-4 < fl(z).imag < float(box[3].mid()) + 1e-4
    nonreal.append((z, ok))
    print('zero', z.real.str(20, radius=False), z.imag.str(20, radius=False), ok, round(time.time() - t0, 1), flush=True)
nonreal = [(z, ins) for z, ins in nonreal if z is not None]
put('nonreal_zeros_upper_half', [{'theta': '%s + %si' % (z.real.str(22, radius=False), z.imag.str(22, radius=False)), 're': float(z.real.mid()), 'im': float(z.imag.mid()),
                                  'theta/2pi': str(fl(z) / (2 * PIF)), 'Newton stayed in the box': ins, '|psi*|': float(abs(P1(z)).mid())} for z, ins in nonreal])

# 4. certificates in ball arithmetic
CR = arb(3) / 100; NR = 1500
AR = [arb(tau[n]) * (-TWO_PI * (n - 1) * CR).exp() for n in range(1, NR + 1)]
def upper(x):   # an arb ball [0, u] with u >= every value of |x|
    return arb(abs(x).upper())
def psi_rig(th):
    # returns balls for psi*, psi*', psi*'' on the acb box th, with the integral and the tail of the series in the radius
    E = (-th * CR).exp()
    S = [acb(0), acb(0), acb(0)]
    for n in range(NR):
        r = 1 / (th + TWO_PI * n)
        t = AR[n] * r
        S[0] += t; t *= r; S[1] -= t; t *= r; S[2] += 2 * t
    g0 = th * E; g1 = E * (1 - CR * th); g2 = -CR * E * (2 - CR * th)
    T = [g0 * S[0], g1 * S[0] + g0 * S[1], g2 * S[0] + 2 * g1 * S[1] + g0 * S[2]]
    am = upper(th); rm = upper(th.real)
    Fmax = CR ** -12 * (TWO_PI * CR - TWO_PI / CR).exp()           # F <= y^-12 e^{2 pi y - 2 pi/y} <= this on (0, c]
    tail = arb(0)                                                   # sum_{n > NR} 2 n^6 e^{-2 pi (n-1) c}/(2 pi (n-1) - |theta| - 1)
    nn = NR + 1
    t1 = 2 * arb(nn) ** 6 * (-TWO_PI * (nn - 1) * CR).exp() / (TWO_PI * (nn - 1) - am - 1)
    ratio = (arb(nn + 1) / nn) ** 6 * (-TWO_PI * CR).exp()
    assert ratio < 1 and (TWO_PI * (nn - 1) - am - 1) > 1
    tail = t1 / (1 - ratio)
    res = []
    for k in range(3):
        Jk = (am * CR ** (k + 1) + k * CR ** k) * (rm * CR).exp() * Fmax
        Tk = arb(2 if k == 2 else 1) * (am + 1) * ((rm + 1) * CR).exp() * tail
        err = Jk + Tk
        res.append(T[k] + acb(arb(0, 1) * err, arb(0, 1) * err))
    return res
cert_real = []
for z in realzeros:
    lo = acb(arb(mp.nstr(z - mp.mpf('1e-20'), 60))); hi = acb(arb(mp.nstr(z + mp.mpf('1e-20'), 60)))
    a = psi_rig(lo)[0].real; b = psi_rig(hi)[0].real
    cert_real.append(bool((a > 0 and b < 0) or (a < 0 and b > 0)))
put('real_zeros_certified', {'all': all(cert_real), 'count': sum(cert_real)})
cert = []
for z, ins in nonreal:
    if not ins:
        cert.append({'theta': str(fl(z)), 'certified': False}); continue
    zc = acb(arb(z.real.mid()), arb(z.imag.mid()))
    v0, v1, _ = psi_rig(zc)
    A = arb(abs(v0).upper()); B = arb(abs(v1).lower())
    good = []
    for e in range(6, 30):
        r = arb(10) ** -e
        box = acb(arb(z.real.mid()) + arb(0, 1) * r, arb(z.imag.mid()) + arb(0, 1) * r)
        C = arb(abs(psi_rig(box)[2]).upper())
        if A + C * r * r * 2 < B * r / 2:   # the box has half side r and contains the disc of radius r
            good.append(e)
    # Rouche at radius 10^-6 gives exactly one zero in that disc, and at the smallest radius 10^-e it lies within 10^-e of z0
    ok = 6 in good
    cert.append({'theta': '%s + %si' % (zc.real.str(25, radius=False), zc.imag.str(25, radius=False)), '|psi*(z0)| <=': A.str(3), '|psi*\'(z0)| >=': B.str(5),
                 'exactly one zero in the disc of radius 1e-6': ok, 'and it lies within 10^-e of z0, e': max(good) if good else None,
                 'radii 10^-e certified, e': good, 'certified': ok})
put('nonreal_zeros_certified', cert)
# the counts: each rectangle should hold its real zeros and the nonreal zeros found in the cells, with their conjugates
cons = []
for d in rect:
    Rr = 2 * PIF * (d['k'] + 0.5)
    nn_ = sum(2 for z, ins in nonreal if abs(fl(z).real) < Rr and fl(z).imag < d['H'])
    cons.append({'k': d['k'], 'H': d['H'], 'zeros by the argument principle': d['zeros (winding + k poles)'], 'real + nonreal found': d['real zeros inside'] + nn_})
put('counts_consistent', {'all': all(c['zeros by the argument principle'] == c['real + nonreal found'] for c in cons) and all(abs(w - 1) < 1e-3 for i, j, w in cells), 'rectangles': cons})
put('largest relative radius on the contours', WORST[0])
out['poles_residue_sign'] = [{'m': m, 'residue sign': -1 if tau[m + 1] > 0 else 1} for m in range(1, MMAX + 1)]
out['seconds'] = time.time() - t0
json.dump(out, open('maxexp_psizeros.json', 'w'), indent=1)
print('done', round(time.time() - t0, 1))
