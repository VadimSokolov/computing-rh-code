# Chapter 13, Section ch:fminterval: certificates, in ball arithmetic (python-flint, Arb), for the zero counts and zeros
# that fm_interval.py computes in floating point.
# Counting. For f analytic on a polygon, the number of zeros inside is (1/2pi) times the change of arg f around it. Each
# edge is cut into segments S; when |f(z) - f(m)| < |f(m)| for every z in a ball B containing S, m the midpoint, f(S)
# lies in the disk of centre f(m) and radius |f(m)|, so arg f changes along S by the principal value of
# arg(f(end)/f(start)). Arb encloses f on B, the segments that fail the test are bisected, and the sum of the principal
# values, divided by 2pi, is a ball that must contain exactly one integer. A zero on the polygon makes the test fail on
# ever shorter segments, and the program stops.
# (1) L(s) = sinh(s)/s + c, c = 1/10. It is even and real on the real axis, so its zeros are symmetric about both axes.
#     On [-6, 6] x [0, 65] (no zero on the real axis, where L >= 1 + c) the count is 20. For |Re s| >= 6 and
#     0 <= Im s <= 65 there is no zero, since |sinh s| >= sinh|Re s| > c(|Re s| + 65) >= c|s|. Each of the 20 zeros
#     is isolated in a box of half width 1e-10 with count 1: a box symmetric about the imaginary axis for the two zeros
#     on it (count 1 there puts the zero on the axis, its mirror image being a zero too), and boxes about x_j + i y_j,
#     j = 1..9, and their mirror images. So the list is complete, every zero is simple, and the printed digits are
#     checked against the boxes.
# (2) L_3(s) = 6 (sinh s - s)/s^3 and L_4(s) = 24 (cosh s - 1 - s^2/2)/s^4. Both are even and real, positive on the
#     real and the imaginary axis, and |L_k(s) - 1| < 1 for |s| <= 1 (the Taylor coefficients are k!/(2j + k)!). The
#     count is 1 on the square [0, R]^2 with the corner x + y < 1 cut off, R = |z_1| + 1/100, and the zero z_1 is
#     isolated in a box of half width 1e-10; so z_1 is the only zero of modulus at most R in the closed first quadrant.
# (3) F_b(z) = int_0^1 (1 - x^b) cos(zx) dx = sin(z)/z - (M(b+1, b+2, iz) + M(b+1, b+2, -iz))/(2(b+1)), M the
#     confluent hypergeometric function, for the eleven a = b + 1 of Table fm:tab:power. F_b is real on the real axis.
#     The count in the box (0.3, 60) x (-12, 12) is 18 for every a. For b > 1, certified sign changes of F_b on the grid
#     of step 0.002 give at least 18 real zeros, so all 18 zeros in the box are real and simple. For b < 1 the count
#     in (0.3, 60) x (1/20, 12) is 9, hence 9 in the mirror image, so none of the 18 lies within 1/20 of the real axis.
#     F_b(2 pi) is enclosed for each a.
# (4) The constants of Section ch:fmexp quoted with fm_interval.py: the first positive root of tan y = y, where
#     sin(y)/y is smallest, and that minimum (bracketed by signs of sin y - y cos y), arccosh(3/2) and sqrt(10) - 2.
# Usage: python3 fm_intcert.py [workers]      writes fm_intcert.json
import json, math, sys, time
from multiprocessing import Pool
from flint import arb, acb, ctx

ctx.prec = 128
I = acb(0, 1)
T0 = time.time()


def ball(z1, z2):
    """A complex ball containing the segment [z1, z2]."""
    m = (z1 + z2)/2
    u = arb(0, 1)                                          # the ball [-1, 1]
    return m + acb(u*abs((z2 - z1).real)/2, u*abs((z2 - z1).imag)/2)


def winding(f, corners, n0=64, minw=2.0**-40):
    """Certified number of zeros of f inside the polygon (counterclockwise corners, given as exact floats)."""
    P = [acb(arb(x), arb(y)) for x, y in corners]
    total = arb(0); segs = 0
    for e in range(len(P)):
        p, q = P[e], P[(e + 1) % len(P)]
        pt = lambda u: p + (q - p)*arb(u)
        val = {}

        def F(u):
            if u not in val:
                val[u] = f(pt(u))
            return val[u]
        stack = [(k/n0, (k + 1)/n0) for k in range(n0 - 1, -1, -1)]
        while stack:
            u, v = stack.pop()
            fm = f(pt((u + v)/2))
            fb = f(ball(pt(u), pt(v)))
            if abs(fb - fm).upper() < abs(fm).lower():
                d = (F(v)/F(u)).arg()
                if d.rad() < 0.5:
                    total += d; segs += 1
                    continue
            if v - u < minw:
                return dict(ok=False, edge=e, at=[u, v], segments=segs)
            m = (u + v)/2
            stack += [(m, v), (u, m)]
    n = total/(2*arb.pi())
    k = int(round(float(n.mid())))
    return dict(ok=bool(abs(n - k) < arb(1)/2), count=k, winding=n.str(8), segments=segs)


def box(x0, x1, y0, y1):
    return [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]


def secant(f, z0, z1, it=60):
    """Floating point refinement of a zero on the midpoints (only to place the boxes)."""
    for _ in range(it):
        f0, f1 = f(z0), f(z1)
        d = f1 - f0
        if abs(d).upper() == 0:
            break
        z2 = acb(z1 - f1*(z1 - z0)/d).mid()
        z0, z1 = z1, z2
        if abs(z1 - z0).upper() < arb(2)**-100:
            break
    return z1


def digits_ok(lo, hi, nd):
    """The interval [lo, hi] rounds to a single value at nd decimals; return it as a string, else None."""
    a = round(float(lo), nd); b = round(float(hi), nd)
    return ('%.*f' % (nd, a)) if a == b else None


def isolate(f, z, w=1e-10, nd=7, symmetric=False):
    x, y = float(z.real.mid()), float(z.imag.mid())
    if symmetric:
        r = winding(f, box(-w, w, y - w, y + w), n0=8)
        return dict(r, y=digits_ok(y - w, y + w, nd))
    r = winding(f, box(x - w, x + w, y - w, y + w), n0=8)
    return dict(r, x=digits_ok(x - w, x + w, nd), y=digits_ok(y - w, y + w, nd))


def atoms():
    c = arb(1)/10
    f = lambda s: (I*s).sinc() + c
    out = dict(count_box=winding(f, box(-6, 6, 0, 65), n0=256))
    assert (arb(6).sinh() > c*(6 + 65)), 'the bound for |Re s| >= 6'
    zs = []
    for y in (3.49906382, 5.679207796):                    # the zeros on the imaginary axis
        z = secant(f, acb(0, y), acb(0, y + 1e-6))
        zs.append(isolate(f, z, symmetric=True))
    for j in range(2, 11):                                 # the nine pairs, started from log(2 c y) + i y
        y = (2*j - 0.5)*math.pi
        s0 = acb(math.log(2*0.1*y), y)
        z = secant(f, s0, s0 + acb(1e-3, 1e-3))
        zs.append(isolate(f, z))
    out['zeros'] = zs
    out['ok'] = bool(out['count_box']['ok'] and out['count_box']['count'] == 20 and all(r['ok'] and r['count'] == 1 for r in zs))
    return out


def kmono(k):
    if k == 3:
        f = lambda s: 6*(s.sinh() - s)/s**3
        z0 = acb(2.7686783, 7.4976763)
    else:
        f = lambda s: 24*(s.cosh() - 1 - s*s/2)/s**4
        z0 = acb(4.5014572, 8.4247845)
    z = secant(f, z0, z0 + acb(1e-6, 1e-6))
    R = math.ceil((float(abs(z).mid()) + 0.01)*100)/100
    poly = [(1, 0), (R, 0), (R, R), (0, R), (0, 1)]
    cnt = winding(f, poly, n0=128)
    iso = isolate(f, z)
    # |L_k(s) - 1| <= k! sum_{j>=1} |s|^(2j)/(2j + k)! < 1 for |s| <= 1
    tail = sum(math.factorial(k)/math.factorial(2*j + k) for j in range(1, 30))
    return dict(k=k, R=R, count_square=cnt, zero=iso, modulus=abs(z).str(10), taylor_tail=tail,
                ok=bool(cnt['ok'] and cnt['count'] == 1 and iso['ok'] and iso['count'] == 1 and tail < 1))


def power(a):
    b = arb(a) - 1
    f = lambda z: z.sinc() - ((I*z).hypgeom_1f1(b + 1, b + 2) + (-I*z).hypgeom_1f1(b + 1, b + 2))/(2*(b + 1))
    out = dict(a=a)
    out['count_box'] = winding(f, box(0.3, 60, -12, 12), n0=256)
    v2pi = f(acb(2*arb.pi()))
    out['F_2pi'] = v2pi.real.str(8)
    out['sign_2pi'] = 1 if v2pi.real > 0 else (-1 if v2pi.real < 0 else 0)
    if float(a) > 2:
        signs = []
        for i in range(29851):
            v = f(acb(arb('0.3') + arb(i)/500)).real
            if v > 0:
                signs.append(1)
            elif v < 0:
                signs.append(-1)
        out['sign_changes'] = sum(1 for s, t in zip(signs, signs[1:]) if s != t)
        out['real_zeros'] = 18 if out['sign_changes'] >= 18 and out['count_box']['count'] == 18 else None
        out['ok'] = bool(out['count_box']['ok'] and out['real_zeros'] == 18 and out['sign_2pi'] == -1)
    else:
        out['count_upper'] = winding(f, box(0.3, 60, 0.05, 12), n0=256)
        up = out['count_upper']
        out['real_zeros'] = 0 if up['ok'] and up['count'] == 9 and out['count_box']['count'] == 18 else None
        out['ok'] = bool(out['count_box']['ok'] and out['real_zeros'] == 0 and out['sign_2pi'] == 1)
    print(json.dumps(out), flush=True)
    return out


def constants():
    g = lambda y: y.sin() - y*y.cos()                     # zero where tan y = y
    lo, hi = arb('4.49340945'), arb('4.49340946')
    ok = bool(g(lo)*g(hi) < 0)
    y = lo.union(hi)
    return dict(tan_root=y.str(12), bracket_ok=ok, sinc_min=(y.sin()/y).str(10), arccosh_3_2=arb(1.5).acosh().str(12),
                sqrt10_minus_2=(arb(10).sqrt() - 2).str(12))


if __name__ == '__main__':
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 15
    A = ['1.1', '1.3', '1.6', '1.9', '1.99', '2.01', '2.2', '2.5', '3', '4', '6']
    with Pool(workers) as pool:
        ra = pool.apply_async(atoms)
        rk = pool.map_async(kmono, [3, 4])
        rp = pool.map_async(power, A, chunksize=1)
        R = dict(atoms=ra.get(), kmonotone=rk.get(), power=rp.get(), constants=constants())
    digits = [z.get('y') for z in R['atoms']['zeros']] + [z.get('x') for z in R['atoms']['zeros'][2:]] + \
        [r['zero'].get(c) for r in R['kmonotone'] for c in ('x', 'y')]
    R['digits_decided'] = all(d is not None for d in digits)
    R['all_ok'] = bool(R['atoms']['ok'] and all(r['ok'] for r in R['kmonotone']) and all(r['ok'] for r in R['power'])
                       and R['constants']['bracket_ok'] and R['digits_decided'])
    R['seconds'] = time.time() - T0
    json.dump(R, open('fm_intcert.json', 'w'), indent=1)
    print('atoms', R['atoms']['ok'], R['atoms']['count_box'], [(z.get('x'), z.get('y')) for z in R['atoms']['zeros']])
    for r in R['kmonotone']:
        print('k monotone', r['k'], r['ok'], r['count_square'], r['zero'].get('x'), r['zero'].get('y'), r['modulus'])
    for r in R['power']:
        print('power a =', r['a'], r['ok'], 'box', r['count_box']['count'], 'real', r['real_zeros'], 'F(2 pi)', r['F_2pi'])
    print('constants', R['constants'])
    print('all certified:', R['all_ok'], '%.0f s' % R['seconds'])
