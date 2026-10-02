# Chapter 13, Section ch:fminterval: zeros of the characteristic functions of symmetric laws on [-1, 1].
# (1) The uniform law with an atom at the origin, (1-p) U(-1,1) + p delta_0 with c = p/(1-p) = 0.1. Its moment generating
#     function is (1-p)(sinh(s)/s + c). The zeros with 0 < Im s < 65: the two on the imaginary axis (sin y/y = -c), the
#     pairs +-x_j + i y_j off it (started from log(2 c y) + i y with y = (2j - 1/2) pi), and an argument principle count
#     on [-6, 6] x [1, 65] that shows the list is complete.
# (2) The symmetric k monotone kernels f_k(x) = (k/2)(1-|x|)_+^(k-1) for k = 3, 4. Their moment generating functions are
#     6 (sinh s - s)/s^3 and 24 (cosh s - 1 - s^2/2)/s^4. The zeros in the first quadrant with Im s < 40, ordered by
#     height, with argument principle counts on [0, 15] x [0, 40] and on the boxes below the first zero.
# (3) The densities proportional to 1 - |x|^b on [-1, 1], b = a - 1: zeros of F_b(z) = int_0^1 (1 - x^b) cos(zx) dx in the
#     box (0.3, 60) x (-12, 12) by the argument principle, and real zeros in (0.3, 60) by sign changes on a grid of step
#     0.002 (the real zeros near 2 pi m split by about 0.1 when b = 1.01). F_b(z) = sin(z)/z - (M(b+1, b+2, iz) +
#     M(b+1, b+2, -iz))/(2(b+1)), with M the confluent hypergeometric function.
# (4) Constants quoted in the text: arccosh(3/2) (spin one at p = 0.6), sqrt(10) - 2 (the moment test for Beta(a, 1)),
#     the minimum of sin(y)/y, and E G_{1/4} = 3 pi/8 - xi(5/2)/xi(1/2) (Section ch:fmexp).
# Usage: python3 fm_interval.py   Writes fm_interval.json.
import json, os, time
from multiprocessing import Pool
import mpmath as mp

mp.mp.dps = 20
t0 = time.time()


def winding(F, corners, n_per_unit=200, maxstep=0.6):
    """Number of zeros of F inside the polygon (counterclockwise corners), by the argument principle.
    Each edge is sampled with n_per_unit points per unit length; a step in arg larger than maxstep is refined."""
    total = mp.mpf(0)
    corners = [mp.mpc(x, y) for (x, y) in corners]    # mp.mpc((x, y)) would read the pair as a mantissa and an exponent
    for (a, b) in zip(corners, corners[1:] + corners[:1]):
        n = max(20, int(abs(b - a)*n_per_unit))
        assert n < 10**5, 'edge too long'
        ts = [mp.mpf(k)/n for k in range(n + 1)]
        prev = F(a)
        for t1, t2 in zip(ts, ts[1:]):
            stack = [(t1, t2)]
            while stack:
                u, v = stack.pop()
                fv = F(a + (b - a)*v)
                d = mp.arg(fv/prev)
                if abs(d) > maxstep and v - u > mp.mpf(10)**-12:
                    m = (u + v)/2
                    stack.append((m, v)); stack.append((u, m))
                    continue
                total += d
                prev = fv
    return int(mp.nint(total/(2*mp.pi))), float(total/(2*mp.pi))


def atoms(c=mp.mpf('0.1')):
    F = lambda s: (mp.sinh(s)/s if s != 0 else mp.mpf(1)) + c
    axis = []
    for lo, hi in [(mp.pi, mp.mpf(4.4934)), (mp.mpf(4.4934), 2*mp.pi)]:
        axis.append(mp.findroot(lambda y: mp.sin(y)/y + c, (lo + mp.mpf('0.01'), hi - mp.mpf('0.01')), solver='bisect'))
    off = []
    for j in range(2, 11):
        y = (2*j - mp.mpf(1)/2)*mp.pi
        z = mp.findroot(F, mp.mpc(mp.log(2*c*y), y))
        off.append(z)
    count = winding(F, [(-6, 1), (6, 1), (6, 65), (-6, 65)], n_per_unit=40)
    ymin = mp.findroot(lambda y: mp.diff(lambda t: mp.sin(t)/t, y), mp.mpf(4.49))
    return dict(c=float(c), imaginary_axis=[mp.nstr(y, 10) for y in axis],
                off_axis=[[mp.nstr(z.real, 8), mp.nstr(z.imag, 10)] for z in off],
                count_box=count, expected=2 + 2*len(off),
                sinc_min=[mp.nstr(ymin, 10), mp.nstr(mp.sin(ymin)/ymin, 10)])


def kmono():
    out = {}
    G = {3: lambda s: (mp.sinh(s) - s)/s**3 if abs(s) > 1e-6 else mp.mpf(1)/6,
         4: lambda s: (mp.cosh(s) - 1 - s**2/2)/s**4 if abs(s) > 1e-6 else mp.mpf(1)/24}
    for k in (3, 4):
        F = G[k]
        roots = []
        for j in range(0, 8):
            # for large s the equation is e^s/2 = s (k = 3) or e^s/2 = s^2/2 (k = 4): iterate s <- log(2 s) + 2 pi i j,
            # respectively s <- 2 log(s) + 2 pi i j, then polish on the exact function
            s = mp.mpc(1, 2*mp.pi*j + 1)
            for _ in range(60):
                s = (mp.log(2*s) if k == 3 else 2*mp.log(s)) + 2j*mp.pi*j
            try:
                z = mp.findroot(F, s)
            except Exception:
                continue
            if z.real < 0: z = -z.conjugate()
            if z.imag < 0: z = z.conjugate()
            if z.imag > 0.5 and z.imag < 40 and all(abs(z - w) > 1e-6 for w in roots): roots.append(z)
        roots.sort(key=lambda z: z.imag)
        cnt = winding(F, [(0, 0), (15, 0), (15, 40), (0, 40)], n_per_unit=40)
        below = winding(F, [(0, 0), (15, 0), (15, float(roots[0].imag) + 2), (0, float(roots[0].imag) + 2)], n_per_unit=60)
        out[k] = dict(roots=[[mp.nstr(z.real, 8), mp.nstr(z.imag, 8)] for z in roots], count_0_40=cnt,
                      count_below_first_plus2=below,
                      positive_on_axis_min=mp.nstr(min(F(mp.mpc(0, t)).real for t in mp.linspace(0.01, 60, 6000)), 6))
    return out


def power(a):
    mp.mp.dps = 20
    b = mp.mpf(a) - 1
    def F(z):
        if z == 0: return 1 - 1/(b + 1)
        return mp.sin(z)/z - (mp.hyp1f1(b + 1, b + 2, 1j*z) + mp.hyp1f1(b + 1, b + 2, -1j*z))/(2*(b + 1))
    zs = mp.linspace(mp.mpf('0.3'), 60, 29851)
    vals = [mp.re(F(z)) for z in zs]
    real = sum(1 for u, v in zip(vals, vals[1:]) if u*v < 0)
    box = winding(F, [(0.3, -12), (60, -12), (60, 12), (0.3, 12)], n_per_unit=40)
    near = []
    for m in (1, 2, 3):
        z0 = 2*mp.pi*m
        near.append(mp.nstr(F(z0).real, 6))
    return dict(a=a, b=float(b), real_zeros=real, box_count=box[0], winding=box[1], F_at_2pi_m=near)


if __name__ == '__main__':
    R = {}
    R['atoms'] = atoms(); print('atoms', R['atoms'], flush=True)
    R['kmonotone'] = kmono(); print('k monotone', R['kmonotone'], flush=True)
    A = [1.1, 1.3, 1.6, 1.9, 1.99, 2.01, 2.2, 2.5, 3, 4, 6]
    with Pool(min(len(A), int(os.environ.get('SLURM_CPUS_PER_TASK', '4')))) as P:
        R['power'] = P.map(power, A)
    for r in R['power']: print('power', r, flush=True)
    xi = lambda s: s*(s - 1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
    mp.mp.dps = 30
    R['constants'] = dict(arccosh_3_2=mp.nstr(mp.acosh(mp.mpf(3)/2), 12), sqrt10_minus_2=mp.nstr(mp.sqrt(10) - 2, 12),
                          EG_quarter=mp.nstr(3*mp.pi/8 - xi(mp.mpf(5)/2)/xi(mp.mpf(1)/2), 12),
                          xi52_over_xi12=mp.nstr(xi(mp.mpf(5)/2)/xi(mp.mpf(1)/2), 12))
    print('constants', R['constants'])
    R['seconds'] = time.time() - t0
    json.dump(R, open('fm_interval.json', 'w'), indent=1, default=str)
