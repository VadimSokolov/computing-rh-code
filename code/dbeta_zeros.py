# dbeta_zeros.py: the zeros of L(s, chi_4) below height 100 and the certified first ordinates used in Section
# sec:dbeta:beta (ch/ch16.tex), Theorem thm:dbeta:minall and Table tab:dbeta:minall (ch/basepoint.tex). Added for the book in
# October 2026, from the first author's note on the Dirichlet beta function in Thorin form (misc/reviews/nick-dbeta/), whose
# zeros and counts are recomputed here independently, in ball arithmetic with Arb (python-flint).
# For f = xi(s, chi_4) = (4/pi)^{(s+1)/2} Gamma((s+1)/2) L(s, chi_4), xi(s, chi_12) = (12/pi)^{s/2} Gamma(s/2) L(s, chi_12) and
# xi(s) = s(s-1)/2 pi^{-s/2} Gamma(s/2) zeta(s), all real on the critical line: (a) sign changes of f(1/2+it) on a grid of
# step 1/32, each sign certified (the real ball excludes zero), refined by bisection to 1e-12, each one a zero on the line;
# (b) the number N(T) of zeros with 0 < Im rho <= T, all multiplicities, by the argument principle: N(T) = (1/pi) times the
# change of arg f along 2 -> 2+iT -> 1/2+iT (f > 0 on the real axis; the rest of the boundary of [-1,2] x [0,T] follows from
# f(1-s) = f(s) and f(conj s) = conj f(s)). Each piece of the path is enclosed in a ball, f is evaluated on the ball, and
# the piece is halved until the enclosure lies in a disc that excludes zero; the change of argument over the piece is then
# the principal argument of the quotient of the end values, summed in ball arithmetic. When N(T) equals the number of sign
# changes, every zero below T is simple and on the critical line.
# Also 2 sum gamma^{-2} over the zeros of L(s, chi_4) below 100 with the main term of its tail, the first five ordinates of
# L(s, chi_12) (below height 14) with their gaps, and L(sigma, chi_4), L(sigma, chi_12) at sigma = 0.1, 0.2, ..., 0.9.
# Run on Hopper (about two minutes): python3 dbeta_zeros.py; writes dbeta_zeros.json in the current folder.
import json
import time

from flint import acb, arb, ctx, dirichlet_char

t0 = time.time()
out = {}


def rec(k, v):
    out[k] = v
    print(k, v, round(time.time() - t0, 1), flush=True)


ctx.prec = 128
CHI4 = dirichlet_char(4, 3)
CHI12 = dirichlet_char(12, 11)
PI = arb.pi()
LOG4PI = (arb(4) / PI).log()
LOG12PI = (arb(12) / PI).log()
LOGPI = PI.log()


def xi4(s):
    return ((s + 1) / 2 * LOG4PI).exp() * ((s + 1) / 2).gamma() * acb.dirichlet_l(s, CHI4)


def xi12(s):
    return (s / 2 * LOG12PI).exp() * (s / 2).gamma() * acb.dirichlet_l(s, CHI12)


def xiz(s):
    return s * (s - 1) / 2 * (-s / 2 * LOGPI).exp() * (s / 2).gamma() * s.zeta()


def sign_on_line(f, t):
    v = f(acb(arb(0.5), arb(t))).real
    if v > 0:
        return 1
    if v < 0:
        return -1
    return 0


def line_zeros(f, T, h=1 / 32):
    n = int(T / h)
    ts = [k * h for k in range(1, n + 1)]
    sg = []
    for t in ts:
        s = sign_on_line(f, t)
        if s == 0:
            raise RuntimeError('sign not certified at %r' % t)
        sg.append(s)
    roots = []
    for k in range(len(ts) - 1):
        if sg[k] != sg[k + 1]:
            a, b, sa = ts[k], ts[k + 1], sg[k]
            while b - a > 1e-12:
                m = (a + b) / 2
                sm = sign_on_line(f, m)
                if sm == 0:
                    break
                if sm == sa:
                    a = m
                else:
                    b = m
            roots.append((a + b) / 2)
    return roots


def ball(z0, z1):
    # complex ball containing the segment [z0, z1] (z0, z1 given as pairs of exact binary floats)
    mr, mi = (z0[0] + z1[0]) / 2, (z0[1] + z1[1]) / 2
    rr, ri = abs(z1[0] - z0[0]) / 2, abs(z1[1] - z0[1]) / 2
    return acb(arb(mr, rr), arb(mi, ri))


def pt(z):
    return acb(arb(z[0]), arb(z[1]))


def seg_arg(f, z0, z1, depth=0):
    # change of arg f along the segment [z0, z1]: accepted when the enclosure of f over the segment lies in a disc
    # |w - m| <= r with |m| > 1.1 r, so that the argument varies by less than 2 arcsin(1/1.1) < pi along the segment
    W = f(ball(z0, z1))
    m = W.mid()
    r = (W.real.rad() ** 2 + W.imag.rad() ** 2).sqrt()
    if abs(m) > 1.1 * r:
        return (f(pt(z1)) / f(pt(z0))).arg(), 1
    if depth > 40:
        raise RuntimeError('piece not resolved near %r' % (z0,))
    zm = ((z0[0] + z1[0]) / 2, (z0[1] + z1[1]) / 2)
    a1, n1 = seg_arg(f, z0, zm, depth + 1)
    a2, n2 = seg_arg(f, zm, z1, depth + 1)
    return a1 + a2, n1 + n2


def count(f, T, h=1 / 16):
    tot = arb(0)
    pieces = 0
    n = int(T / h)
    for k in range(n):
        a, p = seg_arg(f, (2.0, k * h), (2.0, (k + 1) * h))
        tot += a
        pieces += p
    m = int(1.5 / h)
    for k in range(m):
        a, p = seg_arg(f, (2.0 - k * h, float(T)), (2.0 - (k + 1) * h, float(T)))
        tot += a
        pieces += p
    N = tot / PI
    return N, pieces


res = {}
for name, f, T in [('chi4', xi4, 100), ('chi12', xi12, 14), ('xi', xiz, 20)]:
    roots = line_zeros(f, T)
    N, pieces = count(f, T)
    Nint = int(round(float(N.mid())))
    certified = bool((N - Nint).abs_upper() < 0.25)
    gaps = [roots[k + 1] - roots[k] for k in range(len(roots) - 1)]
    res[name] = {'T': T, 'sign_changes': len(roots), 'N(T)_ball': N.str(10), 'N(T)': Nint, 'count_certified': certified,
                 'all_on_line_and_simple': bool(certified and Nint == len(roots)), 'pieces': pieces,
                 'ordinates': ['%.10f' % r for r in roots],
                 'largest_gap': '%.6f' % max(gaps) if gaps else None,
                 'largest_gap_after': '%.6f' % roots[max(range(len(gaps)), key=lambda k: gaps[k])] if gaps else None,
                 'gaps': ['%.6f' % x for x in gaps]}
    rec(name, res[name])

# the main term of the zero count of L(s, chi_4) at T = 100, (T/2pi) log(2T/(pi e)) + 1/8
T = arb(100)
rec('chi4_main_term_100', (T / (2 * PI) * (2 * T / (PI * arb(1).exp())).log() + arb(1) / 8).str(8))
# 2 sum of gamma^{-2} over the ordinates below 100, against V = (log xi)''(1/2, chi_4) = sum over all zeros of gamma^{-2} (under
# GRH), and the main term (1/pi)(log(2T/pi)+1)/T of 2 sum_{gamma>T} gamma^{-2}, from dN(t) ~ (1/2pi) log(2t/pi) dt, at T = 100
g = [arb(x) for x in res['chi4']['ordinates']]
rec('chi4_sum_inverse_squares', {'2 sum_{gamma<100} gamma^-2': (2 * sum(1 / x ** 2 for x in g)).str(8),
                                 'tail main term at T=100': (((2 * T / PI).log() + 1) / (PI * T)).str(4)})
# L(sigma, chi) on the real axis, sigma = 0.1, ..., 0.9, in ball arithmetic (positive, as the positive kernels require)
rec('L_real', {name: {'%.1f' % (k / 10): acb.dirichlet_l(acb(arb(k) / 10), chi).real.str(10) for k in range(1, 10)}
               for name, chi in [('chi4', CHI4), ('chi12', CHI12)]})
json.dump(out, open('dbeta_zeros.json', 'w'), indent=1)
print('done', round(time.time() - t0, 1))
