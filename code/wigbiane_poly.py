# wigbiane_poly.py: the collision bounds of the section on Laguerre's quantity and the collisions of zeros along the de Bruijn
# Newman flow (ch/ch21.tex, after Section sec:ch21:hankel; draft misc/reviews/nick-wigner/merge/wigbiane.tex) tested on
# random polynomials with real zeros. Added for the book in October 2026, from the first author's note on the semicircle law and
# the Thorin framework (misc/reviews/nick-wigner/, version 7, Section 8, which reports 379 polynomials of degree 3 to 8); the
# computation here is independent of the note's scripts.
# For p of degree n with simple real zeros, p_s = e^{s D^2} p = sum_k s^k p^{(2k)}/k! solves d p_s/ds = p_s'' (H_t with
# s = -t). Four families, n = 3, ..., 12, 60 polynomials each: zeros uniform on [-1, 1]; Gaussian; uniform with one zero added
# at distance d v from a random zero (d = 2/(n - 1), v uniform on [0.02, 0.3]: a Lehmer type pair); and a perturbed lattice
# j + 0.45 U_j. The zeros are mapped affinely onto [-1, 1] (collision times scale with the square of the length and Laguerre's
# quantity with its inverse square, so every ratio tested is unchanged). The zeros are followed by numpy.roots with steps at
# most 0.01 g^2 (g the least gap between adjacent real zeros); when g falls below 1/100 of its value at s = 0 the collision is
# located as the double zero p_s = p_s' = 0 by Newton's method in (z, s) at 40 digits (mpmath), with the polynomial expanded
# about the midpoint of the closing gap from its zeros at s = 0; the pair continues as a conjugate pair, and the flow runs until
# at most one real zero is left. On every gap between zeros that are adjacent and real at s = 0, at every step, the script
# computes Lambda = -(log|p_s|)'' at the critical point c (L) and its minimum m over the gap, attained at u, and from
# x_j = 1/(u - z_j) the share P = x_k^2 + x_{k+1}^2 of the two zeros bounding the gap, the rest R = m - P (summed directly),
# S_4 = Lambda''(u)/6, eta_dir = 6 S_4/m^2 - 3, kappa = 6 S_4 - 3 P^2 = 3 (x_k^2 - x_{k+1}^2)^2 + 6 (S_4 of the other zeros),
# W = sum |u - z|^{-4} over the nonreal zeros, B_p = sum b^2/|u - a|^p over the conjugate pairs a +- ib, and
# eta = (36 m B_4 + 108 B_4^2 + 72 (m + 6 B_4)^{1/2} B_5 + 120 B_6)/m^2. For the first collision (all zeros real before it) it
# checks s_1 L(0) >= 1 and s_1 m(0) >= 1 (critical point and Biane forms), s_1 <= s* with rho = sup R and W = 0 (the
# deadline), and, on every gap before s_1, L(s) <= L(0)/(1 - s L(0)) and m(s) <= m(0)/(1 - s m(0)), in double precision at
# every step and again at 40 digits (zeros refined by Newton's method) at the 20 states per polynomial and ratio where the
# double precision ratio comes closest to 1, if within 10^-4 (near a collision numpy.roots loses about half the digits of the
# gap, which matters at that level). For every later
# collision of two zeros adjacent at s = 0 it checks the local bound s_c >= 1/((1 + eta_bar) m(0)), eta_bar the largest eta
# before s_c (and the same with eta_dir), and the deadline s_c <= s* with rho = sup R and W_bar = sup W (and with 6 W_bar
# replaced by the largest -kappa). The suprema are taken over the steps, so they can only be smaller than the true suprema,
# which makes every check stricter. At s = 0 it compares dL/ds and dm/ds, by central differences with step 10^{-12} g^2 at 40
# digits (zeros at s = +-h by Newton's method), with 6 S_4 - 2 S_2^2 - 4 S_3^2/S_2 at c and 6 S_4 - 2 S_2^2 at u. Run from
# any folder: python3 wigbiane_poly.py; writes wigbiane_poly.json there. Uses SLURM_CPUS_PER_TASK processes (default 1).
import json
import math
import multiprocessing as mpr
import os
import time

import mpmath as mp
import numpy as np
from numpy.polynomial import polynomial as Pl

t0 = time.time()
out = {}
DPS = 40
STEP_CAP = 400000
RECHECK = 20


def rec(k, v):
    out[k] = v
    print(k, v, round(time.time() - t0, 1), flush=True)
    json.dump(out, open('wigbiane_poly.json', 'w'), indent=1)


def flow_coeffs(c, s):
    # coefficients (lowest first) of e^{s D^2} p, p with coefficients c (numpy)
    res = np.zeros_like(c)
    d = c.copy()
    k, fac = 0, 1.0
    while len(d) > 0:
        res[:len(d)] += d * (s ** k / fac)
        d = Pl.polyder(d, 2) if len(d) > 2 else np.zeros(0)
        k += 1
        fac *= k
    return res


# 40 digit helpers (mpmath), coefficient lists lowest first
def mp_fromroots(roots):
    c = [mp.mpf(1)]
    for r in roots:
        n = [mp.mpf(0)] * (len(c) + 1)
        for i, a in enumerate(c):
            n[i + 1] += a
            n[i] -= r * a
        c = n
    return c


def mp_der(c):
    return [c[i] * i for i in range(1, len(c))]


def mp_flow(c, s):
    res = [mp.mpf(0)] * len(c)
    d = list(c)
    k, fac = 0, mp.mpf(1)
    while d:
        f = s ** k / fac
        for i in range(len(d)):
            res[i] += d[i] * f
        d = mp_der(mp_der(d)) if len(d) > 2 else []
        k += 1
        fac *= k
    return res


def mp_ev(c, x):
    v = mp.mpf(0)
    for a in reversed(c):
        v = v * x + a
    return v


def double_zero(roots0, zm, s0, g):
    # Newton in (w, s) on p_s(zm + w) = p_s'(zm + w) = 0, p expanded about zm from its zeros at s = 0
    mp.mp.dps = DPS
    zm = mp.mpf(zm)
    q = mp_fromroots([mp.mpf(r) - zm for r in roots0])
    w, s = mp.mpf(0), mp.mpf(s0)
    for _ in range(60):
        cs = mp_flow(q, s)
        d1 = mp_der(cs)
        d2 = mp_der(d1)
        d3 = mp_der(d2)
        F1, F2 = mp_ev(cs, w), mp_ev(d1, w)
        a11, a12, a22 = F2, mp_ev(d2, w), mp_ev(d3, w)
        det = a11 * a22 - a12 * a12
        dw = (a22 * F1 - a12 * F2) / det
        dsn = (-a12 * F1 + a11 * F2) / det
        w -= dw
        s -= dsn
        if abs(dw) < mp.mpf(10) ** (10 - DPS) * g and abs(dsn) < mp.mpf(10) ** (10 - DPS) * abs(s):
            return float(zm + w), float(s), True
    return float(zm + w), float(s), False


def mp_gap_values(zs, ia, ib, u0, c0):
    # L = Lambda(c) and m = Lambda(u) at 40 digits for real zeros zs (mpf), starting from the double precision u0, c0
    g = zs[ib] - zs[ia]
    u = mp.mpf(u0)
    for _ in range(100):
        s3 = s4 = mp.mpf(0)
        for z in zs:
            x = 1 / (u - z)
            x2 = x * x
            s3 += x2 * x
            s4 += x2 * x2
        du = s3 / (3 * s4)
        u += du
        if abs(du) < mp.mpf(10) ** (5 - DPS) * g:
            break
    m = sum(1 / (u - z) ** 2 for z in zs)
    c = mp.mpf(c0)
    for _ in range(100):
        s1 = s2 = mp.mpf(0)
        for z in zs:
            x = 1 / (c - z)
            s1 += x
            s2 += x * x
        dc = s1 / s2
        c += dc
        if abs(dc) < mp.mpf(10) ** (5 - DPS) * g:
            break
    L = sum(1 / (c - z) ** 2 for z in zs)
    return L, m, c


def parts(x, zs):
    # S_1, S_2, S_3, S_4 at the real point x (real parts), zs a list of complex zeros
    s1 = s2 = s3 = s4 = 0j
    for z in zs:
        r = 1 / (x - z)
        r2 = r * r
        s1 += r
        s2 += r2
        s3 += r2 * r
        s4 += r2 * r2
    return s1.real, s2.real, s3.real, s4.real


def lam(x, zs):
    s2 = 0j
    for z in zs:
        r = 1 / (x - z)
        s2 += r * r
    return s2.real


def minimise_gap(lo, hi, zs, u0=None):
    # global minimiser of Lambda on (lo, hi): a grid of 24 points (and the previous minimiser), then safeguarded Newton on
    # Lambda' = -2 S_3 with Lambda'' = 6 S_4
    xs = [lo + (hi - lo) * (0.5 - 0.5 * math.cos(math.pi * j / 25)) for j in range(1, 25)]
    if u0 is not None and lo < u0 < hi:
        xs.append(u0)
        xs.sort()
    vals = [lam(x, zs) for x in xs]
    i = min(range(len(xs)), key=lambda j: vals[j])
    a = xs[i - 1] if i > 0 else lo + 1e-12 * (hi - lo)
    b = xs[i + 1] if i + 1 < len(xs) else hi - 1e-12 * (hi - lo)
    u = xs[i]
    for _ in range(80):
        _, s2, s3, s4 = parts(u, zs)
        l1, l2 = -2 * s3, 6 * s4
        if l1 > 0:
            b = u
        else:
            a = u
        un = u - l1 / l2 if l2 > 0 else (a + b) / 2
        if not (a < un < b):
            un = (a + b) / 2
        if abs(un - u) < 1e-15 * (1 + abs(u)):
            return un
        u = un
    return u


def critical(lo, hi, zs, c0=None):
    # the zero of q = S_1 in (lo, hi); q falls from +inf to -inf
    a, b = lo + 1e-13 * (hi - lo), hi - 1e-13 * (hi - lo)
    c = c0 if (c0 is not None and a < c0 < b) else (lo + hi) / 2
    for _ in range(200):
        q, s2, _, _ = parts(c, zs)
        if q > 0:
            a = c
        else:
            b = c
        cn = c + q / s2 if s2 > 0 else (a + b) / 2
        if not (a < cn < b):
            cn = (a + b) / 2
        if abs(cn - c) < 1e-15 * (1 + abs(c)):
            return cn
        c = cn
    return c


def gap_state(ia, ib, zs, prev=None):
    lo, hi = zs[ia].real, zs[ib].real
    u = minimise_gap(lo, hi, zs, prev['u'] if prev else None)
    xa, xb = 1 / (u - lo), 1 / (u - hi)
    P = xa * xa + xb * xb
    R = S4r = B4 = B5 = B6 = W = 0.0
    for j, z in enumerate(zs):
        if j == ia or j == ib:
            continue
        r = 1 / (u - z)
        r2 = r * r
        R += r2.real
        S4r += (r2 * r2).real
        if z.imag != 0:
            W += abs(r2) ** 2
            if z.imag > 0:
                d = max(abs(u - z.real), 1e-300)
                B4 += z.imag ** 2 / d ** 4
                B5 += z.imag ** 2 / d ** 5
                B6 += z.imag ** 2 / d ** 6
    m = P + R
    kappa = 3 * (xa * xa - xb * xb) ** 2 + 6 * S4r
    if m > 0:
        eta = (36 * m * B4 + 108 * B4 ** 2 + 72 * math.sqrt(max(m + 6 * B4, 0.0)) * B5 + 120 * B6) / m ** 2
        etad = (kappa - 6 * P * R - 3 * R * R) / m ** 2
    else:
        eta = etad = float('inf')
    c = critical(lo, hi, zs, prev['c'] if prev else None)
    L = parts(c, zs)[1]
    return {'u': u, 'm': m, 'P': P, 'R': R, 'S4': (xa ** 4 + xb ** 4 + S4r), 'eta_dir': etad, 'kappa': kappa, 'W': W,
            'eta': eta, 'c': c, 'L': L}


def sstar(m0, rho, Wb):
    root = math.sqrt(9 * rho * rho + 6 * Wb)
    r1, r2 = 3 * rho + root, 3 * rho - root
    if m0 <= r1:
        return None
    if r1 - r2 < 1e-12 * m0:
        return 1 / m0
    return math.log((m0 - r2) / (m0 - r1)) / (r1 - r2)


def make_zeros(fam, n, rng):
    while True:
        if fam == 'uniform':
            z = rng.uniform(-1, 1, n)
        elif fam == 'gauss':
            z = rng.standard_normal(n)
        elif fam == 'lehmer':
            z = rng.uniform(-1, 1, n - 1)
            j = int(rng.integers(0, n - 1))
            z = np.append(z, z[j] + (2.0 / (n - 1)) * rng.uniform(0.02, 0.3) * rng.choice([-1, 1]))
        else:
            z = np.arange(1, n + 1) + 0.45 * rng.uniform(-1, 1, n)
        z = np.sort(z)
        if np.diff(z).min() > 1e-3 * (z[-1] - z[0]):
            return (z - (z[0] + z[-1]) / 2) * (2 / (z[-1] - z[0]))


def match(old, new):
    # assign the roots new to the zeros old by nearest neighbour (steps are small against the gaps)
    res = [0j] * len(old)
    free = list(new)
    for k, z in enumerate(old):
        j = min(range(len(free)), key=lambda i: abs(free[i] - z))
        res[k] = complex(free.pop(j))
    return res


def derivative_checks(z0, gaps0):
    # central differences at 40 digits for dL/ds and dm/ds at s = 0, against the formulas of the evolution lemma
    mp.mp.dps = DPS
    zr = [mp.mpf(float(x)) for x in z0]
    cm = mp_fromroots(zr)
    g0 = min(zr[i + 1] - zr[i] for i in range(len(zr) - 1))
    h = mp.mpf(10) ** -12 * g0 * g0
    roots = {}
    for sg in (1, -1):
        cs = mp_flow(cm, sg * h)
        d1 = mp_der(cs)
        rs = []
        for z in zr:
            for _ in range(12):
                z = z - mp_ev(cs, z) / mp_ev(d1, z)
            rs.append(z)
        roots[sg] = rs
    errL = errm = 0.0
    for (ia, ib), st in gaps0.items():
        vals = {}
        for sg in (1, -1):
            vals[sg] = mp_gap_values(roots[sg], ia, ib, st['u'], st['c'])
        L0, m0, c0 = mp_gap_values(zr, ia, ib, st['u'], st['c'])
        xs = [1 / (c0 - z) for z in zr]
        S2, S3, S4 = sum(x ** 2 for x in xs), sum(x ** 3 for x in xs), sum(x ** 4 for x in xs)
        dL = (vals[1][0] - vals[-1][0]) / (2 * h)
        errL = max(errL, float(abs(dL - (6 * S4 - 2 * S2 ** 2 - 4 * S3 ** 2 / S2)) / S2 ** 2))
        # minimiser at 40 digits for the dm/ds formula
        u = mp.mpf(st['u'])
        for _ in range(100):
            ys = [1 / (u - z) for z in zr]
            du = sum(y ** 3 for y in ys) / (3 * sum(y ** 4 for y in ys))
            u += du
            if abs(du) < mp.mpf(10) ** (5 - DPS) * g0:
                break
        ys = [1 / (u - z) for z in zr]
        m2, m4 = sum(y ** 2 for y in ys), sum(y ** 4 for y in ys)
        dm = (vals[1][1] - vals[-1][1]) / (2 * h)
        errm = max(errm, float(abs(dm - (6 * m4 - 2 * m2 ** 2)) / m2 ** 2))
    return errL, errm


def run_poly(args):
    fam, n, seed = args
    rng = np.random.default_rng(seed)
    z0 = make_zeros(fam, n, rng)
    c = Pl.polyfromroots(z0)
    zs = [complex(z) for z in z0]
    real = [True] * n
    g0 = float(np.diff(z0).min())
    s = 0.0
    smax = 80.0
    gaps = dict(((i, i + 1), {'hist': [], 'alive': True, 'prev': None}) for i in range(n - 1))
    gaps0 = {key: gap_state(key[0], key[1], zs) for key in gaps}
    errL, errm = derivative_checks(z0, gaps0)
    collisions = []
    cand = {'m': [], 'L': []}
    steps = 0
    while s < smax and sum(real) >= 2 and steps < STEP_CAP:
        steps += 1
        order = sorted([k for k in range(n) if real[k]], key=lambda k: zs[k].real)
        rz = [zs[k].real for k in order]
        gl = [rz[q + 1] - rz[q] for q in range(len(rz) - 1)]
        for q in range(len(order) - 1):
            key = (order[q], order[q + 1])
            if key in gaps and gaps[key]['alive']:
                st = gap_state(key[0], key[1], zs, gaps[key]['prev'])
                gaps[key]['prev'] = st
                # compact history: s, m, L, R, eta, eta_dir, W, kappa
                gaps[key]['hist'].append((s, st['m'], st['L'], st['R'], st['eta'], st['eta_dir'], st['W'], st['kappa']))
                if not collisions and s > 0:
                    # states before the first collision where L or m comes within 10^-4 of its bound: kept for a
                    # recheck at 40 digits (near a collision the double precision zeros of numpy.roots are not
                    # accurate enough to test the bound at that level)
                    h0 = gaps[key]['hist'][0]
                    for kind, val, v0 in (('m', st['m'], h0[1]), ('L', st['L'], h0[2])):
                        if s * v0 < 1 and val * (1 - s * v0) / v0 > 1 - 1e-4:
                            cand[kind].append((val * (1 - s * v0) / v0, s, key, list(zs), st['u'], st['c']))
                            if len(cand[kind]) > 200:
                                cand[kind].sort(key=lambda x: -x[0])
                                del cand[kind][RECHECK:]
        q0 = min(range(len(gl)), key=lambda q: gl[q])
        gmin = gl[q0]
        if gmin < 1e-2 * g0:
            ka, kb = order[q0], order[q0 + 1]
            mid = (rz[q0] + rz[q0 + 1]) / 2
            zc, sc, ok = double_zero(z0, mid, s + gmin * gmin / 8, gmin)
            if not (ok and s - 0.01 * gmin * gmin < sc < s + gmin * gmin and abs(zc - mid) < gmin):
                return {'family': fam, 'n': n, 'seed': seed, 'failed': 'double zero at s = %g' % s}
            collisions.append({'ids': (ka, kb), 's': sc, 'pairs': sum(1 for k in range(n) if not real[k]) // 2})
            for key in gaps:
                if gaps[key]['alive'] and (ka in key or kb in key):
                    gaps[key]['alive'] = False
            s = sc + 1e-4 * gmin * gmin
            roots = list(np.roots(flow_coeffs(c, s)[::-1]))
            others = [k for k in range(n) if k not in (ka, kb)]
            new = match([zs[k] for k in others], roots)
            for k, z in zip(others, new):
                zs[k] = z
            for z in new:
                for i_, r in enumerate(roots):
                    if r == z:
                        roots.pop(i_)
                        break
            roots.sort(key=lambda r: r.imag)
            zs[ka], zs[kb] = complex(roots[0]), complex(roots[1])
            real[ka] = real[kb] = False
            for k in range(n):
                if real[k]:
                    zs[k] = complex(zs[k].real)
            continue
        ds = min(0.01 * gmin * gmin, 0.02 * s + 0.001 * g0 * g0)
        for attempt in range(60):
            new = match(zs, list(np.roots(flow_coeffs(c, s + ds)[::-1])))
            if all(abs(new[k].imag) < 1e-7 * (1 + abs(new[k].real)) for k in range(n) if real[k]):
                break
            ds /= 4
        else:
            return {'family': fam, 'n': n, 'seed': seed, 'failed': 'step at s = %g' % s}
        zs = [complex(new[k].real) if real[k] else new[k] for k in range(n)]
        s += ds
    if steps >= STEP_CAP or sum(real) >= 2:
        return {'family': fam, 'n': n, 'seed': seed, 'failed': 'unfinished at s = %g after %d steps' % (s, steps)}
    res = {'family': fam, 'n': n, 'seed': seed, 'collisions': len(collisions), 'steps': steps,
           'checks': {'dL/ds error': errL, 'dm/ds error': errm}}
    c1 = collisions[0]
    H = gaps[c1['ids']]['hist']
    m0, L0 = H[0][1], H[0][2]
    rho = max(x[3] for x in H)
    ss = sstar(m0, rho, 0.0)
    res['first'] = {'s1': c1['s'], 's1 L(0)': c1['s'] * L0, 's1 m(0)': c1['s'] * m0, 's1/s*': c1['s'] / ss if ss else None,
                    '6R(0)/m(0)': 6 * H[0][3] / m0}
    worstL = worstm = 0.0
    for gd in gaps.values():
        Hh = [x for x in gd['hist'] if x[0] < c1['s']]
        if not Hh:
            continue
        L0g, m0g = Hh[0][2], Hh[0][1]
        for x in Hh:
            if x[0] * L0g < 1:
                worstL = max(worstL, x[2] * (1 - x[0] * L0g) / L0g)
            if x[0] * m0g < 1:
                worstm = max(worstm, x[1] * (1 - x[0] * m0g) / m0g)
    res['L ratio'] = worstL
    res['m ratio'] = worstm
    res['rechecked'] = 0
    res['L ratio (40 digits)'] = res['m ratio (40 digits)'] = None
    if cand['m'] or cand['L']:
        mp.mp.dps = DPS
        zr0 = [mp.mpf(float(x)) for x in z0]
        cm = mp_fromroots(zr0)
        base = {}
        for kind in ('m', 'L'):
            cand[kind].sort(key=lambda x: -x[0])
            vals = []
            for (r, s_, key, zsnap, u_, c_) in cand[kind][:RECHECK]:
                if key not in base:
                    base[key] = mp_gap_values(zr0, key[0], key[1], gaps0[key]['u'], gaps0[key]['c'])
                sm = mp.mpf(s_)
                cs = mp_flow(cm, sm)
                d1 = mp_der(cs)
                rts = []
                for z in zsnap:
                    x = mp.mpf(z.real)
                    for _ in range(40):
                        x = x - mp_ev(cs, x) / mp_ev(d1, x)
                    rts.append(x)
                L_, m_, _ = mp_gap_values(rts, key[0], key[1], u_, c_)
                L0_, m0_, _ = base[key]
                v, v0 = (m_, m0_) if kind == 'm' else (L_, L0_)
                vals.append(float(v * (1 - sm * v0) / v0))
                res['rechecked'] += 1
            if vals:
                res[kind + ' ratio (40 digits)'] = max(vals)
    later = []
    for col in collisions[1:]:
        key = col['ids']
        if key not in gaps or not gaps[key]['hist']:
            continue
        H = gaps[key]['hist']
        m0 = H[0][1]
        etab = max(x[4] for x in H)
        etad = max(x[5] for x in H)
        rho = max(x[3] for x in H)
        Wb = max(x[6] for x in H)
        kb = max(0.0, max(-x[7] for x in H))
        upB, upD = sstar(m0, rho, Wb), sstar(m0, rho, kb / 6)
        later.append({'s_c m(0)': col['s'] * m0, 'eta_bar': etab, 'eta_dir_bar': etad,
                      'low ratio (eta)': col['s'] * (1 + etab) * m0,
                      'low ratio (eta_dir)': col['s'] * (1 + max(etad, 0.0)) * m0,
                      'up ratio (W)': col['s'] / upB if upB else None, 'up ratio (kappa)': col['s'] / upD if upD else None,
                      'complex pairs at s_c': col['pairs']})
    res['later'] = later
    return res


if __name__ == '__main__':
    nproc = int(os.environ.get('SLURM_CPUS_PER_TASK', '1'))
    pool = mpr.get_context('fork').Pool(nproc)
    jobs = []
    seed = 7000
    for fam in ('uniform', 'gauss', 'lehmer', 'lattice'):
        for n in range(3, 13):
            for r in range(60):
                seed += 1
                jobs.append((fam, n, seed))
    res = pool.map(run_poly, jobs, chunksize=2)
    rec('polynomials', len(res))
    rec('failures', [r for r in res if 'failed' in r])
    res = [r for r in res if 'failed' not in r]
    firsts = [r['first'] for r in res]
    rec('first collisions', len(firsts))
    rec('least s1 L(0) (critical point form, needs >= 1)', min(f['s1 L(0)'] for f in firsts))
    rec('least s1 m(0) (Biane form, needs >= 1)', min(f['s1 m(0)'] for f in firsts))
    rec('quartiles of s1 L(0)', [float(np.quantile([f['s1 L(0)'] for f in firsts], q)) for q in (0.25, 0.5, 0.75)])
    rec('quartiles of s1 m(0)', [float(np.quantile([f['s1 m(0)'] for f in firsts], q)) for q in (0.25, 0.5, 0.75)])
    sr = [f['s1/s*'] for f in firsts if f['s1/s*'] is not None]
    rec('first collisions with m(0) > 6 rho (deadline defined)', len(sr))
    rec('largest s1/s* (deadline, needs <= 1)', max(sr))
    rec('least s1/s* (how tight the deadline is at worst)', min(sr))
    rec('largest L(s) over its Riccati bound before s1, double precision', max(r['L ratio'] for r in res))
    rec('largest m(s) over its Biane bound before s1, double precision', max(r['m ratio'] for r in res))
    rec('states within 10^-4 of a bound rechecked at 40 digits (at most %d per polynomial and kind)' % RECHECK,
        sum(r['rechecked'] for r in res))
    l40 = [r['L ratio (40 digits)'] for r in res if r['L ratio (40 digits)'] is not None]
    m40 = [r['m ratio (40 digits)'] for r in res if r['m ratio (40 digits)'] is not None]
    rec('largest ratios at 40 digits at the rechecked states (L, m; both need <= 1)',
        [max(l40) if l40 else None, max(m40) if m40 else None])
    wr = max(res, key=lambda r: r['m ratio'])
    rec('polynomial with the largest double precision m ratio, and its ratio at 40 digits',
        {'family': wr['family'], 'n': wr['n'], 'seed': wr['seed'], 'double': wr['m ratio'],
         '40 digits': wr['m ratio (40 digits)'], 's1 m(0)': wr['first']['s1 m(0)']})
    by = {}
    for r in res:
        b = by.setdefault(r['family'], {'count': 0, 'least s1 L(0)': 9.0, 'least s1 m(0)': 9.0, 'all': []})
        b['count'] += 1
        b['least s1 L(0)'] = min(b['least s1 L(0)'], r['first']['s1 L(0)'])
        b['least s1 m(0)'] = min(b['least s1 m(0)'], r['first']['s1 m(0)'])
        b['all'].append(r['first']['s1 m(0)'])
    for b in by.values():
        b['median s1 m(0)'] = float(np.median(b.pop('all')))
    rec('by family', by)
    lat = [x for r in res for x in r['later']]
    rec('later collisions of pairs adjacent at s = 0', len(lat))
    rec('later collisions with at least one conjugate pair present', sum(1 for x in lat if x['complex pairs at s_c'] > 0))
    rec('least s_c (1 + eta_bar) m(0) (local bound, needs >= 1)', min(x['low ratio (eta)'] for x in lat))
    rec('least s_c (1 + max(eta_dir_bar, 0)) m(0) (direct form, needs >= 1)', min(x['low ratio (eta_dir)'] for x in lat))
    rec('least s_c m(0) among later collisions, and how many are below 1',
        [min(x['s_c m(0)'] for x in lat), sum(1 for x in lat if x['s_c m(0)'] < 1)])
    rec('later collisions with eta_dir <= 0 at every step (then s_c m(0) >= 1 follows)',
        sum(1 for x in lat if x['eta_dir_bar'] <= 0))
    rec('largest eta_bar and eta_dir_bar', [max(x['eta_bar'] for x in lat), max(x['eta_dir_bar'] for x in lat)])
    rec('median eta_bar and eta_dir_bar', [float(np.median([x['eta_bar'] for x in lat])),
                                            float(np.median([x['eta_dir_bar'] for x in lat]))])
    uw = [x['up ratio (W)'] for x in lat if x['up ratio (W)'] is not None]
    uq = [x['up ratio (kappa)'] for x in lat if x['up ratio (kappa)'] is not None]
    rec('later: deadline defined (W form, kappa form)', [len(uw), len(uq)])
    rec('largest s_c/s* (W form; kappa form; both need <= 1)', [max(uw) if uw else None, max(uq) if uq else None])
    rec('largest relative errors of dL/ds and dm/ds at s = 0', [max(r['checks']['dL/ds error'] for r in res),
                                                               max(r['checks']['dm/ds error'] for r in res)])
    rec('total steps', sum(r['steps'] for r in res))
    out['seconds'] = time.time() - t0
    json.dump(out, open('wigbiane_poly.json', 'w'), indent=1)
    print('done', round(time.time() - t0, 1))
