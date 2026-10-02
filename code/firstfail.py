# code/firstfail.py (September 2026, from misc/reviews/nick-audit/A-first-failure): Table tab:ch13:first, the first order
# at which a derivative sign of the perturbed heat trace fails. For the closest pair A < B of consecutive zeros with both
# ordinates in each of six height windows, moved to the quadruple 1/2 +- delta +- i g with g = (A+B)/2 as in Theorem
# thm:ch13:perturb (whose pair is the one near 1977.22), it finds the least integer k such that
# c_k^F(t) = (-1)^k W_F^(k)(t) < 0 for some t > 0, for delta = 1/4 and, at the pair near 1977.22, also delta = 1/10.
# Book convention: one zero of each pair {rho, 1-rho}, a = -(rho-1/2)^2, c_k(t) = sum a^k e^{-at}. With t = k/u^2 the
# normalised value S(k,u) = c_k^F(t) u^{-2k} e^k is
#   S = sum_j exp(k psi(r_j)) + 2 exp(k psi(r_+) + 2 t delta^2) cos(2k(delta g/u^2 - arctan(delta/g))),
# psi(r) = log(1+r) - r, r_j = (gamma_j^2 - u^2)/u^2, r_+ = (g^2 + delta^2 - u^2)/u^2, over the line zeros within 40 of g
# other than A and B. The omitted zeros below 3e12 lie on the line [PT21] and add positive terms; their size, and a bound
# for the zeros above g + 40 taken on the line, are recorded (below 1e-1600 relative to S); firstfail_arb.py also
# bounds any zeros off the line above 3e12. Scan: u on [g-8, g+8] in
# steps of 0.002, the local minima refined by bounded Brent minimisation; k over every integer up to 500 and then a
# geometric grid of ratio 1.002 until the first negative minimum; then every integer of the last grid interval. For
# |u - g| > 8 a domination inequality, checked here for the line zeros within 8 of g, shows that each of at least three
# zeros on either side outweighs half the modulus of the quadruple term at every order, so c_k^F > 0 there.
# The scan is computed, not certified; firstfail_arb.py certifies the negative value at the first failing order.
# Also the Gaussian smoothing G(u,t) = sum exp(-2t(gamma~ - u)^2) of Proposition prop:ch12:ray for the same zeros: the
# least t at which its minimum over u is negative, and the least real k at which G(u, k/u^2) < 0 for some u; and the
# first failure time of a model with the line zeros at their mean density (density_model), at these heights and at 3e12.
# Reads zeros_1700.txt; writes firstfail.json. Seven cases in parallel, about five minutes on eight cores.
import json
import os
import time
from multiprocessing import Pool

os.environ.setdefault('OMP_NUM_THREADS', '1')

import mpmath as mp
import numpy as np
from scipy.optimize import brentq, minimize_scalar

STR = open('zeros_1700.txt').read().split()
G = np.array([float(s) for s in STR])
WINDOWS = [(100, 200), (200, 400), (400, 700), (700, 1000), (1000, 1500), (1500, 2000)]
CASES = [(w, '1/4') for w in WINDOWS] + [((1500, 2000), '1/10')]
H, UW, ZW = 0.002, 8.0, 40.0
THR, NMIN, RATIO = 0.05, 5, 1.002   # refine every grid minimum below THR and the NMIN smallest


def pair_index(lo, hi):   # closest pair of consecutive zeros with both ordinates in [lo, hi]
    d = np.diff(G)
    return int(np.argmin(np.where((G[:-1] >= lo) & (G[1:] <= hi), d, 1e9)))


def N_upper(T):   # Trudgian: N(T) <= (T/2pi) log(T/(2 pi e)) + 7/8 + 0.112 log T + 0.278 log log T + 2.510 + 0.2/T
    return (T / (2 * np.pi) * np.log(T / (2 * np.pi * np.e)) + 0.875 + 0.112 * np.log(T) + 0.278 * np.log(np.log(T))
            + 2.51 + 0.2 / T)


class Case:
    def __init__(self, lo, hi, dstr):
        i = pair_index(lo, hi)
        num, den = dstr.split('/')
        self.i, self.dstr, self.delta = i, dstr, int(num) / int(den)
        self.A, self.B, self.g = G[i], G[i + 1], 0.5 * (G[i] + G[i + 1])
        self.line = np.delete(G, [i, i + 1])
        self.Z = self.line[np.abs(self.line - self.g) <= ZW]
        self.U = self.g + np.arange(-int(round(UW / H)), int(round(UW / H)) + 1) * H
        U, Z, g, d = self.U, self.Z, self.g, self.delta
        R = (Z[None, :] - U[:, None]) * (Z[None, :] + U[:, None]) / U[:, None] ** 2
        self.Psi = np.log1p(R) - R
        rp = ((g - U) * (g + U) + d * d) / U ** 2
        self.pq = np.log1p(rp) - rp + 2 * d * d / U ** 2   # log modulus of the quadruple term, per unit k
        self.cq = 2 * (d * g / U ** 2 - np.arctan(d / g))  # its phase, per unit k
        self.D2 = (Z[None, :] - U[:, None]) ** 2

    def grid(self, k):
        return np.exp(k * self.Psi).sum(axis=1) + 2 * np.exp(k * self.pq) * np.cos(k * self.cq)

    def point(self, k, u):
        g, d, Z = self.g, self.delta, self.Z
        r = (Z - u) * (Z + u) / u ** 2
        rp = ((g - u) * (g + u) + d * d) / u ** 2
        return (np.exp(k * (np.log1p(r) - r)).sum()
                + 2 * np.exp(k * (np.log1p(rp) - rp + 2 * d * d / u ** 2)) * np.cos(2 * k * (d * g / u ** 2 - np.arctan(d / g))))

    def refine(self, v, f):   # minimum over the u grid of the values v, local minima refined for the function f
        U = self.U
        loc = np.nonzero((v[1:-1] <= v[:-2]) & (v[1:-1] <= v[2:]))[0] + 1
        cand = set(loc[np.argsort(v[loc])][:NMIN].tolist()) | set(loc[v[loc] < THR].tolist())
        j0 = int(np.argmin(v))
        best = (float(v[j0]), float(U[j0]))
        for j in cand:
            r = minimize_scalar(f, bounds=(U[j] - H, U[j] + H), method='bounded', options=dict(xatol=1e-10))
            if r.fun < best[0]:
                best = (float(r.fun), float(r.x))
        return best

    def mink(self, k):   # minimum of S(k, u) over u in [g-8, g+8]
        return self.refine(self.grid(k), lambda u: self.point(k, u))

    def gauss(self, u, t):   # the Gaussian smoothing G(u, t) of the same zeros
        g, d = self.g, self.delta
        return np.exp(-2 * t * (self.Z - u) ** 2).sum() + 2 * np.exp(-2 * t * ((g - u) ** 2 - d * d)) * np.cos(4 * t * d * (g - u))

    def gmin(self, t=None, k=None):   # minimum over u of G(u, t) at fixed t, or of G(u, k/u^2) at fixed k
        g, d, U = self.g, self.delta, self.U
        T = t if k is None else k / U ** 2
        Tc = T if k is None else T[:, None]
        v = np.exp(-2 * Tc * self.D2).sum(axis=1) + 2 * np.exp(-2 * T * ((g - U) ** 2 - d * d)) * np.cos(4 * T * d * (g - U))
        return self.refine(v, (lambda u: self.gauss(u, t)) if k is None else (lambda u: self.gauss(u, k / u ** 2)))


def dominated(C):
    """For u >= g+8 (u <= g-8) a line zero gamma in (g, g+8] ([g-8, g)) has weight at least half the modulus of the
    quadruple term at every order k >= 0 when P >= 0 (Q >= 0). With a2 = g^2+delta^2 the condition is
    u^2 log(gamma^2/a2) + a2 - gamma^2 >= 2 delta^2, whose left side increases with |u - g| on the side of the zero, so
    it suffices at |u - g| = 8. Above g, log y >= 2(y - 1)/(y + 1) for y = gamma^2/a2 >= 1 and gamma^2 + a2 <= 2u^2 give
    the left side >= P + 2 delta^2. Below g, with e = 1 - gamma^2/a2, the left side at u = g-8 is
    Q + 2 delta^2 + a2 e^2/2 + (g-8)^2 (e + log(1-e)), and the last two terms vanish at e = 0 and increase with e,
    their derivative being e (a2 - (g-8)^2/(1-e)) >= 0 for gamma >= g-8. Three such zeros on each side give c_k^F > 0
    for |u - g| >= 8."""
    g, d = C.g, C.delta
    a2 = g * g + d * d
    up = C.line[(C.line > g) & (C.line <= g + UW)]
    dn = C.line[(C.line < g) & (C.line >= g - UW)]
    P = (up ** 2 - a2) * (2 - (a2 + up ** 2) / (g + UW) ** 2) / 2 - 2 * d * d
    Q = (a2 - dn ** 2) * (a2 + dn ** 2 - 2 * (g - UW) ** 2) / (2 * a2) - 2 * d * d
    return dict(n_above=int(len(up)), n_above_dominating=int((P >= 0).sum()), n_below=int(len(dn)),
                n_below_dominating=int((Q >= 0).sum()), holds=bool((P >= 0).sum() >= 3 and (Q >= 0).sum() >= 3))


def omitted(C, k, u):
    """log10 of the omitted terms at (k, u), relative to u^{2k} e^{-k}: the listed zeros farther than 40 from g, all
    zeros below g-40 (at most N(g-40), each below the weight at (g-40)^2 < k/t), all zeros above g+40 (the weight is
    below its tangent bound there; with N(s) <= s^2 their sum is at most the weight at (g+40)^2 times (g+40)^2 + 1/lam)"""
    t = k / u ** 2
    far = C.line[np.abs(C.line - C.g) > ZW]
    r = (far - u) * (far + u) / u ** 2
    lw = k * (np.log1p(r) - r)
    m = lw.max()
    lo, hi = C.g - ZW, C.g + ZW
    rl, r0 = (lo ** 2 - u ** 2) / u ** 2, (hi ** 2 - u ** 2) / u ** 2
    lam = t - k / hi ** 2
    return dict(listed_far_log10=float((m + np.log(np.exp(lw - m).sum())) / np.log(10)),
                below_bound_log10=float((np.log(N_upper(lo)) + k * (np.log1p(rl) - rl)) / np.log(10)),
                above_bound_log10=float((k * (np.log1p(r0) - r0) + np.log(hi ** 2 + 1 / lam)) / np.log(10)))


def S_mp(C, k, u, dps=40):   # S(k, u) again in mpmath, from the 70 digit ordinates
    mp.mp.dps = dps
    g = (mp.mpf(STR[C.i]) + mp.mpf(STR[C.i + 1])) / 2
    num, den = C.dstr.split('/')
    d = mp.mpf(num) / mp.mpf(den)
    Zs = [mp.mpf(s) for s, x in zip(STR, G) if abs(x - C.g) <= ZW and x != C.A and x != C.B]
    u = mp.mpf(u)
    u2, t = u * u, k / (u * u)
    s = mp.fsum(mp.exp(k * mp.log(z * z / u2) - (z * z - u2) * t) for z in Zs)
    a = mp.mpc(g * g - d * d, -2 * d * g)
    return s + 2 * mp.re(mp.exp(k * mp.log(a / u2) - (a - u2) * t))


def density_model(g, delta):
    """The first failure time of a model in which the line zeros around an isolated quadruple are spread at the mean
    density D = log(g/2pi)/(2pi): the least t at which the quadruple's term 2 e^{x - q^2/(4x)} cos q in G, with
    x = 2 t delta^2 and q = 4 t delta |g - u|, is below minus the background D sqrt(pi/(2t)), for the phase q = pi and for
    the best phase, which solves tan q = -q/(2x) in (pi/2, pi)."""
    D = np.log(g / (2 * np.pi)) / (2 * np.pi)

    def excess(t, best):
        x = 2 * t * delta ** 2
        q = brentq(lambda q: np.tan(q) + q / (2 * x), np.pi / 2 + 1e-12, np.pi) if best else np.pi
        return 2 * np.exp(x - q * q / (4 * x)) * (-np.cos(q)) - D * np.sqrt(np.pi / (2 * t))

    out = {}
    for name, best in (('phase_pi', False), ('best_phase', True)):
        t = 0.05
        while excess(t, best) < 0:
            t += 0.01
        out[name] = brentq(lambda s: excess(s, best), t - 0.01, t, xtol=1e-12)
    return out


def first(fun, x0, step, tol, geometric=False):   # least x with fun(x) < 0, by stepping and bisection
    x = x0
    while fun(x)[0] >= 0:
        x = x * step if geometric else x + step
    lo, hi = (x / step, x) if geometric else (x - step, x)
    while (hi / lo - 1 if geometric else hi - lo) > tol:
        mid = np.sqrt(lo * hi) if geometric else (lo + hi) / 2
        if fun(mid)[0] < 0:
            hi = mid
        else:
            lo = mid
    return hi, fun(hi)


def run(case):
    (lo, hi), dstr = case
    t0 = time.time()
    C = Case(lo, hi, dstr)
    ks, ms, k = [], [], 1
    while True:   # bracket the first failure
        m, u = C.mink(k)
        ks.append(k); ms.append(m)
        if m < 0:
            break
        k = k + 1 if k < 500 else int(np.floor(k * RATIO))
    kfail, kprev = ks[-1], ks[-2]
    ints = list(range(kprev + 1, kfail + 1))   # every integer of the last grid interval
    res = [C.mink(kk) for kk in ints]
    mi = np.array([r[0] for r in res])
    j = int(np.nonzero(mi < 0)[0][0])
    kstar, (m_star, u_star) = ints[j], res[j]
    m_prev, u_prev = res[j - 1] if j > 0 else C.mink(kprev)
    # the Gaussian smoothing: below t_G every G(u, t) >= 0, so the k scan can start at t_G (g-8)^2
    tG, (vG, uG) = first(lambda t: C.gmin(t=t), 0.01, 1e-3, 1e-8)
    kG, (vk, uk) = first(lambda k: C.gmin(k=k), tG * (C.g - UW) ** 2, 1.0002, 1e-9, geometric=True)
    out = dict(window=[lo, hi], delta=dstr, zero_numbers=[C.i + 1, C.i + 2], A=STR[C.i], B=STR[C.i + 1], g=C.g,
               gap=C.B - C.A, n_zeros_within_40=int(len(C.Z)), bracket=[int(kprev), int(kfail)],
               kstar=int(kstar), t_star=kstar / u_star ** 2, u_star=u_star, u_star_minus_g=u_star - C.g,
               min_S_at_kstar=m_star, min_S_at_kstar_minus_1=m_prev, t_at_kstar_minus_1=(kstar - 1) / u_prev ** 2,
               negative_after_kstar=int((mi[j:] < 0).sum()), integers_after_kstar=len(ints) - j,
               kstar_over_g2=kstar / C.g ** 2, kstar_over_glog2=kstar / (C.g * np.log(C.g / (2 * np.pi))) ** 2,
               mp_S_at_kstar=mp.nstr(S_mp(C, kstar, u_star), 12), mp_S_at_kstar_minus_1=mp.nstr(S_mp(C, kstar - 1, u_prev), 12),
               omitted_at_kstar=omitted(C, kstar, u_star), outside_window=dominated(C),
               gauss_t=tG, gauss_u_minus_g=uG - C.g, gauss_k=kG, gauss_t_on_ray=kG / uk ** 2,
               gauss_k_rel_diff=(kG - kstar) / kstar, density_model_t=density_model(C.g, C.delta), seconds=time.time() - t0)
    print(json.dumps(out), flush=True)
    return out


if __name__ == '__main__':
    with Pool(len(CASES)) as pool:
        rows = pool.map(run, CASES)
    model = [dict(g=g, delta=d, **density_model(g, d)) for g in (185.24, 1977.22, 3e12) for d in (0.25, 0.1)]
    json.dump(dict(cases=rows, density_model=model), open('firstfail.json', 'w'), indent=1)
    for m in model:
        print('mean density model at height %.6g, delta %.2f: t* = %.3f with phase pi, %.3f with the best phase'
              % (m['g'], m['delta'], m['phase_pi'], m['best_phase']))
    for r in rows:
        print('%8.2f %4s k* = %9d  t* = %.6f  k*/g^2 = %.3f  min S(k*) = %.3e  min S(k*-1) = %.3e  Gauss k = %.1f (%+.2e)'
              % (r['g'], r['delta'], r['kstar'], r['t_star'], r['kstar_over_g2'], r['min_S_at_kstar'],
                 r['min_S_at_kstar_minus_1'], r['gauss_k'], r['gauss_k_rel_diff']))
