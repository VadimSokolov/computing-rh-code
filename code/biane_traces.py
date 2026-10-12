# biane_traces.py: the Poisson trace of the zeros of zeta against the wave trace of the Liouville operator of Polya's
# function, and the complete monotonicity conditions of the Poisson trace, for Subsection sec:biane:cm of ch/strip.tex
# (Table tab:biane:traces, Figure fig:biane:traces, the first negative orders and the dominant ordinates) and Figure
# fig:biane:cm. Added for the book in October 2026, from the first author's note on Biane's paper on Polya
# (misc/reviews/nick-biane/, Sections 10.4 and 11 there), recomputed independently.
# Theta(v) = sum_{gamma>0} e^{-gamma v} is the Poisson trace K_{1/2}(v) of Section sm:sec:poisson, over the first NZ = 30000
# ordinates from Arb (acb.zeta_zeros, about 25750 at the top), plus the tail int_G^oo e^{-Tv} theta'(T)/pi dT - e^{-Gv} S(G)
# with the Riemann Siegel theta function theta, G midway between the last two ordinates and S(G) = N(G) - theta(G)/pi - 1
# exact; the neglected term v int_G^oo e^{-Tv} S(T) dT is estimated, later in the script, by Trudgian's bound for the integral of S
# (Math. Comp. 80 (2011), Theorem 2.2; not his 2014 bound for S(T), whose proof is incomplete).
# Theta_P(2v) = sum_n e^{-2 nu_n v} over the zeros t_n = 2 nu_n of P(t) = K_{it/2}(2 pi), from the exact phase condition
# Psi(nu_n) = n pi of Section sec:biane:zeros, Psi(nu) = Im log Gamma(1+i nu) - nu log pi - arg 0F1(;1+i nu;pi^2), solved
# by Newton's method in double precision for all t_n up to the same height (checked against the 40 digit zeros of
# biane_polya.py), plus the tail int_G^oo e^{-Tv} dx(T) + e^{-Gv}({x(G)} - 1/2), x(T) = Psi(T/2)/pi.
# D(v) = Theta(v) - Theta_P(2v) is compared with its limit 9/8 and with 9/8 - c v log(1/v), c = 2 pi - 1/(6 pi) - 1/(48 pi);
# the computation of the note (300 ordinates, a plain Weyl tail from the last zero) is repeated for comparison.
# Also: the first order n at which the term of an off line zero in (-1)^n Theta^(n)(v) turns negative, against
# (pi/2 + delta v)/arctan(delta/T); the ordinate that dominates the condition (n, v) among the first 300; and the data of
# the two figures. The code writes v for the variable y of the book. Run from any folder: python3 biane_traces.py; writes biane_traces.json there (about ten minutes).
import json
import math
import time

import mpmath as mp
import numpy as np
from flint import acb, ctx
from scipy.special import loggamma, psi as digamma

t0 = time.time()
out = {}


def rec(k, v):
    out[k] = v
    print(k, v, round(time.time() - t0, 1), flush=True)


mp.mp.dps = 30
EG = float(mp.euler)

# ---------------------------------------------------------------- ordinates of zeta
NZ = 30000
ctx.prec = 64
GAM = np.array([float(z.imag.mid().str(25, radius=False)) for z in acb.zeta_zeros(1, NZ)])
rec('zeta: number of ordinates, the last one', [NZ, GAM[-1]])
assert np.all(np.diff(GAM) > 0)
GZ = 0.5 * (GAM[-1] + GAM[-2])
GZm = GAM[:-1]                         # ordinates below GZ


def theta_rs(T):
    return float(mp.siegeltheta(T))


SG = (len(GZm)) - theta_rs(GZ) / math.pi - 1
rec('zeta: cut G, S(G)', [GZ, SG])

# ---------------------------------------------------------------- zeros of P from the phase, double precision
def F0(nu):
    nu = np.asarray(nu, dtype=float)
    b = 1 + 1j * nu
    term = np.ones_like(b)
    tot = np.ones_like(b)
    for k in range(1, 60):
        term = term * math.pi ** 2 / (k * (b + k - 1))
        tot = tot + term
    return tot


def dargF0(nu, h=1e-6):
    return (np.angle(F0(nu + h)) - np.angle(F0(nu - h))) / (2 * h)


def Psi(nu):
    nu = np.asarray(nu, dtype=float)
    return loggamma(1 + 1j * nu).imag - nu * math.log(math.pi) - np.angle(F0(nu))


def dPsi(nu):
    return digamma(1 + 1j * np.asarray(nu, dtype=float)).real - math.log(math.pi) - dargF0(nu)


# offset of the principal branch at the first zero (Psi(nu_1) must be pi)
nu1 = 9.768770083509977867
off = round((Psi(nu1) - math.pi) / (2 * math.pi))
rec('P: branch offset at nu_1 (multiples of 2 pi), Psi(nu_1)/pi', [off, float((Psi(nu1) - 2 * math.pi * off) / math.pi)])
nmax = int((Psi(GZ / 2 + 50) - 2 * math.pi * off) / math.pi) + 1
n = np.arange(1, nmax + 1, dtype=float)
# initial guess: invert nu log(nu/pi) - nu + pi/4 = n pi
nu = np.maximum(10.0, (n * math.pi) / np.maximum(1.0, np.log(n + 3.0)) + 5)
for _ in range(60):
    g = nu * np.log(nu / math.pi) - nu + math.pi / 4 - n * math.pi
    nu = nu - g / np.log(nu / math.pi)
    nu = np.maximum(nu, 9.0)
for _ in range(12):
    nu = nu - (Psi(nu) - 2 * math.pi * off - n * math.pi) / dPsi(nu)
resid = np.max(np.abs(Psi(nu) - 2 * math.pi * off - n * math.pi))
rec('P: zeros from the phase (to just above the cut), largest phase residual', [len(nu), float(resid)])
assert np.all(np.diff(nu) > 0)
TP = 2 * nu
rec('P: first five t_n = 2 nu_n', [float(x) for x in TP[:5]])
GP = GZ
TPm = TP[TP < GP]
xG = (Psi(GP / 2) - 2 * math.pi * off) / math.pi
rec('P: zeros below the cut, x(G), {x(G)}', [len(TPm), float(xG), float(xG - math.floor(xG))])


def dens_z(T):        # theta'(T)/pi
    return (mp.re(mp.digamma(mp.mpf(1) / 4 + 1j * T / 2)) - mp.log(mp.pi)) / (2 * mp.pi)


def dens_p(T):        # d/dT Psi(T/2)/pi = Psi'(T/2)/(2 pi), Psi'(nu) = log(nu/pi) - (pi^2 - 1/12)/nu^2 + O(nu^-4)
    nuv = T / 2
    return (mp.log(nuv / mp.pi) - (mp.pi ** 2 - mp.mpf(1) / 12) / nuv ** 2) / (2 * mp.pi)


def tail(dens, G, v, corr):
    v, G = mp.mpf(v), mp.mpf(float(G))
    I = mp.quad(lambda T: mp.exp(-T * v) * dens(T), [G, G + 1 / v, G + 5 / v, G + 30 / v, mp.inf])
    return float(I) + float(mp.exp(-G * v)) * float(corr)


def Theta(v):
    return float(np.sum(np.exp(-GZm * v))) + tail(dens_z, GZ, v, -SG)


def ThetaP2(v):
    return float(np.sum(np.exp(-TPm * v))) + tail(dens_p, GP, v, (xG - math.floor(xG)) - 0.5)


def main(v):
    return (math.log(1 / (2 * math.pi * v)) - EG) / (2 * math.pi * v)


c = 2 * math.pi - 1 / (6 * math.pi) - 1 / (48 * math.pi)
VS = [0.05, 0.02, 0.01, 0.005, 0.003, 0.002, 0.001, 5e-4, 2e-4, 1e-4, 1e-5]
rows = []
for v in VS:
    th, tp = Theta(v), ThetaP2(v)
    rows.append({'v': v, 'Theta': th, 'Theta_P(2v)': tp, 'main': main(v), 'Theta/main': th / main(v),
                 'Theta - main': th - main(v), 'Theta_P(2v) - main': tp - main(v), 'D': th - tp,
                 'D/log(1/v)': (th - tp) / math.log(1 / v), '9/8 - c v log(1/v)': 1.125 - c * v * math.log(1 / v),
                 'tail weight e^{-Gv}': math.exp(-GZ * v)})
    print(rows[-1], flush=True)
out['traces'] = rows
rec('c = 2 pi - 1/(6 pi) - 1/(48 pi)', c)
# bound on the neglected term: integrating by parts with int_G^T S, v int_G^oo e^{-Tv} S(T) dT = v^2 int_G^oo e^{-Tv}
# (int_G^T S(t) dt) dT, and Trudgian's bound |int_{t1}^{t2} S(t) dt| <= 2.067 + 0.059 log t2 for t2 > t1 > 168 pi, with
# Jensen's inequality for the concave log, bounds it by v e^{-Gv} (0.059 log(G + 1/v) + 2.067); the record is twice this
rec('bound 2 v e^{-Gv} (0.059 log(G + 1/v) + 2.067) on the neglected term, at v = 1e-3, 1e-4, 1e-5',
    [vv * math.exp(-GZ * vv) * 2 * (0.059 * math.log(GZ + 1 / vv) + 2.067) for vv in (1e-3, 1e-4, 1e-5)])

# the computation of the note: 300 ordinates, the zeros of P below gamma_300 + 20, plain Weyl tails from the last zeros
G300 = GAM[:300]
P312 = TP[TP < G300[-1] + 20]
rec('note: ordinates gamma_300, zeros of P below gamma_300 + 20, zeros of P below gamma_300',
    [float(G300[-1]), len(P312), int(np.sum(TP < G300[-1]))])


def weyl_tail(G, v):
    v = mp.mpf(v)
    return float(mp.quad(lambda T: mp.log(T / (2 * mp.pi)) / (2 * mp.pi) * mp.exp(-T * v), [G, G + 200 / v, mp.inf]))


nrows = []
for v in [0.003, 0.005, 0.01, 0.02]:
    a = float(np.sum(np.exp(-G300 * v))) + weyl_tail(float(G300[-1]), v)
    b = float(np.sum(np.exp(-P312 * v))) + weyl_tail(float(P312[-1]), v)
    nrows.append([v, a / main(v), (a - b) / math.log(1 / v), a - b])
rec('note style: v, Theta/main, D/log(1/v), D', nrows)

# ---------------------------------------------------------------- complete monotonicity conditions
def first_negative(T, d, v):
    phi = math.atan2(d, T)
    nn = 0
    while math.cos(nn * phi - d * v) > 0:
        nn += 1
    return nn, (math.pi / 2 + d * v) / phi


fn = []
for T, d, v in [(50, 0.25, 1.0), (50, 0.25, 0.3), (200, 0.1, 1.0), (14.13, 0.5, 0.5)]:
    nn, pred = first_negative(T, d, v)
    fn.append([T, d, v, nn, round(pred, 2)])
rec('first negative order: T, delta, v, first n with a negative term, (pi/2 + delta v)/arctan(delta/T)', fn)
dom = []
for nn, v in [(20, 1.0), (100, 1.0), (50, 0.5), (100, 5.0), (300, 10.0)]:
    w = nn * np.log(G300) - G300 * v
    k = int(np.argmax(w))
    srt = np.argsort(-w)[:2]
    dom.append([nn, v, nn / v, round(float(G300[k]), 4), [round(float(G300[j]), 4) for j in srt],
                [round(float(w[j]), 5) for j in srt]])
rec('dominant ordinate: n, v, n/v, argmax of n log gamma - gamma v over the first 300, top two, their weights', dom)
H = 3e12
rec('unconditional range H (pi - v) for H = 3e12 at v = 0.1, 1, 3', [H * (math.pi - v) for v in (0.1, 1.0, 3.0)])

# ---------------------------------------------------------------- figure data
vv = np.exp(np.linspace(math.log(1e-4), math.log(0.05), 60))
out['fig_traces'] = {'v': [float(x) for x in vv], 'Theta': [Theta(float(x)) for x in vv],
                     'ThetaP2': [ThetaP2(float(x)) for x in vv]}
TT = np.linspace(5, 200, 800)
prof = {}
for nn, v in [(100, 1.0), (300, 10.0)]:
    lw = nn * np.log(TT) - TT * v
    prof['%d,%g' % (nn, v)] = [float(x) for x in np.exp(lw - lw.max())]
out['fig_cm'] = {'T': [float(x) for x in TT], 'profiles': prof, 'zeros': [float(x) for x in G300[G300 < 200]],
                 'sign': {}}
for T, d, v in [(50, 0.25, 1.0), (14.13, 0.5, 0.5)]:
    phi = math.atan2(d, T)
    ns = np.arange(0, 800)
    out['fig_cm']['sign']['%g,%g,%g' % (T, d, v)] = [float(math.cos(k * phi - d * v)) for k in ns]
out['seconds'] = round(time.time() - t0, 1)
json.dump(out, open('biane_traces.json', 'w'), indent=1)
print('done', out['seconds'])
