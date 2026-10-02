# Reconstruction made in round 2 of the editing (misc/repro/reconstruct/p13_perturb_fig/), not the authors' code.
# It redraws fig/perturb.pdf, Figure fig:ch13:perturb of Chapter ch:finite (ch/ch13.tex), for which code/ has no
# script: code/plots3.py section 8 evaluates |W_F-W|/W in double precision for t in [1e-3, 10^-0.5], where the change
# underflows to 0, and plots it on a log axis from 1e-40 to 1, so its panel is empty. The archived figure (Matplotlib
# 3.10.8, 25 Sep 2026 14:06:33) plots log10 |W_F-W|/W at the 400 points t = logspace(-8, -2) on a linear axis from
# -3000 to 5. This script recomputes the curve in mpmath from the first 1700 zeros, with three versions of W:
#   w1700     W = sum of exp(-gamma^2 t) over the first 1700 zeros (what the caption says, and plots3.py uses)
#   two_term  the two term small time expansion of Chapter ch:heat, (log(1/t)-log(16 pi^2)-gamma_E)/(8 sqrt(pi t))+7/8
#   exact     the explicit formula eq:ch2:Warith, e^{t/4} + archimedean integral - prime sum, evaluated at 40 digits
# and writes
#   perturb.json        the pair, the three curves, and the numbers the caption and the text quote
#   perturb_w1700.pdf   the archived figure redrawn with W from the 1700 zeros
#   perturb_exact.pdf   the same figure with W from the explicit formula
# Run: python3 perturb.py [CORES]   with g_1_400.npy, g_401_1000.npy, g_1001_1700.npy (float64 ordinates from
# mpmath.zetazero, as written by misc/repro/zeros.py) beside it; they are recomputed if absent. About 3 minutes on
# 16 cores. python3 perturb.py --plot redraws the two figures from perturb.json. Drawn with Matplotlib 3.10.8, the
# version that drew the archived figure, perturb_w1700.pdf is byte for byte the archived fig/perturb.pdf apart from
# the creation date (compare_fig.py checks the page content).
import json, sys, time
import numpy as np, mpmath as mp
from multiprocessing import Pool

mp.mp.dps = 40
CORES = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 16
T0 = time.time()


def zero(n):
    mp.mp.dps = 30
    return float(mp.zetazero(n).imag)


def load_zeros():
    try:
        return np.concatenate([np.load(f) for f in ['g_1_400.npy', 'g_401_1000.npy', 'g_1001_1700.npy']]), 'npy'
    except FileNotFoundError:
        with Pool(CORES) as p:
            return np.array(p.map(zero, range(1, 1701), chunksize=10)), 'zetazero'


g, zsrc = load_zeros()
G2 = [mp.mpf(float(x)) ** 2 for x in g]
# the pair, exactly as in code/plots3.py section 8: the closest neighbours among the zeros in (1000, 2000)
gp = g[(g > 1000) & (g < 2000)]; k = int(np.argmin(np.diff(gp))); A_, B_ = gp[k], gp[k + 1]; gg = (A_ + B_) / 2
iA = int(np.nonzero(g == A_)[0][0]) + 1
A, B, G = mp.mpf(float(A_)), mp.mpf(float(B_)), mp.mpf(float(gg))
a = G ** 2 - mp.mpf(1) / 16 - 1j * G / 2          # a = -w^2 with w = 1/4 + i g
LOGPI = mp.log(mp.pi)
# prime powers n <= 1000 with Lambda(n) = log p (for t <= 0.01 the terms with n > 1000 are below e^{-1190})
PP = []
for p in [q for q in range(2, 1001) if all(q % d for d in range(2, int(q ** 0.5) + 1))]:
    n = p
    while n <= 1000:
        PP.append((n, mp.log(p))); n *= p


def dW(t):
    """W_F - W = 2 Re e^{-at} - e^{-A^2 t} - e^{-B^2 t}"""
    t = mp.mpf(t)
    return 2 * mp.re(mp.exp(-a * t)) - mp.exp(-A * A * t) - mp.exp(-B * B * t)


def W1700(t):
    t = mp.mpf(t)
    return mp.fsum(mp.exp(-x * t) for x in G2)


def W2(t):
    t = mp.mpf(t)
    return (mp.log(1 / t) - mp.log(16 * mp.pi ** 2) - mp.euler) / (8 * mp.sqrt(mp.pi * t)) + mp.mpf(7) / 8


def m(r):
    return mp.re(mp.digamma(mp.mpf(1) / 4 + 0.5j * r)) - LOGPI


def arch(t):
    """(1/4pi) int_R e^{-r^2 t} m(r) dr, split at r = 2^k up to 14/sqrt(t) (the rest is below e^{-196})"""
    t = mp.mpf(t); R = 14 / mp.sqrt(t)
    pts = [mp.mpf(0)] + [mp.mpf(2) ** j for j in range(-1, 80) if mp.mpf(2) ** j < R] + [R]
    v, e = mp.quad(lambda r: mp.exp(-r * r * t) * m(r), pts, error=True)
    return v / (2 * mp.pi), e / (2 * mp.pi)


def prime(t):
    t = mp.mpf(t)
    return mp.fsum(L / mp.sqrt(n) * mp.exp(-mp.log(n) ** 2 / (4 * t)) for n, L in PP) / (2 * mp.sqrt(mp.pi * t))


def Wexact(t):
    t = mp.mpf(t)
    ar, er = arch(t)
    return mp.exp(t / 4) + ar - prime(t), er


def point(t):
    """log10 |W_F-W| and log10 of the three versions of W at one t, as strings with 20 digits"""
    d = dW(t); w1 = W1700(t); w2 = W2(t); we, er = Wexact(t)
    s = lambda x: mp.nstr(x, 20)
    return dict(t=float(t), dW=s(d), log10_dW=s(mp.log10(abs(d))), W1700=s(w1), two_term=s(w2), exact=s(we),
                quad_err=float(er))


def rel(t, which):
    d = dW(t)
    w = {'w1700': W1700, 'two_term': W2, 'exact': lambda x: Wexact(x)[0]}[which](t)
    return mp.log10(abs(d) / w)


def golden(f, lo, hi, it=70):
    """maximum of f on [lo, hi] (f unimodal there), by golden section in log10 t"""
    r = (mp.sqrt(5) - 1) / 2
    x1 = hi - r * (hi - lo); x2 = lo + r * (hi - lo); f1 = f(x1); f2 = f(x2)
    for _ in range(it):
        if f1 > f2: hi, x2, f2 = x2, x1, f1; x1 = hi - r * (hi - lo); f1 = f(x1)
        else: lo, x1, f1 = x1, x2, f2; x2 = lo + r * (hi - lo); f2 = f(x2)
    x = (lo + hi) / 2
    return x, f(x)


def maxpoint(which):
    lt, v = golden(lambda u: rel(mp.mpf(10) ** u, which), mp.mpf(-6.6), mp.mpf(-5.6))
    return dict(t=float(mp.mpf(10) ** lt), rel=mp.nstr(mp.mpf(10) ** v, 8), log10_rel=float(v))


def draw(R):
    """the figure, in the style of the archived one (the rcParams of code/deconvplot.py, lplots.py and kmono.py)"""
    import matplotlib
    matplotlib.use('Agg'); import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': 0.25, 'axes.spines.top': False,
                         'axes.spines.right': False})
    tt = np.array(R['t'])
    for name, key in [('perturb_w1700', 'w1700'), ('perturb_exact', 'exact')]:
        fig, ax = plt.subplots(figsize=(6.5, 3.5))
        ax.semilogx(tt, R['log10_rel'][key], lw=1.6)
        ax.set_ylim(-3000, 5); ax.set_xlabel('t'); ax.set_ylabel(r'$\log_{10}|W_F-W|/W$')
        ax.set_title('Moving the pair near height %.0f off the line:\nrelative change in the heat trace'
                     % R['pair']['g'], fontsize=9)
        fig.tight_layout(); fig.savefig(name + '.pdf'); fig.savefig(name + '.png', dpi=200); plt.close()
    return matplotlib.__version__


if __name__ == '__main__':
    if '--plot' in sys.argv:
        print('matplotlib', draw(json.load(open('perturb.json')))); sys.exit()
    tt = np.logspace(-8, -2, 400)
    extra_t = [1e-8, 1e-7, 7e-7, 7.2e-7, 1e-6, 1e-5, 1e-4, 1e-3, 2e-3, 5e-3, 1e-2]
    with Pool(CORES) as p:
        P = p.map(point, list(tt), chunksize=4)
        X = p.map(point, extra_t)
        M = p.map(maxpoint, ['w1700', 'two_term', 'exact'])
    print('curve done %.0f s' % (time.time() - T0), flush=True)
    # checks of the explicit formula against data/arith.json (Wz there is the sum over the 1700 zeros)
    chk = {}
    for t, wz in [(0.002, '1.78338120537913'), (0.005, '0.537191964091041'), (0.01, '0.14969801125378')]:
        we, er = Wexact(t)
        chk[str(t)] = dict(exact=mp.nstr(we, 16), W1700=mp.nstr(W1700(t), 16), arith_json_Wz=wz, quad_err=float(er))
    # the zero of W_F - W near 1/(2 g^2), where the relative change vanishes (a downward spike of the curve that the
    # figure does not resolve), and the times where the curve crosses -1000 and -3000
    t0 = mp.findroot(dW, (mp.mpf('5e-8'), mp.mpf('1.5e-7')), solver='anderson')
    cross = {}
    for lev in [-1000, -3000]:
        u = mp.findroot(lambda u: rel(mp.mpf(10) ** u, 'w1700') - lev, mp.log10((-lev * mp.log(10)) / A ** 2))
        cross[str(lev)] = float(mp.mpf(10) ** u)
    L = lambda key, Q: [float(mp.log10(abs(mp.mpf(q['dW'])) / mp.mpf(q[key]))) for q in Q]
    y = {k: L(key, P) for k, key in [('w1700', 'W1700'), ('two_term', 'two_term'), ('exact', 'exact')]}
    yx = {k: L(key, X) for k, key in [('w1700', 'W1700'), ('two_term', 'two_term'), ('exact', 'exact')]}
    out = dict(
        zeros_source=zsrc, n_zeros=len(g), gamma_1700=float(g[-1]),
        pair=dict(index_A=iA, index_B=iA + 1, A=float(A_), B=float(B_), g=float(gg), gap=float(B_ - A_)),
        grid='numpy.logspace(-8, -2, 400)', t=tt.tolist(), log10_rel=y,
        points=P, extra_points=X, log10_rel_at_extra_points=dict(t=extra_t, **yx),
        maximum=dict(zip(['w1700', 'two_term', 'exact'], M)),
        grid_maximum={k: dict(t=float(tt[int(np.argmax(v))]), log10_rel=float(max(v))) for k, v in y.items()},
        zero_of_change=float(t0), crossings_w1700=cross, arith_json_check=chk,
        max_quad_err=max(q['quad_err'] for q in P), seconds=round(time.time() - T0, 1))
    out['matplotlib'] = draw(out)
    json.dump(out, open('perturb.json', 'w'), indent=1)
    print(json.dumps({k: out[k] for k in ['pair', 'maximum', 'grid_maximum', 'zero_of_change', 'crossings_w1700',
                                          'arith_json_check', 'max_quad_err', 'seconds', 'matplotlib']}, indent=1))
    for j, q in enumerate(X):
        print('t=%-8g log10|dW|=%-13.6f W: 1700 zeros %-13.8g two term %-13.8g exact %-13.8g   log10 rel: %s' % (
            q['t'], float(q['log10_dW']), float(q['W1700']), float(q['two_term']), float(q['exact']),
            ' '.join('%s %.6f' % (k, yx[k][j]) for k in yx)))
