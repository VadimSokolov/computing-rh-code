# Heat trace of the shifted clock from the primes (Section bp:sec:heat of the book).
# With g_t(x) = 1 - exp(-x^2/4t) and k_alpha(t) = pi^{-1} int_0^inf exp(-t theta^2) v_alpha(theta) d theta, for alpha >= 1
#   2 sqrt(pi t) k_alpha(t) = b_alpha + Abar_alpha(t) + int_[1,inf) g_t(log y) y^{-alpha} d(psi(y) - y),
#   Abar_alpha(t) = int_0^inf g_t(x) e^{-(alpha+2)x}/(1 - e^{-2x}) dx = sum_{k>=1} [1/c - sqrt(pi t) e^{c^2 t} erfc(c sqrt t)], c = alpha + 2k.
# Part 1 (floating point): the smallest alpha above which the primes up to N alone certify k_alpha > 0, for N = 1 (no
#   primes), 10^3, 10^5 and 10^7, and the same for the Thorin density v_alpha with no primes.
# Part 2 (ball arithmetic, python-flint): a lower bound for 2 sqrt(pi t) k_1(t), valid for every t > 0, from the primes up
#   to 10^7 and published bounds for psi(y) - y beyond: Buthe 2018 (|psi(y) - y| < 0.94 sqrt y, 11 < y <= 10^19),
#   Broadbent, Kadiri, Lumley, Ng and Wilk 2021 (1.93378e-8 y for y >= 10^19, 6.95e-9 y for y >= e^46), and Johnston
#   and Yang 2023 (8.86 y (log y)^1.514 exp(-0.8288 sqrt(log y)) for log y >= 3000).
# Part 3: the heat kernel form against the Thorin form at alpha = 2 and alpha = 1.
# Writes heat_primes.json. About 90 seconds on 32 cores: python3 heat_primes.py [workers]
import json, math, sys, time
from multiprocessing import Pool
import numpy as np
import mpmath as mp
from scipy.optimize import brentq, minimize_scalar
from scipy.special import erfcx, psi as digamma
from flint import arb, ctx

T0 = time.time()
N = 10**7
WORKERS = int(sys.argv[1]) if len(sys.argv) > 1 else 32
mp.mp.dps = 30

# prime powers up to N with Lambda(n) = log p
sieve = np.ones(N + 1, dtype=bool); sieve[:2] = False
for p in range(2, int(N**0.5) + 1):
    if sieve[p]:
        sieve[p * p::p] = False
primes = np.nonzero(sieve)[0].astype(np.int64)
pn, pp = [primes], [primes]
for p in primes[primes <= int(N**0.5)]:
    q = int(p) * int(p)
    while q <= N:
        pn.append(np.array([q])); pp.append(np.array([p])); q *= int(p)
pn = np.concatenate(pn); pp = np.concatenate(pp); o = np.argsort(pn); pn = pn[o]; pp = pp[o]
xf = np.log(pn.astype(float)); lf = np.log(pp.astype(float))
print(f'{len(pn)} prime powers up to {N}, {len(primes)} primes', flush=True)


def b_float(al):
    s = mp.mpf(al)
    if al == 1:
        return float(1 + mp.euler / 2 - mp.log(4 * mp.pi) / 2)
    return float(1 / s + 1 / (s - 1) - mp.log(mp.pi) / 2 + mp.digamma(s / 2) / 2 + mp.zeta(s, derivative=1) / mp.zeta(s))


# ---------------------------------------------------------------- Part 1: thresholds in floating point
U = np.linspace(-40.0, math.log(60.0), 6001); XU = np.exp(U); DU = U[1] - U[0]


def T_float(al, ts):  # trapezoid rule in u = log x
    h = np.exp(-(al + 2) * XU) / (-np.expm1(-2 * XU)) * XU
    g = -np.expm1(-XU[None, :] ** 2 / (4 * np.asarray(ts)[:, None]))
    f = g * h[None, :]
    return DU * (f.sum(axis=1) - 0.5 * (f[:, 0] + f[:, -1]))


_BINS = {}


def prime_bins(al, Nc):  # exact up to 2000, then bins of width 1e-4 in log n at the weighted mean
    if (al, Nc) in _BINS:
        return _BINS[(al, Nc)]
    if len(_BINS) > 8:
        _BINS.clear()
    _BINS[(al, Nc)] = r = _prime_bins(al, Nc)
    return r


def _prime_bins(al, Nc):
    m = pn <= Nc
    w = lf[m] * np.exp(-al * xf[m]); x = xf[m]; small = pn[m] <= 2000
    if not np.any(~small):
        return w, x
    key = np.floor(x[~small] / 1e-4).astype(np.int64)
    _, inv = np.unique(key, return_inverse=True)
    W = np.bincount(inv, weights=w[~small]); XM = np.bincount(inv, weights=w[~small] * x[~small]) / W
    return np.concatenate([w[small], W]), np.concatenate([x[small], XM])


def F_float(al, Nc, ts):  # b + T + P_N - Q^inf; for Nc < 2 no primes
    ts = np.atleast_1d(np.asarray(ts, float))
    out = b_float(al) + T_float(al, ts) - (1 / (al - 1) - np.sqrt(np.pi * ts) * erfcx((al - 1) * np.sqrt(ts)))
    if Nc >= 2:
        W, X = prime_bins(al, Nc)
        for i in range(0, len(ts), 50):
            tt = ts[i:i + 50]
            out[i:i + 50] += W @ (-np.expm1(-X[:, None] ** 2 / (4 * tt[None, :])))
    return out


TG = np.exp(np.linspace(math.log(1e-6), math.log(1e6), 481))


def minF(al, Nc, refine=True):
    v = F_float(al, Nc, TG); i = int(np.argmin(v))
    if not refine:
        return float(v[i]), float(TG[i])
    lo, hi = math.log(TG[max(i - 1, 0)]), math.log(TG[min(i + 1, len(TG) - 1)])
    r = minimize_scalar(lambda u: float(F_float(al, Nc, [math.exp(u)])[0]), bounds=(lo, hi), method='bounded',
                        options=dict(xatol=1e-6))
    return min(float(r.fun), float(v[i])), math.exp(float(r.x))


ladder = []
for Nc in [1, 10**3, 10**5, 10**7]:
    a = brentq(lambda al: minF(al, Nc)[0], 1.05, 5.0, xtol=1e-6)
    grid = np.linspace(a + 1e-3, 10.0, 40)
    above = min(minF(al, Nc, refine=False)[0] for al in grid)
    tmin = minF(a, Nc)[1]
    ladder.append(dict(N=Nc, alpha_N=a, argmin_t_at_threshold=tmin, min_over_alpha_grid_above=above))
    print(f'heat trace, primes up to {Nc}: threshold {a:.4f} (argmin t {tmin:.3g}); min over 40 alpha in (threshold, 10]: {above:.3g}', flush=True)
b10 = b_float(10.0)
print(f'alpha >= 10: b_10 = {b10:.5f} > 1/9 = {1/9:.5f}: {b10 > 1/9}', flush=True)


def thorin_noprime(al):  # min over theta of b + A - pole, no primes
    th = np.concatenate([np.linspace(1e-3, 50, 5000), np.linspace(50, 2000, 2000)])
    A = 0.5 * (digamma((al + 2 + 1j * th) / 2).real - digamma((al + 2) / 2))
    v = b_float(al) + A - th**2 / ((al - 1) * ((al - 1) ** 2 + th**2))
    i = int(np.argmin(v)); return float(v[i]), float(th[i])


a_thorin = brentq(lambda al: thorin_noprime(al)[0], 2.0, 6.0, xtol=1e-6)
thorin_above = min(thorin_noprime(al)[0] for al in np.linspace(a_thorin + 1e-3, 10.0, 40))
print(f'Thorin density, no primes: threshold {a_thorin:.4f} (argmin theta {thorin_noprime(a_thorin)[1]:.3f}); '
      f'min over 40 alpha in (threshold, 10]: {thorin_above:.3g}', flush=True)

# ---------------------------------------------------------------- Part 2: certificate at alpha = 1 in ball arithmetic
ctx.prec = 96
M = 10**4; H = 5e-4
PI = arb.pi(); X = arb(N).log(); b1 = 1 + arb.const_euler() / 2 - (4 * PI).log() / 2
bounds = [M]
while bounds[-1] < N:
    nxt = max(bounds[-1] + 1, int(M * math.exp(H * len(bounds))))
    bounds.append(min(nxt, N))
bounds = np.array(bounds, dtype=np.int64)
blk = np.searchsorted(bounds, pn, side='left') - 1           # block j holds bounds[j] < n <= bounds[j+1]
logp = {}
W_ex, X2_ex, Wb = [], [], [arb(0) for _ in range(len(bounds) - 1)]
psiN = arb(0)
for n, p, j in zip(pn.tolist(), pp.tolist(), blk.tolist()):
    lp = logp.get(p)
    if lp is None:
        lp = arb(p).log(); logp[p] = lp
    psiN += lp
    if n <= M:
        W_ex.append(lp / n); X2_ex.append(arb(n).log() ** 2)
    else:
        Wb[j] += lp / n
X2_b = [arb(int(b)).log() ** 2 for b in bounds[:-1]]         # g_t at the left end of each block: a lower bound
WALL = W_ex + Wb; X2ALL = X2_ex + X2_b
print(f'ball arithmetic: {len(W_ex)} exact terms, {len(Wb)} blocks; psi(N) = {psiN.str(12)}; {time.time() - T0:.0f} s', flush=True)

# tail: |int_(N,inf) g_t(log y) y^-1 d(psi(y) - y)| <= g_t(log N) |psi(N) - N|/N + (1 - g_t(log N)) sup_{x >= log N} eps(x)
#   + int_{log N}^inf eps(x) dx <= sup_{x >= log N} eps(x) + int_{log N}^inf eps(x) dx, since |psi(N) - N|/N <= eps(log N)
L19 = 19 * arb(10).log()
E_BKLNW1, E_BKLNW2 = arb(193378) / 10**13, arb(695) / 10**11
eps_sup = arb(94) / 100 * (-X / 2).exp()
assert eps_sup > E_BKLNW1
A_, B_, C_ = arb(886) / 100, arb(1514) / 1000, arb(8288) / 10000; u0 = arb(3000).sqrt(); Cp = C_ - (2 * B_ + 1) / u0
assert Cp > 0
J = 2 * A_ * u0 ** (2 * B_ + 1) * (-C_ * u0).exp() / Cp        # int_3000^inf A x^B e^{-C sqrt x} dx, with x = u^2
int_eps = (arb(188) / 100 * ((-X / 2).exp() - (-L19 / 2).exp()) + E_BKLNW1 * (46 - L19)
           + E_BKLNW2 * (3000 - 46) + J)
bnd = abs(psiN - N) / N
assert bnd < eps_sup
tau = eps_sup + int_eps
print(f'tau = {tau.str(6)}: sup eps {eps_sup.str(4)}, int eps {int_eps.str(4)} (Johnston-Yang part {J.str(3)}); boundary term {bnd.str(4)} <= sup eps', flush=True)

_P, _Q, _T = {}, {}, {}


def P_low(t):
    if t not in _P:
        c = -1 / (4 * arb(t)); s = arb(0)
        for w, x2 in zip(WALL, X2ALL):
            s += w * (1 - (c * x2).exp())
        _P[t] = s
    return _P[t]


def Q_up(t):
    if t not in _Q:
        a = arb(t); _Q[t] = X - (PI * a).sqrt() * (X / (2 * a.sqrt())).erf()
    return _Q[t]


def T_low(t, K=1000):  # partial sum of positive terms; terms with k > K or c sqrt t > 25 are dropped
    if t not in _T:
        a = arb(t); st = a.sqrt(); spt = (PI * a).sqrt(); s = arb(0)
        for k in range(1, K + 1):
            z = (1 + 2 * k) * st
            if z > 25:
                break
            s += arb(1) / (1 + 2 * k) - spt * (z * z).exp() * z.erfc()
        _T[t] = s
    return _T[t]


TARGET = arb(21) / 1000


def certify(cell):
    t0, t1 = cell; work = [(t0, t1)]; worst = None; wt = None; nc = 0
    while work:
        a, b = work.pop()
        lb = b1 + T_low(b) + P_low(b) - Q_up(a) - tau
        if lb > TARGET:
            nc += 1; lo = float(lb.lower())
            if worst is None or lo < worst:
                worst, wt = lo, (a, b)
            continue
        if b / a < 1 + 1e-7:
            return dict(ok=False, cell=(a, b), lb=lb.str(6), cells=nc)
        r = (b / a) ** 0.125; pts = [a * r**i for i in range(8)] + [b]
        work.extend((pts[i], pts[i + 1]) for i in range(8))
    return dict(ok=True, worst=worst, where=wt, cells=nc)


TMIN, TMAX = 1e-10, 1e8
edges = np.exp(np.linspace(math.log(TMIN), math.log(TMAX), 4 * WORKERS * 8 + 1))
edges[0], edges[-1] = TMIN, TMAX                                # pin the outer edges to the tail cutoffs, so the cells
cells = [(float(edges[i]), float(edges[i + 1])) for i in range(len(edges) - 1)]   # and the two tail bounds share endpoints
assert cells[0][0] == TMIN and cells[-1][1] == TMAX and all(cells[i][1] == cells[i + 1][0] for i in range(len(cells) - 1))
with Pool(WORKERS) as pool:
    res = pool.map(certify, cells, chunksize=1)
lb_small = b1 + T_low(TMIN) + P_low(TMIN) - X - tau             # t <= TMIN: T and P decrease in t, Q <= log N
lb_large = b1 - Q_up(TMAX) - tau                                # t >= TMAX: T, P >= 0 and Q decreases in t
ok = all(r['ok'] for r in res) and lb_small > TARGET and lb_large > TARGET
good = [r for r in res if r['ok']]
w = min(good, key=lambda r: r['worst']) if good else None
cert = dict(target=0.021, all_ok=bool(ok), cells=int(sum(r['cells'] for r in res)), t_range=[TMIN, TMAX],
            min_lower_bound=min(([w['worst']] if w else []) + [float(lb_small.lower()), float(lb_large.lower())]),
            min_cell=w['where'] if w else None, lb_small_t=lb_small.str(6), lb_large_t=lb_large.str(6),
            failures=[r for r in res if not r['ok']][:5],
            b1=b1.str(15), tau=tau.str(8), tau_boundary=bnd.str(6), tau_sup_eps=eps_sup.str(6), tau_int_eps=int_eps.str(6),
            tau_johnston_yang=J.str(4), psiN=psiN.str(15), exact_terms=len(W_ex), blocks=len(Wb), block_width=H,
            seconds=time.time() - T0)
print('certificate at alpha = 1:', {k: cert[k] for k in ['all_ok', 'cells', 'min_lower_bound', 'min_cell']}, flush=True)

# float profile of the certified lower bound (without cell losses), to locate its minimum
tp = np.exp(np.linspace(math.log(1e-3), math.log(1e6), 400))
W1, X1 = prime_bins(1.0, N)
prof = []
for t in tp:
    P = float(W1 @ -np.expm1(-X1**2 / (4 * t)))
    Qn = float(X.mid()) - math.sqrt(math.pi * t) * math.erf(float(X.mid()) / (2 * math.sqrt(t)))
    prof.append(b_float(1) + float(T_float(1.0, [t])[0]) + P - Qn - float(tau.mid()))
i = int(np.argmin(prof))
cert['profile_min'] = dict(t=float(tp[i]), value=float(prof[i]))
print(f'float profile of b1 + T + P - Q - tau: min {prof[i]:.6f} at t = {tp[i]:.4g}', flush=True)

# ---------------------------------------------------------------- Part 3: heat kernel form against the Thorin form


def v_mp(al, th):
    s = mp.mpc(al, th)
    return mp.re(1 / s + 1 / (s - 1) - mp.log(mp.pi) / 2 + mp.digamma(s / 2) / 2 + mp.zeta(s, derivative=1) / mp.zeta(s))


def thorin_side(al, t):  # 2 sqrt(pi t) k_alpha(t) from the Thorin form; [0, 1e-6] by the value b_alpha at 0
    top = math.sqrt(80 / t); e = mp.mpf('1e-6')
    I = mp.quad(lambda th: mp.exp(-t * th * th) * v_mp(al, th), mp.linspace(e, top, 12)) + e * b_float(al)
    return float(2 * mp.sqrt(mp.pi * t) * I / mp.pi)


checks = []
for al, t in [(2.0, 0.5), (2.0, 3.0), (1.0, 1.0), (1.0, 12.0), (1.0, 100.0)]:
    W_, X_ = prime_bins(al, N)
    P = float(W_ @ -np.expm1(-X_**2 / (4 * t)))
    if al == 1.0:
        Qv = float(X.mid()) - math.sqrt(math.pi * t) * math.erf(float(X.mid()) / (2 * math.sqrt(t)))
    else:
        Qv = 1 / (al - 1) - math.sqrt(math.pi * t) * erfcx((al - 1) * math.sqrt(t))
    heat = b_float(al) + float(T_float(al, [t])[0]) + P - Qv
    th = thorin_side(al, t)
    checks.append(dict(alpha=al, t=t, heat_kernel_form=heat, thorin_form=th, difference=heat - th))
    print(f'alpha {al} t {t}: heat kernel form {heat:.10f}, Thorin form {th:.10f}, difference {heat - th:.2e}', flush=True)

out = dict(N=N, ladder=ladder, b10=b10, thorin_noprime_threshold=a_thorin, thorin_noprime_min_above=thorin_above,
           certificate=cert, checks=checks,
           seconds=time.time() - T0)
json.dump(out, open('heat_primes.json', 'w'), indent=1)
print(f'done in {time.time() - T0:.0f} s', flush=True)
