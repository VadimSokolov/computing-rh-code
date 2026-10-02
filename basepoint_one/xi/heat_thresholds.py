# Table 9.7 (bp:tab:heat), rows one to five: certificates in ball arithmetic (python-flint arb) for the thresholds that
# heat_primes.py locates in floating point.
# Heat trace rows. With g_t(x) = 1 - exp(-x^2/4t), the right side of (bp:eq:heatlb) is
#   F_N(alpha, t) = b_alpha + Abar_alpha(t) + sum_{n<=N} Lambda(n) n^-alpha g_t(log n) - L_alpha(t),
#   Abar_alpha(t) = sum_{k>=1} [1/c - sqrt(pi t) e^{c^2 t} erfc(c sqrt t)], c = alpha + 2k (every term positive),
#   L_alpha(t) = int_0^inf g_t(x) e^{-(alpha-1)x} dx = 1/(alpha-1) - sqrt(pi t) e^{(alpha-1)^2 t} erfc((alpha-1) sqrt t),
#   b_alpha = 1/alpha + 1/(alpha-1) - log(pi)/2 + psi(alpha/2)/2 + zeta'(alpha)/zeta(alpha).
# Abar and the prime sum decrease in alpha and in t, L decreases in alpha and in t, so on a cell [a1, a2] x [t1, t2]
#   F >= b([a1, a2]) + Abar_low(a2, t2) + P_low(a2, t2) - L(a1, t1),
# with b bounded below on [a1, a2] by its mean value form, Abar_low the terms k <= 40 and a lower bound for the rest,
# and P_low the prime powers up to 10^4 exactly and the rest in blocks l0 < log n <= l1 of ratio e^{1/250}, on which
# g_t(log n) >= g_t(l0) + (log n - l0) min(g_t'(l0), g_t'(l1)), since g_t'(x) = x e^{-x^2/4t}/2t is unimodal; the block
# sums of Lambda(n) n^-a2 and of Lambda(n) n^-a2 (log n - l0) are computed exactly. The cells cover
# [alpha_hi, 10] x [0, infinity]: t = 0 enters through
# L(a1, 0) = 1/(a1 - 1), and t = infinity through Abar = P = 0. Where this corner bound fails on a finite time cell
# [t1, t2], the mean value form in t is tried before the cell is split:
#   F >= b([a1, a2]) + Abar_low(a2, t2) + P_low(a2, tm) - L(a1, tm) - r sup_{[t1, t2]} |d/dt (P(a2, t) - L(a1, t))|,
# tm the midpoint and r the half width, with d/dt g_t(x) = -(x^2/4t^2) e^{-x^2/4t} enclosed for t in [t1, t2] (and, in a
# block, log n in [l0, l1]) and dL/dt = -E'(z)/(2 sqrt t), z = (a1 - 1) sqrt t,
# E'(z) = sqrt(pi) e^{z^2} erfc(z)(1 + 2z^2) - 2z. For alpha >= 10 the book's argument applies (b_alpha increases,
# b_10 > 1/9 >= the pole term), and b_10 > 1/9 is checked here in ball arithmetic.
# Below: at alpha_lo the right side is negative at one time t*, certified with every prime power up to N, the terms of
# Abar up to k = 200 and the bound sum_{k>200} <= 1/(8 t (alpha + 401)^2) for the rest.
# Thorin row. G(alpha, theta) = b_alpha + A_alpha(theta) - theta^2/((alpha-1)((alpha-1)^2 + theta^2)),
# A_alpha(theta) = (Re psi((alpha+2+i theta)/2) - psi((alpha+2)/2))/2, increasing in theta; the pole term decreases in
# alpha and increases in theta. Cells cover [alpha_hi, 10] x [0, 2000]; beyond 2000, A >= A(alpha, 2000) and the pole
# term is at most 1/(alpha-1).
# The printed thresholds are certified to four decimals: alpha_lo = printed - 5e-5 fails and every alpha >= alpha_hi
# succeeds, alpha_hi being printed + 5e-5, or the three decimal bound of the table when that is smaller (1.284 for
# N = 10^7), less 1e-12 against the rounding of float end points; so the threshold lies in (alpha_lo, alpha_hi], and the
# bounds of the table are proved. The alpha cells have widths FAC times the distance to the estimated threshold: the
# corner bounds lose about the width times |dP/dalpha| + |dL/dalpha|, each near 1/(alpha - 1)^2, while F itself grows
# with slope 0.18 to 0.76, so FAC must shrink as the threshold approaches 1.
# Usage: python3 heat_thresholds.py [workers]      writes heat_thresholds.json
import json, math, sys, time
from multiprocessing import Pool
import numpy as np
from flint import arb, acb, acb_series, ctx

T0 = time.time()
ctx.prec = 128
WORKERS = int(sys.argv[1]) if len(sys.argv) > 1 else 32
NMAX = 10**7
M = 10**4                    # exact prime powers up to M
H = 4e-3                     # block ratio e^H above M
FAC = {1: 0.04, 10**3: 0.025, 10**5: 0.015, 10**7: 0.01, 'thorin': 0.04}
PRINTED = {1: '3.1900', 10**3: '1.5320', 10**5: '1.3667', 10**7: '1.2840'}
ESTIMATE = {1: 3.190001491035211, 10**3: 1.5320409692773362, 10**5: 1.3666608788793266, 10**7: 1.2839912918396332}
THORIN_PRINTED, THORIN_ESTIMATE = '3.5748', 3.574776332980113
TMIN, TMAX = 1e-8, 1e8
PI = arb.pi()

# prime powers up to NMAX with Lambda(n) = log p
sieve = np.ones(NMAX + 1, dtype=bool); sieve[:2] = False
for p in range(2, int(NMAX**0.5) + 1):
    if sieve[p]:
        sieve[p*p::p] = False
primes = np.nonzero(sieve)[0].astype(np.int64)
pn, pp = [primes], [primes]
for p in primes[primes <= int(NMAX**0.5)]:
    q = int(p)*int(p)
    while q <= NMAX:
        pn.append(np.array([q])); pp.append(np.array([p])); q *= int(p)
pn = np.concatenate(pn); pp = np.concatenate(pp); o = np.argsort(pn); pn = pn[o]; pp = pp[o]
print(f'{len(pn)} prime powers up to {NMAX}', flush=True)
bounds = [M]
while bounds[-1] < NMAX:
    bounds.append(min(max(bounds[-1] + 1, int(M*math.exp(H*len(bounds)))), NMAX))
bounds = np.array(bounds, dtype=np.int64)
LOGP = {}
for p in np.unique(pp).tolist():
    LOGP[p] = arb(p).log()


def prime_data(N):
    """Exact terms (Lambda, log n) for n <= min(N, M); blocks (indices of prime powers, j) for the nonempty blocks
    bounds[j] < n <= bounds[j+1]; ends[j] = log bounds[j]."""
    m = pn <= N
    n_, p_ = pn[m], pp[m]
    ex = [(LOGP[p], LOGP[p]*round(math.log(n)/math.log(p))) for n, p in zip(n_.tolist(), p_.tolist()) if n <= M]
    blocks = []; ends = []
    if N > M:
        cut = np.searchsorted(n_, bounds, side='right')
        for j in range(len(bounds) - 1):
            if cut[j + 1] > cut[j]:
                blocks.append((np.arange(cut[j], cut[j + 1]), j))
        ends = [arb(int(b)).log() for b in bounds[:blocks[-1][1] + 2].tolist()] if blocks else []
    return n_, p_, ex, blocks, ends


_W = {}
_LN = {}


def logn(N):
    """log n for the prime powers n <= N, as balls (cached per process)."""
    if N not in _LN:
        _LN.clear()
        _LN[N] = [arb(int(n)).log() for n in PD[N][0].tolist()]
    return _LN[N]


def coef(N, a):
    """Lambda(n) n^-a for the exact terms; for the boundaries k of the blocks, A[k] = sum Lambda(n) n^-a over the block
    starting at k and B[k] = log(bounds[k]) sum Lambda(n) n^-a (log n - log bounds[k-1]) over the block ending at k
    (zero for empty blocks); and the sum of all A. Exactly in ball arithmetic, cached per (N, a)."""
    key = (N, a.str(30))
    if key not in _W:
        n_, p_, ex, blocks, ends = PD[N]
        cex = [lp*(-a*ln).exp() for lp, ln in ex]
        A = [arb(0)]*len(ends); B = [arb(0)]*len(ends); W1 = {}
        if blocks:
            LN = logn(N)
            for idx, j in blocks:
                w = arb(0); w1 = arb(0)
                for i in idx.tolist():
                    v = LOGP[int(p_[i])]*(-a*LN[i]).exp()
                    w += v; w1 += v*(LN[i] - ends[j])
                A[j] = w; B[j + 1] = w1*ends[j + 1]; W1[j] = w1
        Asum = sum(A, arb(0))
        if len(_W) > 4:
            _W.clear()
        _W[key] = (cex, A, B, Asum, W1)
    return _W[key]


def g(t, x):
    """g_t(x) = 1 - exp(-x^2/4t); t = inf gives 0."""
    if t is None:
        return arb(0)
    return 1 - (-(x*x)/(4*t)).exp()


def P_low(N, a, t):
    """Lower bound of sum_{n<=N} Lambda(n) n^-a g_t(log n) at the point t (> 0) or t = None (infinity): 0."""
    if N < 2 or t is None:
        return arb(0)
    n_, p_, ex, blocks, ends = PD[N]
    cex, A, B, Asum, W1 = coef(N, a)
    q = 1/(4*t); h = 2*q
    s = sum(cex, arb(0))
    for c, (lp, ln) in zip(cex, ex):
        s -= c*(-(ln*ln)*q).exp()
    if not blocks:
        return s
    E = [(-(l*l)*q).exp() for l in ends]
    if 2*t < ends[0]*ends[0]:
        # every block lies right of the peak sqrt(2t) of g_t', so the minimum of g_t' on a block is at its right end:
        # sum_j A_j (1 - E_j) + h sum_j W1_j l_{j+1} E_{j+1} = Asum + sum_k E_k (h B_k - A_k)
        s += Asum
        for e, a_, b_ in zip(E, A, B):
            s += e*(h*b_ - a_)
        return s
    for idx, j in blocks:
        x = ends[j]*E[j]; y = ends[j + 1]*E[j + 1]
        mn = x.lower() if x.lower() < y.lower() else y.lower()
        s += A[j]*(1 - E[j]) + W1[j]*h*mn
    return s


def P_exact(N, a, t):
    if N < 2:
        return arb(0)
    n_, p_ = PD[N][:2]
    LN = logn(N)
    s = arb(0)
    for i, p in enumerate(p_.tolist()):
        s += LOGP[p]*(-a*LN[i]).exp()*g(t, LN[i])
    return s


def b_ball(a):
    """b_alpha on the ball a (alpha > 1)."""
    z = acb_series([acb(a), acb(1)], prec=2).zeta().coeffs()
    return 1/a + 1/(a - 1) - PI.log()/2 + (a/2).digamma()/2 + (z[1]/z[0]).real


def b_deriv(a):
    """b'(alpha) on the ball a: -1/a^2 - 1/(a-1)^2 + psi'(a/2)/4 + (zeta'/zeta)'(a)."""
    z = acb_series([acb(a), acb(1)], prec=3).zeta().coeffs()
    s = acb_series([acb(a)/2, acb(1)/2], prec=3).lgamma().coeffs()          # s[2] = psi'(a/2)/8
    return -1/a**2 - 1/(a - 1)**2 + 2*s[2].real + (2*z[2]/z[0] - (z[1]/z[0])**2).real


def b_cell(a1f, a2f):
    """Lower bound of b_alpha on [a1f, a2f] by the mean value form b(m) - r sup|b'|; enclosing b on the whole ball
    instead loses about ten times the width, the poles of 1/(alpha-1) and zeta'/zeta cancelling."""
    m = arb((a1f + a2f)/2)
    r = max(arb(a2f) - m, m - arb(a1f))
    return arb((b_ball(m) - abs(b_deriv(m + r*arb(0, 1)))*r).lower())


ZBIG = 20


def term(c, t, st, spt, side):
    """int_0^inf g_t(x) e^{-cx} dx = (1 - E(z))/c, E(z) = sqrt(pi) z e^{z^2} erfc(z), z = c sqrt t. For z > ZBIG the
    enveloping asymptotic series, 1 - 1/(2z^2) <= E(z) <= 1 - 1/(2z^2) + 3/(4z^4), replaces erfc, which Arb encloses
    only loosely at very large arguments."""
    z = c*st
    if z > ZBIG:
        u = 1/(2*c**3*t)
        return u - 3/(4*c**5*t*t) if side < 0 else u
    return 1/c - spt*(z*z).exp()*z.erfc()


def Abar_low(a, t, K1=40):
    """Lower bound of Abar_a(t): the terms k <= K1, and for k > K1 the bound term_k >= 1/(2 t c^3) - 3/(4 t^2 c^5),
    summed by sum_{k>K1} f(k) >= int_{K1+1}^inf f + f(K1+1)/2 for the convex decreasing f = 1/(2 t c^3), and
    sum_{k>K1} 3/(4 t^2 c^5) <= int_{K1}^inf."""
    if t is None:
        return arb(0)
    st = t.sqrt(); spt = (PI*t).sqrt(); s = arb(0)
    for k in range(1, K1 + 1):
        s += term(a + 2*k, t, st, spt, -1)
    c1 = a + 2*K1 + 2; c0 = a + 2*K1
    tail = 1/(8*t*c1**2) + 1/(4*t*c1**3) - 3/(32*t*t*c0**4)
    if tail > 0:
        s += tail
    return s


def Abar_up(a, t, K1=200):
    """Upper bound of Abar_a(t): the terms k <= K1, and for k > K1, term_k <= 1/(2 t c^3), summed by the midpoint rule
    for the convex f: sum_{k>K1} f(k) <= int_{K1+1/2}^inf f = 1/(8 t (a + 2 K1 + 1)^2)."""
    st = t.sqrt(); spt = (PI*t).sqrt(); s = arb(0)
    for k in range(1, K1 + 1):
        s += term(a + 2*k, t, st, spt, 1)
    return s + 1/(8*t*(a + 2*K1 + 1)**2)


def L(a, t, side=1):
    """The pole term L_a(t) (side = 1: an upper bound, side = -1: a lower bound); L(a, 0) = 1/(a - 1)."""
    c = a - 1
    if t is None:
        return arb(0)
    if t == 0:
        return 1/c
    return term(c, t, t.sqrt(), (PI*t).sqrt(), side)


SQ = {}


def squares(N):
    """(log n)^2 for the exact terms, and a ball containing (log n)^2 for every n of each block (cached)."""
    if N not in SQ:
        n_, p_, ex, blocks, ends = PD[N]
        SQ[N] = ([ln*ln for _, ln in ex], [ends[j].union(ends[j + 1])**2 for _, j in blocks])
    return SQ[N]


def dP_dt(N, a, T):
    """Enclosure of d/dt sum_{n<=N} Lambda(n) n^-a g_t(log n) for every t in the ball T. The weights are positive, so a
    block contributes its weight sum times an enclosure of d/dt g_t over log n in [l0, l1]."""
    if N < 2:
        return arb(0)
    n_, p_, ex, blocks, ends = PD[N]
    cex, A, B, Asum, W1 = coef(N, a)
    ex2, X2 = squares(N)
    q = 1/(4*T)
    s = arb(0)
    for c, x2 in zip(cex, ex2):
        s -= c*x2*(-x2*q).exp()
    for (idx, j), x2 in zip(blocks, X2):
        s -= A[j]*x2*(-x2*q).exp()
    return 4*q*q*s


def dL_dt(a, T):
    """Enclosure of d/dt L_a(t) for every t in the ball T, or None where erfc would be enclosed too loosely."""
    st = T.sqrt(); z = (a - 1)*st
    if not z < ZBIG:
        return None
    Ep = PI.sqrt()*(z*z).exp()*z.erfc()*(1 + 2*z*z) - 2*z
    return -Ep/(2*st)


def lb_heat_mv(N, a1, a2, b, t1, t2):
    """Mean value form of the lower bound on [a1, a2] x [t1, t2], 0 < t1 < t2 < infinity; None if not available."""
    T = t1.union(t2); tm = T.mid(); r = T.rad()
    dL = dL_dt(a1, T)
    if dL is None:
        return None
    dG = dP_dt(N, a2, T) - dL
    return b + Abar_low(a2, t2) + P_low(N, a2, tm) - L(a1, tm) - r*abs(dG).upper()


def lb_heat(N, a1, a2, b, t1, t2):
    """Lower bound of F_N on [a1, a2] x [t1, t2]; t1 may be 0 and t2 None (infinity)."""
    return b + Abar_low(a2, t2) + P_low(N, a2, t2) - L(a1, t1)


def cover_heat(task):
    """Cover [a1, a2] x [0, inf] for F_N; splits the alpha cell if a time cell cannot be certified."""
    N, a1f, a2f, depth = task
    a1, a2 = arb(a1f), arb(a2f)
    b = b_cell(a1f, a2f)
    edges = [arb(0)] + [arb(10)**(e/8) for e in range(-64, 65)] + [None]       # 0, 1e-8 .. 1e8, infinity
    work = [(edges[i], edges[i + 1]) for i in range(len(edges) - 1)]
    cells = 0; worst = None; wcell = None
    while work:
        t1, t2 = work.pop()
        lb = lb_heat(N, a1, a2, b, t1, t2)
        if not lb > 0 and not (t1 == 0 or t2 is None):
            lb2 = lb_heat_mv(N, a1, a2, b, t1, t2)
            if lb2 is not None and lb2 > 0:
                lb = lb2
        if lb > 0:
            cells += 1
            lo = float(lb.lower().str(15, radius=False))
            if worst is None or lo < worst:
                worst, wcell = lo, (a1f, a2f, float(t1.str(10, radius=False)), None if t2 is None else float(t2.str(10, radius=False)))
            continue
        if t1 == 0 or t2 is None or t2/t1 < 1 + 1e-6:
            if depth >= 14:
                return dict(ok=False, N=N, cell=(a1f, a2f, str(t1), str(t2)), lb=lb.str(6), cells=cells)
            mid = (a1f + a2f)/2
            r1 = cover_heat((N, a1f, mid, depth + 1)); r2 = cover_heat((N, mid, a2f, depth + 1))
            ok = r1['ok'] and r2['ok']
            ws = [r for r in (r1, r2) if r.get('worst') is not None]
            w = min(ws, key=lambda r: r['worst']) if ws else None
            # the rest of this alpha cell is redone by the halves: return their union
            return dict(ok=ok, N=N, cells=cells + r1['cells'] + r2['cells'], worst=w['worst'] if w else None,
                        where=w['where'] if w else None, split=True, fail=None if ok else (r1 if not r1['ok'] else r2))
        r = (t2/t1)**(arb(1)/8); pts = [t1*r**i for i in range(8)] + [t2]
        work.extend((pts[i], pts[i + 1]) for i in range(8))
    return dict(ok=True, N=N, cells=cells, worst=worst, where=wcell)


def alpha_cells(est, a_hi, fac, top=10.0):
    """Cells [a_i, a_{i+1}] from a_hi to top, of width fac times the distance to the estimated threshold."""
    pts = [a_hi]
    while pts[-1] < top:
        d = min(max(1e-9, fac*(pts[-1] - est)), 0.25)
        pts.append(min(top, pts[-1] + d))
    return [(pts[i], pts[i + 1]) for i in range(len(pts) - 1)]


# ---------------------------------------------------------------- floating point search for the failing time
def F_float_point(N, al, t):
    """Float value of F_N(al, t) by the ball functions at 64 bits (used only to locate t*)."""
    ctx.prec = 64
    a = arb(al); tt = arb(t)
    v = b_ball(a) + Abar_up(a, tt) + (P_float(N, al, t)) - L(a, tt)
    ctx.prec = 128
    return float(v.mid().str(15, radius=False))


_PF = {}


def P_float(N, al, t):
    if N < 2:
        return arb(0)
    if (N, al) not in _PF:
        m = pn <= N
        _PF[(N, al)] = (np.log(pp[m].astype(float))*np.exp(-al*np.log(pn[m].astype(float))), np.log(pn[m].astype(float)))
    w, x = _PF[(N, al)]
    return arb(float(w @ -np.expm1(-x**2/(4*t))))


def neg_heat(N, al_lo, t_guess):
    ts = [t_guess*math.exp(u) for u in np.linspace(-0.3, 0.3, 241)]
    vals = [F_float_point(N, al_lo, t) for t in ts]
    i = int(np.argmin(vals)); tstar = ts[i]
    a = arb(al_lo); t = arb(tstar)
    up = b_ball(a) + Abar_up(a, t) + P_exact(N, a, t) - L(a, t, -1)
    return dict(alpha=al_lo, t=tstar, float_min=vals[i], upper=up.str(8), negative=bool(up < 0))


# ---------------------------------------------------------------- Thorin form, no primes
def A_th(a, th):
    return ((acb(a + 2, th)/2).digamma().real - ((a + 2)/2).digamma())/2


def pole(a, th):
    c = a - 1
    return th*th/(c*(c*c + th*th))


def cover_thorin(task):
    a1f, a2f, depth = task
    a1, a2 = arb(a1f), arb(a2f)
    ab = a1.union(a2)                     # contains [a1, a2]: both ends are exact doubles
    b = b_cell(a1f, a2f)
    TH = 2000
    work = [(arb(2*i), arb(2*i + 2)) for i in range(TH//2 - 1, -1, -1)]
    cells = 0; worst = None; wcell = None
    # tail theta >= TH: A >= A(alpha, TH), pole <= 1/(alpha - 1)
    lbt = b + A_th(ab, arb(TH)) - 1/(a1 - 1)
    if not lbt > 0:
        return dict(ok=False, cell=(a1f, a2f, 'tail'), lb=lbt.str(6), cells=0)
    while work:
        th1, th2 = work.pop()
        lb = b + A_th(ab, th1) - pole(a1, th2)
        if lb > 0:
            cells += 1
            lo = float(lb.lower().str(15, radius=False))
            if worst is None or lo < worst:
                worst, wcell = lo, (a1f, a2f, float(th1.str(10, radius=False)))
            continue
        if th2 - th1 < arb(1)/2**30:
            if depth >= 16:
                return dict(ok=False, cell=(a1f, a2f, str(th1)), lb=lb.str(6), cells=cells)
            mid = (a1f + a2f)/2
            r1 = cover_thorin((a1f, mid, depth + 1)); r2 = cover_thorin((mid, a2f, depth + 1))
            ws = [r for r in (r1, r2) if r.get('worst') is not None]
            w = min(ws, key=lambda r: r['worst']) if ws else None
            return dict(ok=r1['ok'] and r2['ok'], cells=cells + r1['cells'] + r2['cells'],
                        worst=w['worst'] if w else None, where=w['where'] if w else None, split=True)
        m = (th1 + th2)/2
        work.extend([(m, th2), (th1, m)])
    return dict(ok=True, cells=cells, worst=worst, where=wcell)


def neg_thorin(al_lo):
    ctx.prec = 64
    a = arb(al_lo); b = b_ball(a)
    ths = np.linspace(0.01, 20, 4000)
    vals = [float((b + A_th(a, arb(th)) - pole(a, arb(th))).mid().str(15, radius=False)) for th in ths]
    i = int(np.argmin(vals))
    lo, hi = max(0.01, ths[i] - 0.01), ths[i] + 0.01
    fine = np.linspace(lo, hi, 2001)
    v2 = [float((b + A_th(a, arb(th)) - pole(a, arb(th))).mid().str(15, radius=False)) for th in fine]
    j = int(np.argmin(v2)); th = fine[j]
    ctx.prec = 128
    a = arb(al_lo); x = b_ball(a) + A_th(a, arb(th)) - pole(a, arb(th))
    return dict(alpha=al_lo, theta=float(th), upper=x.str(8), negative=bool(x < 0))


# alpha_hi: printed + 5e-5, or the three decimal bound of the table if smaller, less 1e-12
A_HI = {N: min(float(pr) + 5e-5, math.ceil(ESTIMATE[N]*1000)/1000) - 1e-12 for N, pr in PRINTED.items()}
A_HI['thorin'] = min(float(THORIN_PRINTED) + 5e-5, math.ceil(THORIN_ESTIMATE*1000)/1000) - 1e-12
PD = {}
for N in PRINTED:
    PD[N] = prime_data(N)
print(f'prime data ready, {time.time() - T0:.0f} s', flush=True)

if __name__ == '__main__':
    out = dict(rows=[])
    b10 = b_ball(arb(10))
    out['b10'] = b10.str(10); out['b10_gt_1_9'] = bool(b10 > arb(1)/9)
    print('b_10 =', b10.str(10), '> 1/9:', bool(b10 > arb(1)/9), flush=True)
    tasks = []
    for N, pr in PRINTED.items():
        tasks += [(N, a1, a2, 0) for a1, a2 in alpha_cells(ESTIMATE[N], A_HI[N], FAC[N])]
    th_tasks = [(a1, a2, 0) for a1, a2 in alpha_cells(THORIN_ESTIMATE, A_HI['thorin'], FAC['thorin'])]
    tstar = {1: 0.04062, 10**3: 1.658, 10**5: 4.150, 10**7: 7.518}      # floating point minimisers near the threshold
    with Pool(WORKERS) as pool:
        rh = pool.map_async(cover_heat, sorted(tasks, key=lambda x: -x[0]), chunksize=1)
        rt = pool.map_async(cover_thorin, th_tasks, chunksize=1)
        rn = pool.starmap_async(neg_heat, [(N, float(pr) - 5e-5, tstar[N]) for N, pr in PRINTED.items()], chunksize=1)
        rnt = pool.apply_async(neg_thorin, (float(THORIN_PRINTED) - 5e-5,))
        RH, RT, RN, RNT = rh.get(), rt.get(), rn.get(), rnt.get()
    for (N, pr), neg in zip(PRINTED.items(), RN):
        rs = [r for r in RH if r['N'] == N]
        ok = all(r['ok'] for r in rs)
        w = min((r for r in rs if r.get('worst') is not None), key=lambda r: r['worst'])
        row = dict(N=N, printed=pr, alpha_hi=A_HI[N], alpha_lo=float(pr) - 5e-5, positive_above=ok,
                   alpha_cells=len(rs), cells=sum(r['cells'] for r in rs), min_lower_bound=w['worst'], at=w['where'],
                   negative_below=neg['negative'], neg=neg, failures=[r for r in rs if not r['ok']][:3])
        out['rows'].append(row)
        print(f"primes up to {N}: positive on [{row['alpha_hi']:.6f}, 10] x [0, inf]: {ok} ({row['cells']} cells, "
              f"min lower bound {w['worst']:.3g}); negative at alpha {row['alpha_lo']:.5f}, t {neg['t']:.4g}: "
              f"{neg['negative']} ({neg['upper']})", flush=True)
    okt = all(r['ok'] for r in RT)
    w = min((r for r in RT if r.get('worst') is not None), key=lambda r: r['worst'])
    out['thorin'] = dict(printed=THORIN_PRINTED, alpha_hi=A_HI['thorin'], positive_above=okt, cells=sum(r['cells'] for r in RT),
                         min_lower_bound=w['worst'], at=w['where'], negative_below=RNT['negative'], neg=RNT,
                         failures=[r for r in RT if not r['ok']][:3])
    print(f"Thorin form, no primes: positive on [{A_HI['thorin']:.6f}, 10] x [0, inf]: {okt}; "
          f"negative at {float(THORIN_PRINTED) - 5e-5:.5f}, theta {RNT['theta']:.4f}: {RNT['negative']} ({RNT['upper']})")
    out['all_ok'] = bool(out['b10_gt_1_9'] and okt and RNT['negative'] and all(r['positive_above'] and r['negative_below'] for r in out['rows']))
    out['seconds'] = time.time() - T0
    json.dump(out, open('heat_thresholds.json', 'w'), indent=1)
    print(f"all certified: {out['all_ok']}; {time.time() - T0:.0f} s")
