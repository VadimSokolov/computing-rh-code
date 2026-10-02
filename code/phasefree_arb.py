# code/phasefree_arb.py (October 2026): certificates, in ball arithmetic (python-flint, Arb), for the phase free
# thresholds t0(u) of Section ch:lessons that phasefree.py computes in floating point. With
#   m(r) = Re psi(1/4 + ir/2) - log pi,   A(u,t) = (1/2pi) int_R e^{-2t(r-u)^2} m(r) dr,   pol(u,t) = 2 e^{t/2-2tu^2} cos 2tu,
#   S(t) = sum_{n>=2} Lambda(n) n^{-1/2} e^{-(log n)^2/8t} = 2 sqrt(2 pi t) (e^{t/2} + (1/4pi) int_R e^{-2tr^2} m(r) dr - W(2t)),
# the bound |cos(u log n)| <= 1 proves G(u,t) + G(-u,t) > 0 wherever D(u,t) = sqrt(2 pi t)(pol + A)(u,t) - S(t) > 0.
# For each height u it proves D(u,t) > 0 for every 0 < t <= t_lo and D(u,t_hi) < 0, so the first time t0(u) at which
# the bound fails lies in (t_lo, t_hi].
# Facts used.
#  (a) m is even and increasing in |r|: Re psi(x+iy) = psi(x) + sum_{n>=0} y^2/((n+x)((n+x)^2+y^2)). The same series
#      gives psi(1/4) - log pi <= m(r) <= 1.13 + log max(|r|/2, 1), so |m(r)| <= 5.4 + log(1 + 2|r|), which bounds the
#      tails |r - u| > L = 12/sqrt t of the integrals (at most 1e-120 here).
#  (b) S increases with t (termwise), so S(t) <= S(t') for t <= t'. For t <= 1/1600 the prime sum itself is tiny:
#      (log n)^2 >= log 2 log n and Lambda(n) n^{-1/2} <= log(n)/sqrt(n) <= 1 give S(t) <= sum_{n>=2} n^{-c}
#      <= 2^{-c}(1 + 2/(c-1)), c = log 2/(8t). (The explicit formula is not used there: W(2t) is not small.)
#  (c) sqrt(2 pi t) A(u,t) = E m(u + sZ)/2 with s = 1/(2 sqrt t). Since m is even and increasing in |r|, and |u + sZ|
#      increases stochastically with u >= 0, E m(u + sZ) increases with u >= 0. It need not increase with s at fixed u
#      (m is concave at large r), but E m(u + sZ) >= E m(sZ) >= E m(s_a Z) for u >= 0 and s >= s_a, the second since
#      |sZ| >= |s_a Z| pointwise. With t_a the double nearest 1/1600 and s_a = 1/(2 sqrt t_a), about 20, this gives
#      sqrt(2 pi t) A >= M/2, M = E m(s_a Z), for t <= t_a, and sqrt(2 pi t)|pol| <= sqrt(2 pi) e^{t/2 - 1/2}/u there
#      (the maximum of sqrt t e^{-2tu^2} is at t = 1/(4u^2)).
#  (d) dA/dt = -(1/2pi) int 2(r-u)^2 e^{-2t(r-u)^2} m(r) dr, and m >= 0 for |r| >= 7 (m(7) > 0 is checked), so for
#      u > 7, dA/dt <= k(t) = (|m(0)|/2pi) 28 (u+7)^2 e^{-2t(u-7)^2}, decreasing in t. On a cell [t_j, t_j+1],
#      D >= sqrt(2 pi t_j)(A(u, t_j+1) - k(t_j)(t_j+1 - t_j)) - sqrt(2 pi t_j+1) 2 e^{t_j+1/2 - 2 t_j u^2} - S(t_j+1),
#      by (b) and A >= 0. The cells cover [1/1600, t_lo] and are bisected until this bound is positive.
#  (e) W(2t) = sum_gamma e^{-2t gamma~^2} over one zero of each pair: the zeros below H = 3e12 lie on the line [PT21] and
#      give positive terms, those above have |e^{-2t gamma~^2}| <= e^{t/2} e^{-2t gamma^2}; with N(s) <= s^2 (s >= 14)
#      and partial summation, 0 - e^{t/2}(1 + 2tH^2)e^{-2tH^2}/(2t) <= W(2t) <= e^{t/2}(1 + 392t)e^{-392t}/(2t).
# The integrals are computed by acb.integral (Gauss-Legendre with rigorous error bounds, the integrand analytic in
# |Im r| < 1/2, the poles of the two digamma factors being at r = +-i(1/2 + 2k)).
# The item cmint checks only the positivity part at u = 2.78e8 < sqrt(10^17/1.29), for Theorem thm:ch12:cmint; by (c),
# E m(u + sZ) increases with u, and at these t the bounds of (c) and (d) for the polar term and k(t) decrease with u,
# so the proof of positivity up to t = 1.29 holds at every larger height as well.
# Also S(6), with 2 pi e^{2 S(6)}, which only approximates the height beyond which D(u, 6) > 0 (sqrt(2 pi t) A is
# log(u/2 pi)/2 only up to O(u^-2)), and the item t6, which proves G(u, 6) > 0 for every u >= 1e213 whatever the zeros
# are. G(u,t) >= D(u,t)/sqrt(2 pi t) - B(u,t), where |G(-u,t)| <= B(u,t) = e^{t/2}((u+14)^2 + 1/(2t))e^{-2t(u+14)^2}
# (partial summation with N(s) <= s^2, Section sec:ch12:ray), and |pol| <= 2 e^{t/2 - 2tu^2}. At u = 1e213 it checks
# sqrt(12 pi)(A(u,6) - 2 e^{3 - 12u^2} - B(u,6)) - S(6) > 0; the left side increases with u, since A does by (c) and
# the two envelopes decrease, so the bound holds at every larger u.
# Reads firstfail.json (the heights of the moved pairs); writes phasefree_arb.json.
import json
import math
import sys
from multiprocessing import Pool

from flint import acb, arb, ctx

H = arb(3)*arb(10)**12


def setup(prec):
    ctx.prec = prec
    global PI, LOGPI, Q, I
    PI = arb.pi(); LOGPI = PI.log(); Q = acb(1)/4; I = acb(0, 1)


def m_acb(z):
    w = I*z/2
    return ((Q + w).digamma() + (Q - w).digamma())/2 - LOGPI


def gauss_m_integral(u, t, L):
    """int_{-L}^{L} e^{-2t s^2} m(u + s) ds, rigorous, plus the tail bound for |s| > L from (a)."""
    f = lambda s, analytic: (-2*t*s*s).exp()*m_acb(u + s)
    tol = arb(2)**(-70)
    v = acb.integral(f, -L, L, rel_tol=tol, abs_tol=tol)
    st = (2*t).sqrt()
    tail = (arb('5.4') + (1 + 2*abs(u)).log() + (1 + 2*L).log())*(PI/(2*t)).sqrt()*(L*st).erfc() + (-2*t*L*L).exp()/(2*t)
    return v.real + arb(0, 1)*tail


def W_bounds(t):
    lo = -(t/2).exp()*(1 + 2*t*H*H)*(-2*t*H*H).exp()/(2*t)
    hi = (t/2).exp()*(1 + 392*t)*(-392*t).exp()/(2*t)
    return lo, hi


_S = {}


def S_ball(tf):
    """S(t) at the float t (exactly), as a ball containing the true value."""
    if tf not in _S:
        t = arb(tf)
        L = arb(12)/t.sqrt()
        arch = gauss_m_integral(arb(0), t, L)/(4*PI)
        wl, wh = W_bounds(t)
        c = 2*(2*PI*t).sqrt()
        lo = c*((t/2).exp() + arch - wh); hi = c*((t/2).exp() + arch - wl)
        _S[tf] = lo.union(hi)
    return _S[tf]


def A_ball(u, tf):
    t = arb(tf)
    L = arb(12)/t.sqrt()
    return gauss_m_integral(u, t, L)/(2*PI)


def D_ball(u, tf):
    t = arb(tf)
    pol = 2*(t/2 - 2*t*u*u).exp()*(2*t*u).cos()
    return (2*PI*t).sqrt()*(pol + A_ball(u, tf)) - S_ball(tf)


def run(item):
    name, us, t_lo, t_hi = item
    digits = int(math.log10(float(us))) + 1
    setup(128 + int(3.33*digits))
    u = arb(us)
    out = dict(u=us)
    # float location of t0 by bisection on the midpoints (only to place the bracket)
    a, b = 0.05, 6.0
    for _ in range(40):
        c = (a + b)/2
        if float(D_ball(u, c).mid().str(20, radius=False)) > 0:
            a = c
        else:
            b = c
    out['t0_float'] = (a + b)/2
    if t_lo is None:
        t_lo = math.floor(((a + b)/2 - 1e-4)*1e4)/1e4
        t_hi = math.ceil(((a + b)/2 + 1e-4)*1e4)/1e4
    out['t_lo'], out['t_hi'] = t_lo, t_hi
    # (c) small t
    ta = 1/1600
    M = gauss_m_integral(arb(0), arb(ta), arb(12)/arb(ta).sqrt())*(2*arb(ta)/PI).sqrt()          # E m(20 Z)
    c = arb(2).log()/(8*arb(ta))
    S_small = arb(2)**(-c)*(1 + 2/(c - 1))
    small = M/2 - (2*PI).sqrt()*(arb(ta)/2 - arb(1)/2).exp()/u - S_small
    out['small_t'] = dict(M=M.str(8), S_bound=S_small.str(3), bound=small.str(8), ok=bool(small > 0))
    # (d) cells on [ta, t_lo]
    m0 = m_acb(acb(0)).real; m7 = m_acb(acb(7)).real
    assert m7 > 0 and u > 7
    kap = lambda tt: abs(m0)/(2*PI)*28*(u + 7)**2*(-2*tt*(u - 7)**2).exp()
    pts = [ta]
    while pts[-1] < t_lo:
        pts.append(min(t_lo, pts[-1]*1.1))
    stack = [(pts[i], pts[i + 1]) for i in range(len(pts) - 2, -1, -1)]
    cells = 0; worst = None; fail = None
    while stack:
        x, y = stack.pop()
        X, Y = arb(x), arb(y)
        lb = ((2*PI*X).sqrt()*(A_ball(u, y) - kap(X)*(Y - X)) - (2*PI*Y).sqrt()*2*(Y/2 - 2*X*u*u).exp() - S_ball(y))
        if lb > 0:
            cells += 1
            v = float(lb.lower().str(15, radius=False))
            if worst is None or v < worst[0]:
                worst = (v, x, y)
            continue
        if y - x < 1e-9:
            fail = (x, y, lb.str(6)); break
        z = math.sqrt(x*y)
        stack += [(z, y), (x, z)]
    out['cells'] = cells; out['worst'] = worst; out['fail'] = fail
    Dhi = D_ball(u, t_hi)
    out['D_at_t_hi'] = Dhi.str(8)
    out['positive_up_to_t_lo'] = bool(out['small_t']['ok'] and fail is None)
    out['negative_at_t_hi'] = bool(Dhi < 0)
    out['ok'] = out['positive_up_to_t_lo'] and (name == 'cmint' or out['negative_at_t_hi'])
    print(json.dumps(out), flush=True)
    return name, out


def s6():
    setup(192)
    S = S_ball(6.0)
    lg = (2*PI).log()/arb(10).log() + 2*S/arb(10).log()
    return dict(S6=S.str(15), log10_height=lg.str(10))


def t6():
    us = '1e213'
    setup(128 + int(3.33*214))
    u = arb(us); t = arb(6)
    c = (2*PI*t).sqrt()
    arch = c*A_ball(u, 6.0)
    pol = 2*(t/2 - 2*t*u*u).exp()
    B = (t/2).exp()*((u + 14)**2 + 1/(2*t))*(-2*t*(u + 14)**2).exp()
    S = S_ball(6.0)
    lb = arch - c*(pol + B) - S
    return dict(u=us, t=6, arch=arch.str(15), S6=S.str(15), polar_reflected_below_1e_100=bool(c*(pol + B) < arb(10)**-100),
                margin=lb.str(10), ok=bool(lb > 0))


if __name__ == '__main__':
    heights = [repr(c['g']) for c in json.load(open('firstfail.json'))['cases'] if c['delta'] == '1/4']
    items = [('quad', h, None, None) for h in heights]
    items += [('1e3', '1e3', 0.445, 0.455), ('3e12', '3e12', 1.715, 1.725), ('1e50', '1e50', 3.635, 3.645),
              ('cmint', '2.78e8', 1.29, 1.35)]
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else len(items)
    with Pool(workers) as pool:
        r6 = pool.apply_async(s6)
        r7 = pool.apply_async(t6)
        res = pool.map(run, items, chunksize=1)
        S6 = r6.get(); T6 = r7.get()
    out = dict(rows=[dict(name=n, **o) for n, o in res], S6=S6, t6=T6)
    out['all_ok'] = all(o['ok'] for _, o in res) and T6['ok']
    json.dump(out, open('phasefree_arb.json', 'w'), indent=1)
    for n, o in res:
        print('%-6s u = %-22s t0 in (%s, %s]: positive up to t_lo %s (%d cells), D(t_hi) = %s'
              % (n, o['u'], o['t_lo'], o['t_hi'], o['positive_up_to_t_lo'], o['cells'], o['D_at_t_hi']))
    print('S(6) =', S6['S6'], ' log10 of 2 pi e^{2 S(6)} =', S6['log10_height'])
    print('t = 6, u = %s: sqrt(12 pi) A = %s, S(6) = %s, margin %s: G(u, 6) > 0 for every u >= 1e213 %s'
          % (T6['u'], T6['arch'], T6['S6'], T6['margin'], T6['ok']))
    print('all certified:', out['all_ok'])
