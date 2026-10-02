# code/firstfail_arb.py (September 2026, from misc/reviews/nick-audit/A-first-failure): certificates, in ball arithmetic
# (python-flint, Arb, 256 bits), for the negative values of Table tab:ch13:first. For each case of firstfail.json, the
# closest pair A < B of consecutive zeros in a window moved to the quadruple 1/2 +- delta +- i g, g = (A+B)/2, it proves
# c_k^F(t) < 0 at the first failing order k found by firstfail.py and at the shortest decimal rounding t of its failure
# time for which the certificate succeeds, where c_k^F(t) = (-1)^k W_F^(k)(t) = sum a^k e^{-at} over one zero of each
# pair {rho, 1-rho} of F, a = -(rho-1/2)^2. That no smaller order fails is computed by firstfail.py, not certified here.
#
# After division by the positive normaliser (k/t)^k e^{-k}, c_k^F(t) is split into four parts.
#   near  the zeros with g - 40 < gamma < g + 40 other than A and B, from the 70 digit strings of zeros_1700.txt (each
#         given the radius 1e-60), terms gamma^{2k} e^{-gamma^2 t}; the cuts g +- 40 are checked to fall between zeros;
#   quad  2 Re a^k e^{-at} with a = g^2 - delta^2 - 2 i delta g;
#   below the zeros with 0 < gamma <= g - 40, all on the line: at most N(g - 40) of them, and since x^k e^{-xt} increases
#         on (0, k/t) and (g - 40)^2 < k/t, each term is at most the weight at x = (g - 40)^2;
#   above the zeros with gamma >= g + 40. With x0 = (g + 40)^2 > k/t and lam = t - k/x0 > 0, concavity of k log x - xt
#         gives x^k e^{-xt} <= x0^k e^{-x0 t} e^{-lam (x - x0)} for x >= x0. The zeros up to H = 3e12 lie on the line
#         [PT21] and contribute at most x0^k e^{-x0 t} (x0 + 1/lam), by partial summation against N(s) <= s^2. A zero above
#         H may lie off the line, but |a| <= gamma^2 + 1/4 and Re a >= gamma^2 - 1/4, so its term is at most e^{t/2} times
#         the weight at gamma^2 + 1/4, and all of them together at most e^{t/2} x0^k e^{-x0 t} e^{-lam (H^2 - x0)} (H^2 + 1/lam).
# The certificate is the ball inequality near + quad + below + above < 0, with below and above upper bounds.
# N(T) is bounded by Trudgian's |N(T) - (T/2pi) log(T/(2 pi e)) - 7/8| <= 0.112 log T + 0.278 log log T + 2.510 + 0.2/T
# (T >= e) [Trudgian14]. For s >= 100 it gives N(s) <= s^2: log s <= s/20 and log log s <= log s there, so the bound is at
# most s^2/(40 pi) + 0.02 s + 3.4. zeros_1700.txt (zeros_hp.py) holds all zeros up to 2197.26, simple and on the line.
# Run from the folder of zeros_1700.txt and firstfail.json; writes firstfail_arb.json. Under a minute on one core.
import json

from flint import acb, arb, ctx

ctx.prec = 256
STR = open('zeros_1700.txt').read().split()
ZERO = [arb(s, '1e-60') for s in STR]
FLT = [float(s) for s in STR]
H = arb(3) * arb(10) ** 12
W = 40
pi = arb.pi()


def N_upper(T):   # Trudgian's upper bound for N(T), T >= e
    return (T / (2 * pi) * (T / (2 * pi * arb(1).exp())).log() + arb(7) / 8 + arb('0.112') * T.log()
            + arb('0.278') * T.log().log() + arb('2.510') + arb('0.2') / T)


def certify(i, delta, k, t):
    A, B = ZERO[i], ZERO[i + 1]
    g = (A + B) / 2
    lo, hi = g - W, g + W
    near_idx = [j for j in range(len(ZERO)) if j not in (i, i + 1) and abs(FLT[j] - (FLT[i] + FLT[i + 1]) / 2) <= W]
    j0, j1 = min(near_idx), max(near_idx)
    ok_cuts = bool(ZERO[j0 - 1] < lo and lo < ZERO[j0] and ZERO[j1] < hi and hi < ZERO[j1 + 1]
                   and j0 <= i and i + 1 <= j1 and len(near_idx) == j1 - j0 - 1)
    kk = arb(k)
    L0 = kk * (kk / t).log() - kk                       # log of the normaliser (k/t)^k e^{-k}
    lw = lambda x: kk * x.log() - x * t - L0           # log of the normalised weight x^k e^{-xt}
    near = sum((lw(ZERO[j] ** 2).exp() for j in near_idx), arb(0))
    a = acb(g * g - delta * delta, -2 * delta * g)
    quad = 2 * (kk * a.log() - a * t - L0).exp().real
    ok_below = bool(lo * lo < kk / t)
    below = N_upper(lo) * lw(lo * lo).exp()
    x0 = hi * hi
    lam = t - kk / x0
    ok_above = bool(lam > 0 and hi > 100)
    above = lw(x0).exp() * (x0 + 1 / lam + (t / 2 - lam * (H * H - x0)).exp() * (H * H + 1 / lam))
    total = near + quad + below + above
    ok = ok_cuts and ok_below and ok_above
    return dict(ok=ok, certified=bool(total < 0) and ok, total=total, near=near, quad=quad, below=below, above=above,
                n_near=len(near_idx), zero_range=[j0 + 1, j1 + 1])


def s(x, n=10):
    return x.str(n, radius=True)


out = []
for c in json.load(open('firstfail.json'))['cases']:
    i = c['zero_numbers'][0] - 1
    assert STR[i] == c['A'] and STR[i + 1] == c['B']
    num, den = c['delta'].split('/')
    delta = arb(num) / arb(den)
    row = dict(window=c['window'], delta=c['delta'], zero_numbers=c['zero_numbers'], k=c['kstar'], certified=False)
    for d in range(2, 16):
        tstr = ('%.' + str(d) + 'f') % c['t_star']
        r = certify(i, delta, c['kstar'], arb(tstr))
        if r['certified']:
            row.update(t=tstr, certified=True, total=s(r['total']), near=s(r['near'], 15), quad=s(r['quad'], 15),
                       below=s(r['below'], 5), above=s(r['above'], 5), n_near=r['n_near'], near_zero_numbers=r['zero_range'])
            break
    out.append(row)
    print(json.dumps(row), flush=True)
json.dump(dict(prec=ctx.prec, all_certified=all(r['certified'] for r in out), cases=out),
          open('firstfail_arb.json', 'w'), indent=1)
print('all certified:', all(r['certified'] for r in out))
