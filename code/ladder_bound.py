# ladder_bound.py: the one range proof of Theorem thm:ch12:idsd and of the first rungs, the paragraph after the two range
# method in Section ch:heat (ch/ch12.tex). Added for the book in October 2026, from the first author's note on the semicircle
# law and the Thorin framework (misc/reviews/nick-wigner/, Theorem 3 there). Let H >= 100, let every zero with ordinate in
# (0, H] lie on the line, and let 0 <= k <= H/8. A quadruple 1/2 +- delta +- i gamma has gamma > H and phase
# |2 delta gamma t - k phi| < gamma t + 1, so its term in c_k(t) is positive unless gamma t > pi/2 - 1 > 0.57, and then at most
# 2 (gamma^2/99.9)^k e^{-0.285 gamma} times the term gamma_1^{2k} e^{-gamma_1^2 t} of the first zero. Hence at every t > 0
#   c_k(t) >= (1 - S_k(H)) gamma_1^{2k} e^{-gamma_1^2 t},  S_k(H) = sum over the zeros above height H of (gamma^2/99.9)^k e^{-0.285 gamma}.
# Every inequality below is evaluated in ball arithmetic (python-flint, Arb) and counts only if it holds for the whole ball.
#  (1) Constants: gamma_1^2 > 199.79 and pi/2 - 1 > 0.57; for gamma >= 100, (gamma^2 + 1/4)/gamma_1^2 <= gamma^2/99.9, that is
#      gamma^2 (gamma_1^2 - 99.9) >= 24.975, and gamma^2 - 1/4 - gamma_1^2 >= gamma^2/2; both sides of each differ by a function
#      that increases with gamma, so each is checked at gamma = 100.
#  (2) N(T+1) - N(T) < 2 log T for T >= 100. The explicit Riemann von Mangoldt formula [HSW22, (1.7), (2.7) and Corollary 1.4]
#      gives |N(T) - M(T) - 7/8| <= R(T) for T >= e, with M(T) = (T/2pi) log(T/(2 pi e)), R(T) = 2.5167 + 1/(50T) for
#      T <= X0 = 30 610 046 000 and R(T) = 0.1095 log T + 0.2042 log log T + 3.0305 + 1/(50T) < 0.3 log T for T > X0, so
#      N(T+1) - N(T) <= M(T+1) - M(T) + R(T) + R(T+1). While T + 1 <= X0 the margin 2 log T - M(T+1) + M(T) - R(T) - R(T+1) is
#      positive at T = 100 and increases: M'(T+1) - M'(T) = log(1 + 1/T)/(2 pi) <= 1/(2 pi T) and R decreases there, so the
#      derivative of the margin is at least (2 - 1/(2 pi))/T > 0. For T + 1 > X0, R(T) + R(T+1) < 0.6 log(T+1) and
#      M(T+1) - M(T) <= log(T+1)/(2 pi), and (1/(2 pi) + 0.6) log(T+1) < 2 log T, at X0 and beyond since 2/T > 0.76/(T+1).
#      That R(T) < 0.3 log T for T > X0 follows from its value at X0 and 0.3 - 0.1095 - 0.2042/log X0 > 0. (The first version
#      used Trudgian's bound 0.112 log T + 0.278 log log T + 2.510 + 0.2/T, whose published proof is incomplete [BMOR21, HSW22].)
#  (3) With f(T) = (T^2/99.9)^k e^{-0.285 T} and c = 0.285 - 2k/H > 0, f(T) <= f(H) e^{-c (T - H)} for T >= H, and (H + j, H + j + 1]
#      holds fewer than 2 log(H + j) <= 2 (log(H + 1) + j) zeros, so
#      S_k(H) <= 2 f(H) [log(H + 1)/(1 - e^{-c}) + e^{-c}/(1 - e^{-c})^2],
#      evaluated as log10 of the bound for (k, H) = (0, 100), (1, 100) and (10^10, 3 * 10^12).
# Run from any folder: python3 ladder_bound.py; writes ladder_bound.json there.
import json
from flint import arb, acb, ctx

ctx.prec = 200
out = {}


def rec(k, v):
    out[k] = v
    print(k, v, flush=True)


def bstr(x, n=12):
    return x.str(n, radius=True)


pi = arb.pi()
g1 = acb.zeta_zeros(1, 1)[0].imag
g1sq = g1 * g1
rec('(1) constants', {
    'gamma_1': bstr(g1, 20), 'gamma_1^2': bstr(g1sq, 15),
    'gamma_1^2 > 199.79': bool(g1sq > arb('199.79')),
    'pi/2 - 1 > 0.57': bool(pi / 2 - 1 > arb('0.57')),
    '100^2 (gamma_1^2 - 99.9) >= 24.975': bool(arb(10000) * (g1sq - arb('99.9')) >= arb('24.975')),
    '100^2/2 >= 1/4 + gamma_1^2': bool(arb(5000) >= arb('0.25') + g1sq)})


def M(T):
    return T / (2 * pi) * (T / (2 * pi * arb(1).exp())).log()


X0 = arb(30610046000)  # |S(T)| <= 2.5167 for 0 <= T <= X0 [HSW22, (1.7)]


def R(x):  # [HSW22, (1.7), (2.7) and Corollary 1.4]
    if x <= X0:
        return arb('2.5167') + 1 / (50 * x)
    assert x > X0
    return arb('0.1095') * x.log() + arb('0.2042') * x.log().log() + arb('3.0305') + 1 / (50 * x)


def margin(T):
    T = arb(T)
    return 2 * T.log() - (M(T + 1) - M(T)) - R(T) - R(T + 1)


dconst = 2 - 1 / (2 * pi)
m100 = margin(100)
RX0 = arb('0.1095') * X0.log() + arb('0.2042') * X0.log().log() + arb('3.0305') + 1 / (50 * X0)
slope = arb('0.3') - arb('0.1095') - arb('0.2042') / X0.log()
beyond = 2 * X0.log() - (1 / (2 * pi) + arb('0.6')) * (X0 + 1).log()
rec('(2) N(T+1) - N(T) < 2 log T for T >= 100', {
    'margin at T = 100': bstr(m100), 'margin at T = 100 positive': bool(m100 > 0),
    'derivative constant': bstr(dconst), 'derivative constant positive': bool(dconst > 0),
    'second form of R at X0': bstr(RX0), 'below 0.3 log X0': bool(RX0 < arb('0.3') * X0.log()),
    '0.3 - 0.1095 - 0.2042/log X0 positive': bool(slope > 0),
    '2 log X0 - (1/(2 pi) + 0.6) log(X0 + 1)': bstr(beyond), 'beyond X0 positive': bool(beyond > 0 and 1 / (2 * pi) + arb('0.6') < 2),
    'margin at T = 1e3, 1e6, 1e12 (for comparison)': [bstr(margin(10**e), 8) for e in (3, 6, 12)]})


def log10_S_bound(k, H):
    k = arb(k)
    H = arb(H)
    c = arb('0.285') - 2 * k / H
    logf = k * (H * H / arb('99.9')).log() - arb('0.285') * H
    ec = (-c).exp()
    br = (H + 1).log() / (1 - ec) + ec / (1 - ec) ** 2
    return (arb(2).log() + logf + br.log()) / arb(10).log(), c, bool(8 * k <= H)


res = {}
for (k, H, claim) in ((0, 100, '-10.5'), (1, 100, '-8.5'), (10**10, 3 * 10**12, '-1.4e11')):
    L10, c, small_k = log10_S_bound(k, H)
    res['k = %s, H = %s' % (k, H)] = {'log10 of the bound for S_k(H)': bstr(L10, 8), 'c': bstr(c, 8), 'c > 0': bool(c > 0),
                                      'k <= H/8': small_k, 'S_k(H) < 10^' + claim: bool(L10 < arb(claim))}
rec('(3) S_k(H)', res)
json.dump(out, open('ladder_bound.json', 'w'), indent=1)
print('done')
