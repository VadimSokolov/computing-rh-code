# window_orders.py: certificate for Theorem thm:ch12:window in Section sec:ch12:reach (ch/ch12.tex). Added for the book in
# October 2026, from the first case of the proof of the unconditional strip in the authors' note on the reciprocal Xi ladder
# (misc/reviews/nick-ladder/, item 163 of the authors' report). With H = 3e12, the height to which RH is verified [PT21],
# the theorem says that c_k(t) = (-1)^k W^(k)(t) > 0 for every t >= 1.29 and every integer k <= t(H - 23)^2: with
# u = sqrt(k/t), a zero on the line in [u, u + 15] dominates every zero above H, whatever their number, by the Gaussian
# factor of the weight v^{2k} e^{-t v^2}, which peaks at v = u. The proof uses three numerical facts.
#  1. Every interval [u, u + 15] with 0 <= u <= H - 23 contains the ordinate of a zero. For u <= 1000 this follows from the
#     zeros below height 1015, computed by Arb (acb.zeta_zeros: isolated, counted by Turing's method, refined in ball
#     arithmetic): gamma_1 < 15 and consecutive ordinates differ by less than 7. For u >= 1000 it follows from the explicit
#     Riemann von Mangoldt formula of Hasanalizade, Shen and Wong [HSW22, (1.7), (2.7) and Corollary 1.4],
#     |N(T) - (T/2pi) log(T/(2 pi e)) - 7/8| <= E(T) for T >= e, with E(T) = 2.5167 + 1/(50T) for T <= X0 = 30 610 046 000
#     and E(T) = 0.1095 log T + 0.2042 log log T + 3.0305 + 1/(50T) < 0.3 log T for T > X0. (The first version used
#     Trudgian's bound 0.112 log T + 0.278 log log T + 2.510 + 0.2/T, whose published proof is incomplete [BMOR21, HSW22].)
#     It gives N(u + 15) - N(u) >= g(u) = (15/2pi) log(u/2pi) - E(u) - E(u + 15), since the main term has derivative
#     (1/2pi) log(T/2pi), increasing in T. While u + 15 <= X0, E decreases, so g increases, and g(1000) > 4 gives
#     N(u + 15) - N(u) > 4. For u + 15 > X0, E(u) and E(u + 15) are below 0.3 log(u + 15), so g(u) exceeds
#     (15/2pi) log(u/2pi) - 0.6 log(u + 15), which increases (15/2pi > 0.6) and is above 38 at u = X0 - 15. That E(T) < 0.3 log T
#     for T > X0 follows from its value at X0 and from 0.3 - 0.1095 - 0.2042/log X0 > 0, which makes the difference increase.
#  2. t(D^2 - 2*15^2 - 1/2) > log(H^2 + 1/t + H/(tD)) at t = 1.29 and D = 23. The left side increases and the right side
#     decreases in t and in D, so the inequality holds for every t >= 1.29 and every D >= 23; with D = H - u it says that
#     the zero in [u, u + 15], which contributes at least f(u) e^{-450t}, outweighs the zeros above H, which contribute at
#     most e^{t/2} f(u) e^{-tD^2} (H^2 + 1/t + H/(tD)) in modulus (partial summation with N(r) <= r^2).
#  3. 1.29 (H - 23)^2 > 1.16e25, the range of orders at every t > 0 once Theorem thm:ch12:cmint(b) covers t <= 1.29.
# Every inequality is evaluated in ball arithmetic (python-flint, Arb) and counts only if it holds for the whole ball.
# Run from any folder: python3 window_orders.py; writes window_orders.json there. A few seconds.
import json
import time

from flint import arb, acb, ctx

ctx.prec = 128
pi = arb.pi()
out, holds = {}, []
T0 = time.time()


def s(x, n=12):
    return x.str(n, radius=False) if isinstance(x, arb) else x


def check(name, cond, **vals):
    holds.append(bool(cond))
    out[name] = dict(holds=bool(cond), **{k: s(v) for k, v in vals.items()})
    print(name, bool(cond), round(time.time() - T0, 1), flush=True)


H, h, D, t1, U = arb(3) * arb(10) ** 12, arb(15), arb(23), arb('1.29'), arb(1000)

# 1. Gaps: the zeros below height 1015, and the explicit Riemann von Mangoldt formula above u = 1000
n = 670
Z = acb.zeta_zeros(1, n)
half = arb(1) / 2
assert all(z.real == half for z in Z)
G = [z.imag for z in Z]
gaps = [b - a for a, b in zip(G[:-1], G[1:])]
gmax = max(gaps, key=lambda x: x.upper())
check('zeros_below_1015', G[0] < h and G[-1] > U + h and all(g < 7 for g in gaps),
      zeros=n, gamma_1=G[0], last_ordinate=G[-1], largest_gap=gmax)


X0 = arb(30610046000)  # |S(T)| <= 2.5167 for 0 <= T <= X0 [HSW22, (1.7)]


def E(T):  # the error term for N(T), T >= e [HSW22, (1.7), (2.7) and Corollary 1.4]
    T = arb(T)
    if T <= X0:
        return arb('2.5167') + 1 / (50 * T)
    assert T > X0
    return arb('0.1095') * T.log() + arb('0.2042') * T.log().log() + arb('3.0305') + 1 / (50 * T)


g1000 = h / (2 * pi) * (U / (2 * pi)).log() - E(U) - E(U + h)
EX0 = arb('0.1095') * X0.log() + arb('0.2042') * X0.log().log() + arb('3.0305') + 1 / (50 * X0)
slope = arb('0.3') - arb('0.1095') - arb('0.2042') / X0.log()  # T (0.3 log T - E(T))' exceeds this for T >= X0
gX = h / (2 * pi) * ((X0 - h) / (2 * pi)).log() - arb('0.6') * X0.log()  # the lower bound for g at u = X0 - 15
check('count_gap_above_1000', g1000 > 4 and EX0 < arb('0.3') * X0.log() and slope > 0 and
      arb('2.5167') + 1 / (50 * U) < arb('0.3') * X0.log() and h / (2 * pi) > arb('0.6') and gX > 38,
      g_at_1000=g1000, E_at_X0=EX0, three_tenths_log_X0=arb('0.3') * X0.log(), slope_factor=slope,
      g_lower_at_X0_minus_15=gX, slope_of_main_term=h / (2 * pi))

# 2. The zeros above H against the zero in [u, u + 15], at t = 1.29 and D = H - u = 23, with the two numbers the proof
# prints: t(D^2 - 450.5) >= 1.29 * 78.5 > 101, and log(H^2 + 1 + H/23) < 57.5, which bounds the right side for t >= 1.29
lhs = t1 * (D * D - 2 * h * h - half)
rhs = (H * H + 1 + H / D).log()
check('domination_at_t_1.29_D_23', D * D - 2 * h * h - half == arb(157) / 2 and lhs > 101 and rhs < arb('57.5') and 1 / t1 < 1,
      lhs=lhs, rhs=rhs, margin=lhs - rhs)

# 3. The range of orders at every t > 0
K = t1 * (H - D) ** 2
check('orders_at_every_t', K > arb('1.16e25'), K=K)

out['all_hold'] = all(holds)
out['seconds'] = round(time.time() - T0, 1)
json.dump(out, open('window_orders.json', 'w'), indent=1)
print('all inequalities hold:', all(holds), '| checks:', len(holds))
print(json.dumps(out, indent=1))
