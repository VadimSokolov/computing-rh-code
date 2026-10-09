# Certificate for Proposition bp:prop:valpha: v_alpha(theta) > b_alpha = v_alpha(0) for every theta != 0 and every
# basepoint alpha >= 1, with no hypothesis on the zeros, and for every alpha > 7/8 if no zero has real part greater
# than 7/8 (Part E).  Ball arithmetic (python-flint).
#
# v_alpha(theta) = Re(xi'/xi)(alpha + i theta) = sum_rho (alpha-beta)/((alpha-beta)^2 + (theta-gamma)^2), and
# b_alpha' = sum_rho (gamma^2-(alpha-beta)^2)/((alpha-beta)^2+gamma^2)^2 > 0 for alpha < gamma_1, so b increases.
# Near theta = 0, v_alpha is strictly convex on |theta| <= gamma_1 - alpha/sqrt3 (each Cauchy kernel is convex
# where |theta-gamma| >= (alpha-beta)/sqrt3), so for alpha <= 9 the work is beyond theta0(alpha) = gamma_1 - alpha/sqrt3.
#
# Part A, 1 <= alpha <= 3: the zeros of v1min.py (all on the line, every theta in [12.4, 2000] within D = 3.4437 of
# an ordinate) give v_alpha(theta) >= r(alpha) = (alpha-1/2)/((alpha-1/2)^2 + D^2) on [theta0, 2000], and the
# window count of the proof (at least 20 zeros in (theta-15, theta+15], sum of alpha-beta over them = n(alpha-1/2))
# gives v_alpha(theta) > w(alpha) = 20(alpha-1/2)/(alpha^2+225) for theta >= 2000.  r, w and b all increase on
# [1,3], so on a cell [a_k, a_{k+1}] it suffices that r(a_k) > b(a_{k+1}) and w(a_k) > b(a_{k+1}).
#
# Part B, 3 <= alpha <= 7: by the Euler product, Re(zeta'/zeta)(alpha+i theta) >= (zeta'/zeta)(alpha), and the
# pole terms Re(1/(alpha+i theta) + 1/(alpha-1+i theta)) are positive, so
#   v_alpha(theta) - b_alpha >= L(alpha,theta) = (Re psi((alpha+i theta)/2) - psi(alpha/2))/2 - 1/alpha - 1/(alpha-1),
# which increases with theta.  On a cell, theta0(alpha) >= th = theta0(a_{k+1}), and |dL/dalpha| <= M(a_k) with
# M(a) = 1/a + 3/a^2 + 1/(a-1)^2 (from psi'(x) <= 1/x + 1/x^2), so L(a_k, th) - M(a_k)(a_{k+1}-a_k) > 0 suffices.
#
# Part C, alpha >= 7: the Levy form gives, with q_c(theta) = theta^2/(c(c^2+theta^2)) and the prime terms dropped,
#   v_alpha(theta) - b_alpha >= sum_{k>=1} q_{alpha+2k}(theta) - q_{alpha-1}(theta) >= Phi(theta^2),
#   Phi(s) = log(1+s/A^2)/4 - s/(B(B^2+s)),  A = alpha+2,  B = alpha-1,
# since q_c decreases in c and (1/2) int_A^oo q_c dc = log(1+theta^2/A^2)/4.  Phi(0) = 0, and Phi'(s) >= 0 exactly when
# (B^2+s)^2 >= 4B(A^2+s), the two sides differing by a function that increases in s once B >= 2.
# alpha >= 9: p(alpha) = B^3 - 4A^2 is positive at 9 and increasing (p' > 0 at 9, p'' > 0), so Phi(s) > 0 for s > 0.
# 7 <= alpha <= 9: convexity covers |theta| <= theta0(9); beyond, s >= s9 = theta0(9)^2, and with 6 <= B <= 8, A <= 11,
# (36+s9)^2 > 32(121+s9) makes Phi increase on [s9, oo), where Phi(s9) >= log(1+s9/121)/4 - s9/(6(36+s9)) > 0.
#
# Part D: b_alpha, the mean E T_alpha^circ = (log xi)''(alpha)/2 and the coefficient c(alpha) = -(xi'/xi)''(alpha)/2 of
# theta^2 in v_alpha - b_alpha, from the Taylor series of log xi at alpha.
#
# Part E, 7/8 <= alpha <= 1, for the conditional part of the proposition: if no zero has real part greater than 7/8
# (the half plane announced in October 2026), every zero has alpha - beta in (0, alpha) for alpha > 7/8, which is all
# that Part A uses.  theta0(alpha) >= theta0(1) > 13.5, so convexity covers (0, theta0(alpha)], the zeros of v1min.py
# cover [theta0(alpha), 2000] and the window count covers theta >= 2000; r, w and b increase on [7/8, 1] as on [1, 3]
# (|alpha - beta| < 1 < gamma_1).  b(alpha) for alpha < 1 comes from the same formula, zeta having no zeros on (0, 1).
#
# Writes valpha_min.json.  Run on Hopper: bash misc/tools/hopper_run.sh -c 1 -t 60 -g valpha_min.json valpha_min.py
import json, time
from flint import arb, acb, acb_series, ctx

ctx.prec = 160
t0 = time.time()
G1 = acb.zeta_zero(1).imag              # gamma_1 as a certified ball (Arb isolates the zero and counts by Turing's method)
D = arb('3.4437')
SQ3 = arb(3).sqrt()
LOGPI = arb.pi().log()
H = 256                                   # cells of width 1/256, exact in binary


def b(alpha):                             # b_alpha = (xi'/xi)(alpha), alpha > 1
    s = acb(alpha)
    z0, z1 = acb_series([s, acb(1)], prec=2).zeta().coeffs()[:2]
    return (1 / s + 1 / (s - 1) - LOGPI / 2 + (s / 2).digamma() / 2 + z1 / z0).real


def b1():                                 # b_1 = 1 + gamma_E/2 - log(4 pi)/2, the first Keiper Li coefficient
    return 1 + arb.const_euler() / 2 - (4 * arb.pi()).log() / 2


def L(alpha, theta):
    z = acb(alpha, theta) / 2
    return (z.digamma().real - (arb(alpha) / 2).digamma()) / 2 - 1 / arb(alpha) - 1 / (arb(alpha) - 1)


def lower(x):                             # rigorous lower end of a ball, as a float
    return float(x.mid() - x.rad())


out = {'prec_bits': ctx.prec, 'cell': '1/256'}

# Part A
grid = [arb(1) + arb(k) / H for k in range(0, 2 * H + 1)]        # 1 .. 3
bv = [b1()] + [b(a) for a in grid[1:]]
r = lambda a: (a - arb('0.5')) / ((a - arb('0.5')) ** 2 + D ** 2)
w = lambda a: 20 * (a - arb('0.5')) / (a ** 2 + 225)
mr = min(lower(r(grid[k]) - bv[k + 1]) for k in range(2 * H))
mw = min(lower(w(grid[k]) - bv[k + 1]) for k in range(2 * H))
mono = all(lower(bv[k + 1] - bv[k]) > 0 for k in range(2 * H))
out['A'] = {'range': [1, 3], 'min_r_minus_b': mr, 'min_w_minus_b': mw, 'b_increasing_on_grid': mono,
            'ok': mr > 0 and mw > 0}
out['b'] = {str(a): str(v) for a, v in [(1, bv[0])] + [(x / 4, b(arb(x) / 4)) for x in (5, 6, 7, 8, 10, 12, 16, 20, 24, 28)]}

# Part B
res, first_fail = [], None
k = 0
while True:
    ak, ak1 = arb(3) + arb(k) / H, arb(3) + arb(k + 1) / H
    th = G1 - ak1 / SQ3
    M = 1 / ak + 3 / ak ** 2 + 1 / (ak - 1) ** 2
    lb = lower(L(ak, th) - M / H)
    res.append(lb)
    if lb <= 0:
        first_fail = 3 + k / H
        break
    k += 1
    if k > 6 * H:
        break
out['B'] = {'from': 3, 'certified_to': 3 + len([x for x in res if x > 0]) / H, 'first_failing_cell_start': first_fail,
            'min_margin_on_3_7': min(res[:4 * H]), 'margin_at_3': res[0], 'margin_at_7': res[4 * H - 1] if len(res) >= 4 * H else None}

# Part C
p = lambda a: (a - 1) ** 3 - 4 * (a + 2) ** 2
dp = lambda a: 3 * (a - 1) ** 2 - 8 * (a + 2)
th9 = G1 - 9 / SQ3
s9 = th9 ** 2
mono = (36 + s9) ** 2 - 32 * (121 + s9)
phi9 = (1 + s9 / 121).log() / 4 - s9 / (6 * (36 + s9))
cvals = [lower(x) for x in (p(arb(9)), dp(arb(9)), mono, phi9)]
out['C'] = {'p_at_9': cvals[0], 'dp_at_9': cvals[1], 'theta0_9': float(th9.mid()), 'monotone_margin': cvals[2],
            'Phi_lower_on_7_9': cvals[3], 'ok': all(x > 0 for x in cvals)}


# Part D
def logxi_coeffs(alpha, n=4):             # Taylor coefficients of log xi(alpha + x) - log(1/2), alpha > 1
    s = acb_series([acb(alpha), acb(1)], prec=n)
    return (s.log() + (s - 1).log() - s * (LOGPI / 2) + (s / 2).lgamma() + s.zeta().log()).coeffs()


out['D'] = {}
for a in ('1.25', '1.5', '1.75', '2', '3', '5', '7', '9'):
    cf = logxi_coeffs(arb(a))
    out['D'][a] = {'b': str(cf[1].real), 'mean_T_circ': str(cf[2].real), 'c': str((-3 * cf[3]).real)}

# Part E
grid = [arb(7) / 8 + arb(k) / H for k in range(0, H // 8 + 1)]   # 7/8 .. 1
bv = [b(a) for a in grid[:-1]] + [b1()]
mr = [lower(r(grid[k]) - bv[k + 1]) for k in range(H // 8)]
mw = [lower(w(grid[k]) - bv[k + 1]) for k in range(H // 8)]
mono = all(lower(bv[k + 1] - bv[k]) > 0 for k in range(H // 8))
th1 = lower(G1 - 1 / SQ3 - arb('13.5'))
out['E'] = {'range': ['7/8', 1], 'hypothesis': 'no zero with real part greater than 7/8',
            'min_r_minus_b': min(mr), 'at_alpha_r': 0.875 + mr.index(min(mr)) / H,
            'min_w_minus_b': min(mw), 'at_alpha_w': 0.875 + mw.index(min(mw)) / H,
            'b_increasing_on_grid': mono, 'theta0_1_minus_13.5': th1,
            'b_7/8': str(bv[0]), 'r_7/8': str(r(grid[0])), 'w_7/8': str(w(grid[0])),
            'ok': min(mr) > 0 and min(mw) > 0 and mono and th1 > 0}

out['theta0'] = {'alpha=1': float((G1 - 1 / SQ3).mid()), 'alpha=3': float((G1 - 3 / SQ3).mid()),
                 'alpha=7': float((G1 - 7 / SQ3).mid())}
out['seconds'] = round(time.time() - t0, 1)
print(json.dumps(out, indent=1))
json.dump(out, open('valpha_min.json', 'w'), indent=1)
