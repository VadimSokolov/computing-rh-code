# Added for the book (not part of the authors' package): a certificate, in Arb ball arithmetic, for the numbers used in
# the proof of Theorem ag:thm:explicit (Section ch:almost) and in the partial summation of Theorem sm:thm:strip (3)
# (Section ch:strip). Every comparison is directed: an Arb comparison x < y is True only if it holds for every point of
# the two balls, and every check below must print True.
#
# Theorem ag:thm:explicit (3) needs every theta in [0, T0/2], T0 = 3e12, to lie within gamma_1 of the ordinate of a zero
# on the critical line with ordinate below T0. Three ranges:
#   [0, gamma_1]      the ordinate gamma_1 itself;
#   [gamma_1, 2000]   the first 1800 zeros, from acb.zeta_zeros (FLINT's acb_dirichlet_zeta_zeros: the consecutive zeros
#                     from the first, isolated in order on the critical line with Gram's law and Rosser's rule in the
#                     ranges where they hold, as in v1min.py and Proposition bp:prop:v1min), reach above 2000, and
#                     consecutive ordinates differ by less than 6.8874, half of which is less than gamma_1;
#   [2000, T0/2]      the window (theta - H, theta + H] with H = 14 < gamma_1 holds at least
#                     n(theta) = (H/pi) log((theta - H)/(2 pi)) - 2 E(theta + H) zeros, where
#                     E(T) = 0.112 log T + 0.278 log log T + 2.51 + 0.2/T (T >= e) is Trudgian's bound for
#                     |N(T) - (T/2 pi) log(T/(2 pi e)) - 7/8|, as in the proof of Proposition bp:prop:v1min;
#                     n increases on [2000, oo), because n'(theta) >= (H/pi - 0.224 - 0.556/log(theta + H))/(theta + H)
#                     and the factor increases with theta, and n(2000) > 17; the zero found lies on the critical line
#                     by Platt and Trudgian, since theta + H < T0.
# The last step of (3): for 6e-10 < a < 1/2, a/(pi(a^2 + gamma_1^2)) > 6e-10/(pi(1/4 + gamma_1^2)) > 9.5e-13, and the
# bound of part (2) is 2 Sigma0/pi < 2(1.48e-12)/pi < 9.43e-13. Here Sigma0 is the sum of gamma^-2 over the zeros above T0;
# by partial summation against N(t) = (t/2 pi) log(t/(2 pi e)) + 7/8 + R(t), |R| <= E, it differs from the main term
# (log(T0/2 pi) + 1)/(2 pi T0) by at most E(T0)/T0^2 + 2 int_{T0}^oo E(t) t^-3 dt, which is bounded in closed form with
# log log t <= log log T0 + (log t - log T0)/log T0.
#
# Theorem sm:thm:strip (3): N(t) <= (t/2 pi) log t for t >= gamma_1 follows from Trudgian's bound once
# 7/8 + E(t) <= (t/2 pi) log(2 pi e) there. It holds at t = gamma_1, and for t >= gamma_1 the left side has slope at most
# (0.112 + 0.278/log gamma_1)/gamma_1, less than the slope log(2 pi e)/(2 pi) of the right side. Then
# sum_{gamma > T0} gamma^-2 <= (1 + log T0)/(pi T0) <= log T0/T0 for T0 >= 2, since (1 + log 2)/pi < log 2.
#
# Run on Hopper: bash misc/tools/hopper_run.sh -t 30 -g explicit_bounds.json riesz_program/code/explicit_bounds.py
# (about ten seconds); the output goes to riesz_program/data/explicit_bounds.json.
import json, time
import flint
from flint import arb, acb, ctx

t_start = time.time()
ctx.prec = 128
pi = arb.pi()
half = arb(1) / 2
out = {'python_flint': flint.__version__, 'prec_bits': ctx.prec}
checks = {}


def s(x, d=20):
    return x.str(d, radius=True)


def E(T):
    """Trudgian's bound for |N(T) - (T/2 pi) log(T/(2 pi e)) - 7/8|, valid for T >= e."""
    T = arb(T)
    return arb('0.112') * T.log() + arb('0.278') * T.log().log() + arb('2.51') + arb('0.2') / T


def n(theta, H):
    """Lower bound for the number of zeros with ordinate in (theta - H, theta + H], counted with multiplicity."""
    theta, H = arb(theta), arb(H)
    return H / pi * ((theta - H) / (2 * pi)).log() - 2 * E(theta + H)


# the first 1800 zeros
zs = []
for start in range(1, 1801, 200):
    batch = acb.zeta_zeros(start, 200)
    assert len(batch) == 200
    for z in batch:
        assert z.real == half and z.real.rad() == 0       # exactly on the critical line
        zs.append(z.imag)
g1 = zs[0]
gaps = [zs[i + 1] - zs[i] for i in range(len(zs) - 1)]
imax = max(range(len(gaps)), key=lambda i: float(gaps[i].mid()))
out['zeros'] = {'count': len(zs), 'gamma_1': s(g1, 25), 'gamma_1800': s(zs[-1], 20),
                'largest_gap': s(gaps[imax], 25), 'largest_gap_between': [s(zs[imax], 15), s(zs[imax + 1], 15)]}
checks['all_real_parts_exactly_one_half'] = True
checks['gaps_positive'] = all(g > 0 for g in gaps)
checks['gamma_1800_gt_2000'] = bool(zs[-1] > 2000)
checks['every_gap_lt_6.8874'] = all(g < arb('6.8874') for g in gaps)
checks['half_of_6.8874_lt_gamma_1'] = bool(arb('6.8874') / 2 < g1)

# the window count with H = 14 above height 2000
H = 14
T0 = arb(3) * arb(10) ** 12
n2000 = n(2000, H)
factor = arb(H) / pi - arb('0.224') - arb('0.556') / arb(2000 + H).log()
out['window'] = {'H': H, 'n_2000': s(n2000, 20), 'derivative_factor_at_2000': s(factor, 15)}
checks['H_lt_gamma_1'] = bool(arb(H) < g1)
checks['n_2000_gt_17'] = bool(n2000 > 17)
checks['derivative_factor_at_2000_positive'] = bool(factor > 0)
checks['E_increasing_beyond_2000_minus_H'] = bool(arb('0.112') * (2000 - H) > arb('0.2'))   # E'(T) > 0 if 0.112 T > 0.2
checks['T0_over_2_plus_H_lt_T0'] = bool(T0 / 2 + H < T0)

# Sigma0 and the last step of part (3)
L = T0.log()
LL = L.log()
main = ((T0 / (2 * pi)).log() + 1) / (2 * pi * T0)
int_E = (arb('0.112') * (2 * L + 1) / (4 * T0 ** 2) + arb('0.278') * (LL / (2 * T0 ** 2) + 1 / (4 * L * T0 ** 2))
         + arb('2.51') / (2 * T0 ** 2) + arb('0.2') / (3 * T0 ** 3))
err = E(T0) / T0 ** 2 + 2 * int_E
lb = arb('6e-10') / (pi * (arb(1) / 4 + g1 ** 2))
out['Sigma0'] = {'main_term': s(main, 15), 'error_bound': s(err, 10), 'two_times_1.48e-12_over_pi': s(2 * arb('1.48e-12') / pi, 15)}
out['part3_lower_bound_at_6e-10'] = s(lb, 15)
out['threshold_9.43e-13_pi_gamma1_sq'] = s(arb('9.43e-13') * pi * g1 ** 2, 15)
checks['Sigma0_main_term_matches_1.4797e-12'] = bool(abs(main - arb('1.4797e-12')) < arb('5e-17'))
checks['Sigma0_error_lt_2e-24'] = bool(err < arb('2e-24'))
checks['Sigma0_lt_1.48e-12'] = bool(main + err < arb('1.48e-12'))
checks['two_times_1.48e-12_over_pi_lt_9.43e-13'] = bool(2 * arb('1.48e-12') / pi < arb('9.43e-13'))
checks['lower_bound_at_6e-10_gt_9.5e-13'] = bool(lb > arb('9.5e-13'))

# Theorem sm:thm:strip (3): N(t) <= (t/2 pi) log t for t >= gamma_1, and the constant C = 4
lhs = arb(7) / 8 + E(g1)
rhs = g1 * (2 * pi * arb(1).exp()).log() / (2 * pi)
lhs_slope = (arb('0.112') + arb('0.278') / g1.log()) / g1
rhs_slope = (2 * pi * arb(1).exp()).log() / (2 * pi)
out['strip'] = {'seven_eighths_plus_E_at_gamma_1': s(lhs, 12), 'rhs_at_gamma_1': s(rhs, 12),
                'lhs_slope_bound': s(lhs_slope, 10), 'rhs_slope': s(rhs_slope, 10),
                'C4_bound_T0_3e12_u_le_1': s(8 * L / T0, 10)}
checks['strip_N_bound_at_gamma_1'] = bool(lhs < rhs)
checks['strip_slopes'] = bool(lhs_slope < rhs_slope)
checks['strip_constant_T0_ge_2'] = bool((1 + arb(2).log()) / pi < arb(2).log())

out['checks'] = checks
out['all_checks_true'] = all(checks.values())
out['elapsed_s'] = round(time.time() - t_start, 1)
print(json.dumps(out, indent=1))
json.dump(out, open('explicit_bounds.json', 'w'), indent=1)
if not out['all_checks_true']:
    raise SystemExit('a check failed: ' + ', '.join(k for k, v in checks.items() if not v))
