# Added for the book (not part of the authors' package): the Thorin density of the basepoint one clock T_1
# attains its minimum only at theta = 0, i.e. v_1(theta) = Re (xi'/xi)(1 + i theta) > v_1(0) = b_1 for theta != 0
# (Proposition bp:prop:v1min, Section ch:basepoint).
#
# Every inequality that the proof uses is checked here as a directed comparison of Arb balls (python-flint): an Arb
# comparison returns True only if it holds for every point of the balls. No ball is converted to a float before it is
# compared, and the decimal constants of the text enter as balls, arb('3.4437'), never as floats. The proof treats
# three ranges of theta >= 0 (v_1 is even):
#   [0, 13.5]     v_1 is strictly convex there, because the Cauchy kernel (1-beta)/((1-beta)^2 + x^2) is convex for
#                 |x| > (1-beta)/sqrt 3 and every ordinate has |gamma| >= gamma_1, with gamma_1 - 13.5 > 0.63 > 1/sqrt 3;
#   [13.5, 2000]  the first 1800 zeros come from acb.zeta_zeros, FLINT's consecutive indexed zeros from the first one,
#                 isolated on the critical line by Gram's law and Rosser's rule, so that none is missed; every real part
#                 is exactly 1/2 (an exact ball), gamma_1800 > 2000, and consecutive ordinates differ by at most
#                 gamma_2 - gamma_1 < 6.8874, so every theta lies within 3.4437 of an ordinate, and a zero on the line at
#                 distance x < 4.6259 contributes (1/2)/(1/4 + x^2) > b_1, more than 0.041 at x = 3.4437;
#   [2000, oo)    the explicit Riemann von Mangoldt formula of Hasanalizade, Shen and Wong [HSW22, (1.7), (2.7) and
#                 Corollary 1.4], |N(T) - (T/2 pi) log(T/(2 pi e)) - 7/8| <= E(T) for T >= e, with E(T) = 2.5167 + 1/(50T)
#                 for T <= X0 = 30 610 046 000 and E(T) = 0.1095 log T + 0.2042 log log T + 3.0305 + 1/(50T) < 0.3 log T
#                 beyond, gives at least n(theta) = (H/pi) log((theta - H)/(2 pi)) - E(theta - H) - E(theta + H) zeros in
#                 the window (theta - H, theta + H], H = 15. While theta + H <= X0, n increases (E decreases there) and
#                 n(2000) > 22.44; for theta + H > X0, n(theta) > (H/pi) log((theta - H)/(2 pi)) - 0.6 log(theta + H) > 90.
#                 So the window holds at least 23 zeros; the reflection rho -> 1 - conj(rho) permutes them, so their
#                 widths 1 - beta add up to half their number, and each of their terms in v_1 exceeds (1 - beta)/(1 + H^2);
#                 hence v_1(theta) > 23/452 > 0.05 > b_1 wherever the zeros lie. Backlund's classical bound
#                 E(T) = 0.137 log T + 0.443 log log T + 4.35 + 0.2/T (T >= 200), which increases, gives instead
#                 n(theta) = (H/pi) log((theta - H)/(2 pi)) - 2 E(theta + H), n(2000) > 14.89 and 15/452 > 0.033 > b_1.
#                 (The first version used Trudgian's bound 0.112 log T + 0.278 log log T + 2.510 + 0.2/T, with
#                 n(2000) > 19.62 and 20/452 > 0.044; its published proof is incomplete [BMOR21, HSW22].)
# The script also checks the remark on the Stirling term (its error is below 0.01/T for T >= 14, against 1/(50T) in E, and
# would move the counts by less than 2e-5), and it computes the numbers quoted after the proof: E T_1^circ = (1/2)(log xi)''(1), v_1''(0),
# values of v_1, the increase of v_1 on a grid of (0, 13.5], and the smallest value of v_1 on the grid of step 0.02 of
# [13.5, 2000]. These are checks of quoted numbers, not steps of the proof.
# Run from riesz_program/code (after riesz_program/setup.sh): python3 v1min.py. It writes v1min.json, which is
# riesz_program/data/v1min.json, prints it, and stops with an error if any check fails. It takes a few minutes.
import json, os, platform, time
import flint
from flint import arb, acb, acb_series, ctx

START = time.time()
ctx.prec = 128
half = arb(1) / 2
pi = arb.pi()


def s(x, d=20):
    """An Arb ball as a string, midpoint and radius."""
    return x.str(d, radius=True)


proof, quoted = {}, {}                      # name -> bool, each the value of a directed Arb comparison


def check(table, name, cond):
    assert isinstance(cond, bool), name
    table[name] = cond
    if not cond:
        print('CHECK FAILED:', name, flush=True)


out = {'python_flint': flint.__version__, 'precision_bits': ctx.prec,
       'slurm_job_id': os.environ.get('SLURM_JOB_ID'), 'host': platform.node()}

# b_1 = v_1(0) = sum 1/rho = (xi'/xi)(1)
b1 = 1 + arb.const_euler() / 2 - (4 * pi).log() / 2
out['b1'] = s(b1, 25)
check(quoted, 'b_1 = -(1/2) log pi + (1/2) psi(3/2) + gamma_E (the two forms overlap)',
      (-(pi.log()) / 2 + (arb(3) / 2).digamma() / 2 + arb.const_euler()).overlaps(b1))
check(quoted, '0.02309570896605 < b_1 < 0.02309570896615 (b_1 = 0.0230957089661 in the text)',
      arb('0.02309570896605') < b1 < arb('0.02309570896615'))


def v1(th):
    """Re (xi'/xi)(1 + i th) for real th > 0, as an Arb ball, from the arithmetic form 1/(1 + th^2) - (1/2) log pi
    + (1/2) Re psi((1 + i th)/2) + Re (zeta'/zeta)(1 + i th)."""
    c = acb_series([acb(1, th), 1], prec=2).zeta().coeffs()
    return 1 / (1 + th * th) - pi.log() / 2 + acb(half, th / 2).digamma().real / 2 + (c[1] / c[0]).real


# The first 1800 zeros, in order from the first (acb_dirichlet_zeta_zeros: consecutive indexed nontrivial zeros).
zs, on_line = [], True
for start in range(1, 1801, 200):
    batch = acb.zeta_zeros(start, 200)
    assert len(batch) == 200
    for z in batch:
        on_line = on_line and bool(z.real == half) and bool(z.real.rad() == 0)
        zs.append(z.imag)
gaps = [zs[i + 1] - zs[i] for i in range(len(zs) - 1)]
gap1 = gaps[0]
i2 = max(range(1, len(gaps)), key=lambda i: gaps[i].upper())        # the largest later gap, for the record only

# Range 1, [0, 13.5]: convexity
check(proof, 'gamma_1 > 14.134725, so every ordinate has |gamma| >= 14.134725', zs[0] > arb('14.134725'))
check(quoted, 'gamma_1 < 14.134726', zs[0] < arb('14.134726'))
check(proof, 'gamma_1 - 13.5 > 0.63', zs[0] - arb('13.5') > arb('0.63'))
check(proof, '0.63 > 1/sqrt 3', arb('0.63') > 1 / arb(3).sqrt())

# Range 2, [13.5, 2000]: the nearest zero
check(proof, 'every one of the 1800 zeros has real part exactly 1/2 (an exact ball)', on_line)
check(proof, 'the ordinates increase strictly (every gap > 0)', all(g > 0 for g in gaps))
check(proof, 'gamma_1800 > 2000', zs[-1] > 2000)
check(quoted, '2304.3 < gamma_1800 < 2304.4 (the zeros reach height 2304.3)', arb('2304.3') < zs[-1] < arb('2304.4'))
check(proof, 'every later gap < gamma_2 - gamma_1', all(g < gap1 for g in gaps[1:]))
check(proof, 'gamma_2 - gamma_1 < 6.8874', gap1 < arb('6.8874'))
check(quoted, 'gamma_2 - gamma_1 > 6.8873', gap1 > arb('6.8873'))
check(proof, 'every gap / 2 < 3.4437', all(g / 2 < arb('3.4437') for g in gaps))
check(proof, 'gamma_1 - 13.5 < 3.4437', zs[0] - arb('13.5') < arb('3.4437'))
check(quoted, 'gamma_1 - 13.5 < (gamma_2 - gamma_1)/2, so the largest distance from [13.5, 2000] to an ordinate is half the first gap',
      zs[0] - arb('13.5') < gap1 / 2)
xstar = (1 / (2 * b1) - arb(1) / 4).sqrt()                          # (1/2)/(1/4 + x^2) = b_1 at x = x_star
check(proof, 'x_star > 4.6259, so (1/2)/(1/4 + x^2) > b_1 for x < 4.6259', xstar > arb('4.6259'))
lb2 = half / (arb(1) / 4 + arb('3.4437') ** 2)
check(proof, '(1/2)/(1/4 + 3.4437^2) > 0.041', lb2 > arb('0.041'))
check(proof, '0.041 > b_1', arb('0.041') > b1)
out['zeros'] = {
    'zeros_used': len(zs), 'first_ordinate': s(zs[0], 30), 'second_ordinate': s(zs[1], 30), 'largest_ordinate': s(zs[-1], 20),
    'max_gap_gamma2_minus_gamma1': s(gap1, 30), 'max_gap_between': [s(zs[0], 20), s(zs[1], 20)],
    'largest_later_gap': s(gaps[i2], 15), 'largest_later_gap_between': [s(zs[i2], 15), s(zs[i2 + 1], 15)],
    'max_distance_to_nearest_zero_on_13.5_2000': s(gap1 / 2, 25), 'x_star': s(xstar, 20),
    'range2_lower_bound_at_distance_3.4437': s(lb2, 20), 'range2_lower_bound_at_half_the_first_gap': s(half / (arb(1) / 4 + (gap1 / 2) ** 2), 20)}

# Range 3, [2000, oo): the window count and the reflection
H = arb(15)
X0 = arb(30610046000)                     # |S(T)| <= 2.5167 for 0 <= T <= X0 [HSW22, (1.7)]


def E_hsw(T, st=None):
    """The bound E(T) for |N(T) - (T/2 pi) log(T/(2 pi e)) - 7/8| of [HSW22, (1.7), (2.7) and Corollary 1.4], T >= e, with
    the Stirling term st/T: st = 1/50 from (2.7) by default, st = 0.01 for the sharper bound of the remark below."""
    T = arb(T)
    c = arb(1) / 50 if st is None else arb(st)
    if T <= X0:
        return arb('2.5167') + c / T
    assert T > X0
    return arb('0.1095') * T.log() + arb('0.2042') * T.log().log() + arb('3.0305') + c / T


def n_hsw(theta, st=None):
    """n(theta) = (H/pi) log((theta - H)/(2 pi)) - E(theta - H) - E(theta + H), a lower bound for the number of zeros
    with ordinate in (theta - H, theta + H], counted with multiplicity."""
    theta = arb(theta)
    return H / pi * ((theta - H) / (2 * pi)).log() - E_hsw(theta - H, st) - E_hsw(theta + H, st)


def n_backlund(theta):
    """With Backlund's E(T) = 0.137 log T + 0.443 log log T + 4.35 + 0.2/T (T >= 200), which increases,
    n(theta) = (H/pi) log((theta - H)/(2 pi)) - 2 E(theta + H)."""
    theta = arb(theta)
    T = theta + H
    E = arb('0.137') * T.log() + arb('0.443') * T.log().log() + arb('4.35') + arb('0.2') / T
    return H / pi * ((theta - H) / (2 * pi)).log() - 2 * E


out['stirling_term_in_E'] = '1/50'
n2000 = n_hsw(2000)
EX0 = arb('0.1095') * X0.log() + arb('0.2042') * X0.log().log() + arb('3.0305') + 1 / (50 * X0)
slope = arb('0.3') - arb('0.1095') - arb('0.2042') / X0.log()
nX0 = H / pi * ((X0 - 2 * H) / (2 * pi)).log() - arb('0.6') * X0.log()   # the lower bound at theta + H = X0
tail = arb(23) / (2 * (1 + H * H))                                         # 23/452
check(proof, 'HSW22: n(2000) > 22.44', n2000 > arb('22.44'))
check(proof, 'HSW22: n(2000) > 22, so the window holds at least 23 zeros while theta + H <= X0, where n increases', n2000 > 22)
check(proof, 'HSW22: the second form of E at X0 is below 0.3 log X0', EX0 < arb('0.3') * X0.log())
check(proof, 'HSW22: 0.3 - 0.1095 - 0.2042/log X0 > 0, so E(T) < 0.3 log T for every T > X0', slope > 0)
check(proof, 'HSW22: 2.5167 + 1/(50 * 1985) < 0.3 log X0, which bounds E(theta - H) when theta - H <= X0 < theta + H',
      arb('2.5167') + arb(1) / (50 * 1985) < arb('0.3') * X0.log())
check(proof, 'HSW22: (H/pi) log((theta - H)/(2 pi)) - 0.6 log(theta + H) > 90 at theta + H = X0', nX0 > 90)
check(proof, 'HSW22: H/pi > 0.6, so that this lower bound increases with theta', H / pi > arb('0.6'))
check(proof, '1 + H^2 = 226', 1 + H * H == 226)
check(proof, '23/452 > 0.05', tail > arb('0.05'))
check(proof, '0.05 > b_1', arb('0.05') > b1)
nB = n_backlund(2000)
derivB = H / pi - arb('0.274') - arb('0.886') / arb(2015).log()
check(proof, 'Backlund: n(2000) > 14.89', nB > arb('14.89'))
check(proof, 'Backlund: n(2000) > 14, so the window holds at least 15 zeros', nB > 14)
check(proof, 'Backlund: H/pi - 0.274 - 0.886/log 2015 > 0, so n increases on [2000, oo)', derivB > 0)
check(proof, 'Backlund: E increases for T >= 200 (0.137 > 0.2/200)', arb('0.137') > arb('0.2') / 200)
check(proof, 'Backlund: 15/452 > 0.033 > b_1', arb(15) / 452 > arb('0.033') and arb('0.033') > b1)
out['range3'] = {
    'HSW22': {'E': '2.5167 + 1/(50T) for T <= X0 = 30610046000; 0.1095 log T + 0.2042 log log T + 3.0305 + 1/(50T) beyond',
              'n_2000': s(n2000, 25), 'n_2000_with_0.01_over_T': s(n_hsw(2000, '0.01'), 25),
              'second_form_of_E_at_X0': s(EX0, 15), 'three_tenths_log_X0': s(arb('0.3') * X0.log(), 15),
              'lower_bound_at_theta_plus_H_eq_X0': s(nX0, 15), 'zeros_in_window_at_least': 23,
              'lower_bound_for_v1_beyond_2000': '23/452', 'lower_bound_value': s(tail, 15)},
    'Backlund': {'constants': ['0.137', '0.443', '4.35'], 'n_2000': s(nB, 25), 'derivative_numerator_at_2000': s(derivB, 15),
                 'zeros_in_window_at_least': 15, 'lower_bound_for_v1_beyond_2000': '15/452'},
    'Trudgian_2014_not_used': 'the first version used 0.112 log T + 0.278 log log T + 2.510 + 0.2/T (n(2000) > 19.62, '
                              'at least 20 zeros), whose published proof is incomplete [BMOR21, HSW22]'}

# The remark on the Stirling term: the error of Stirling's formula in the smooth part theta(T)/pi + 1 of N(T), where
# theta(T) = Im log Gamma(1/4 + iT/2) - (T/2) log pi, is below 0.01/T for T >= 14. With z = 1/4 + iT/2, the bound
# 1/(90|z|^3) for the remainder of Stirling's series after 1/(12z) (DLMF 5.11(ii), Re z > 0) and the bounds
# x - x^3/3 <= arctan x <= x and y - y^2/2 <= log(1 + y) <= y give
# T |theta(T) - (T/2) log(T/(2 pi)) + T/2 + pi/8| <= 1/48 + (1/24 + 4/45)/T^2, which decreases in T, so its value at
# T = 14, divided by pi, bounds T times the error for every T >= 14. The grid evaluation only illustrates the bound.
stir_bound = (arb(1) / 48 + (arb(1) / 24 + arb(4) / 45) / 196) / pi
check(proof, '(1/48 + (1/24 + 4/45)/14^2)/pi < 0.01', stir_bound < arb('0.01'))
change = (arb(1) / 50 - arb('0.01')) * (arb(1) / 1985 + arb(1) / 2015)   # n with 0.01/T minus n with 1/(50T), at theta = 2000
check(proof, '(1/50 - 0.01)(1/1985 + 1/2015) < 2e-5, and the difference decreases with theta', change < arb('0.00002'))
check(proof, 'the two values of n(2000) differ by that amount', (n_hsw(2000, '0.01') - n2000 - change).contains(0))


def stirling_error(T):
    """T |theta(T) - (T/2) log(T/(2 pi)) + T/2 + pi/8| / pi as an Arb ball, at 200 bits because of the cancellation."""
    err = acb(arb(1) / 4, T / 2).lgamma().imag - T / 2 * arb.pi().log() - (T / 2 * (T / (2 * arb.pi())).log() - T / 2 - arb.pi() / 8)
    return abs(T * err / arb.pi())


ctx.prec = 200
grid = [arb(280 + k) / 20 for k in range(19720)] + [arb(10) ** (3 + arb(k) / 50) for k in range(301)]  # step 0.05 on [14, 1000), then up to 1e9
errs = [stirling_error(T) for T in grid]
iw = max(range(len(errs)), key=lambda i: errs[i].upper())
check(quoted, 'T times the Stirling error / pi is below the bound at every grid point of [14, 1e9]', all(e < stir_bound for e in errs))
ctx.prec = 128
out['stirling_check'] = {'error_times_T_bound_for_T_ge_14': s(stir_bound, 10), 'error_times_T_max_on_grid': s(errs[iw], 10),
                         'max_at_T': s(grid[iw], 10), 'one_over_48pi': s(1 / (48 * pi), 10), 'count_change_0.01_vs_1_over_50': s(change, 10)}

# The paragraph after the proof: E T_1^circ = (1/2)(log xi)''(1) = -(1/2) sum rho^-2 and v_1''(0) = -(log xi)'''(1),
# Taylor coefficients at x = 0 of log xi(1 + x) = log xi(-x) = log(1 + x) + (x/2) log pi + log Gamma(1 - x/2) + log(-zeta(-x)),
# a form that avoids the pole of zeta at 1.
L = (acb_series([1, 1], prec=4).log() + acb_series([0, pi.log() / 2], prec=4) + acb_series([1, -half], prec=4).lgamma()
     + (-acb_series([0, -1], prec=4).zeta()).log()).coeffs()
check(quoted, 'log xi(1) overlaps -log 2', L[0].real.overlaps(-arb(2).log()))
check(quoted, "(xi'/xi)(1) overlaps b_1", L[1].real.overlaps(b1))
ET1c, v1pp = L[2].real, -6 * L[3].real
# a second route to E T_1^circ: sum rho^-2 = 1 + gamma_E^2 + 2 gamma_1 - pi^2/8, with the Stieltjes constant gamma_1
check(quoted, 'E T_1^circ overlaps (pi^2/8 - 1 - gamma_E^2 - 2 gamma_1)/2', ET1c.overlaps((pi ** 2 / 8 - 1 - arb.const_euler() ** 2 - 2 * acb.stieltjes(1).real) / 2))
check(quoted, '0.02307715 < E T_1^circ < 0.02307725 (0.0230772 in the text)', arb('0.02307715') < ET1c < arb('0.02307725'))
check(quoted, "v_1''(0) > 0", v1pp > 0)
out['E_T1circ'] = s(ET1c, 20)
out['v1_second_derivative_at_0'] = s(v1pp, 20)

# Values of v_1, and the increase of v_1 on the grid of step 0.25 of (0, 13.5] (range 1 needs no computation)
out['v1_at'] = {t: s(v1(arb(t))) for t in ['1', '2', '5', '10', '13.5', '14.134725', '17.6', '100', '1000', '2000']}
d1, d5 = v1(arb(1)) - b1, v1(arb(5)) - b1
out['v1_minus_b1_at'] = {'1': s(d1), '5': s(d5)}
check(quoted, '1.115e-4 < v_1(1) - b_1 < 1.125e-4 (1.12e-4 in the text)', arb('0.0001115') < d1 < arb('0.0001125'))
check(quoted, '3.305e-3 < v_1(5) - b_1 < 3.315e-3 (3.31e-3 in the text)', arb('0.003305') < d5 < arb('0.003315'))
vals1 = [v1(arb(k) / 4) for k in range(1, 55)]
check(quoted, 'b_1 < v_1(0.25) < v_1(0.5) < ... < v_1(13.5)', vals1[0] > b1 and all(vals1[i] < vals1[i + 1] for i in range(len(vals1) - 1)))

# The smallest value of v_1 on the grid theta = 13.5 + 0.02 k of [13.5, 2000]
best, kbest, above = None, None, True
for k in range(99326):
    v = v1(arb(1350 + 2 * k) / 100)
    above = above and bool(v > b1)
    if best is None or v.upper() < best.upper():
        best, kbest = v, k
tbest = '%d.%02d' % divmod(1350 + 2 * kbest, 100)
check(quoted, 'every value of v_1 on the grid of step 0.02 of [13.5, 2000] exceeds b_1', above)
check(quoted, '0.113255 < min < 0.113265 on the grid, at theta = 17.50 (0.11326 at 17.5 in the text)',
      arb('0.113255') < best < arb('0.113265') and tbest == '17.50')
out['profile_min_on_13.5_2000'] = {'value': s(best), 'theta': tbest, 'grid_points': 99326}

out['checks_of_the_proof'] = proof
out['checks_of_quoted_numbers'] = quoted
# Midpoints under the key names of the earlier version of this script, which misc/repro/v1_collect.py reads as the
# reference of its independent floating point check with mpmath. No comparison of this script uses them.
out['for_v1_collect'] = 'x_star, max_gap, max_gap_between, max_distance_to_nearest_zero and profile_min_on_13.5_T0 are midpoints, for misc/repro/v1_collect.py'
out['x_star'] = float(xstar.mid())
out['max_gap'] = float(gap1.mid())
out['max_gap_between'] = [float(zs[0].mid()), float(zs[1].mid())]
out['max_distance_to_nearest_zero'] = float((gap1 / 2).mid())
out['profile_min_on_13.5_T0'] = {'value': float(best.mid()), 'theta': float(tbest)}
out['all_checks_pass'] = all(proof.values()) and all(quoted.values())
out['elapsed_s'] = round(time.time() - START, 1)
json.dump(out, open('v1min.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
assert out['all_checks_pass'], 'a check failed'
