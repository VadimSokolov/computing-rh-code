# pakes_sizebias.py: the size biased clock of Proposition pk:prop:sizebias (Section pk:sec:form, ch/pakes.tex).
# Added for the book in October 2026, for the merge of the first author's chapter on Pakes, the Lambert law and the Thorin
# condition (misc/reviews/nick-pakes/).
# The size biased clock is T + A with A independent of T and of density W(t)/E T; under RH, A = E/V with E standard
# exponential and P(V = gamma^2) = gamma^{-2}/E T. Computed here:
#   1. E T = sum_k 1/a_k = xi''(1/2)/(2 xi(1/2)) and Var T = sum_k 1/a_k^2 = (E T)^2 - 2 xi''''(1/2)/(24 xi(1/2)), from the
#      Taylor coefficients of xi at 1/2 (no zero enters), at 60 and at 80 digits; E A = Var T/E T;
#   2. the weights gamma_n^{-2}/E T for n <= 6 and their partial sums, from the 70 digit ordinates of zeros_1700.txt
#      (written by zeros_hp.py with Arb) and, independently, from mpmath's zetazero;
#   3. the sum of the weights over the first 1700 ordinates, and the same sum completed by the tail of the smooth count
#      N_0(t) = theta(t)/pi + 1 above gamma_1700, with the correction (N_0(gamma_1700) - 1700)/gamma_1700^2; under RH, which
#      holds below gamma_1700, the weights of all the zeros add up to one.
# Run on Hopper (seconds), from a folder containing zeros_1700.txt: python3 pakes_sizebias.py; writes pakes_sizebias.json there.
import json
import time

import mpmath as mp

t0 = time.time()
out = {}


def xi(s):
    return s * (s - 1) * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s) / 2


def moments(dps):
    with mp.workdps(dps):
        c = mp.taylor(xi, mp.mpf(1) / 2, 4)
        m0 = c[2] / c[0]
        var = m0 ** 2 - 2 * c[4] / c[0]
        return +c[0], m0, var


with mp.workdps(80):
    xi60, m0_60, var60 = moments(60)
    xi80, m0, var = moments(80)
    out['xi(1/2)'] = mp.nstr(xi80, 30)
    out['E T'] = mp.nstr(m0, 30)
    out['Var T'] = mp.nstr(var, 30)
    out['E T, 60 against 80 digits'] = mp.nstr(abs(m0 - m0_60), 3)
    out['Var T, 60 against 80 digits'] = mp.nstr(abs(var - var60), 3)
    out['E A = Var T / E T'] = mp.nstr(var / m0, 20)

    g = [mp.mpf(s) for s in open('zeros_1700.txt').read().split()]
    out['number of ordinates read'] = len(g)
    w = [x ** -2 / m0 for x in g]
    out['weights gamma_n^{-2}/E T, n <= 6'] = [mp.nstr(x, 10) for x in w[:6]]
    out['partial sums, n <= 6'] = [mp.nstr(mp.fsum(w[:n]), 10) for n in range(1, 7)]
    with mp.workdps(30):
        zz = [mp.im(mp.zetazero(n)) for n in range(1, 7)]
    out['max abs difference of gamma_1..gamma_6, zetazero against zeros_1700.txt'] = mp.nstr(
        max(abs(a - b) for a, b in zip(zz, g[:6])), 3)
    s1700 = mp.fsum(w)
    out['sum of the weights, n <= 1700'] = mp.nstr(s1700, 15)
    T = g[-1]
    smooth_tail = mp.quad(lambda t: mp.siegeltheta(t, 1) / mp.pi / t ** 2, [T, 2 * T, 10 * T, mp.inf])
    corr = (mp.siegeltheta(T) / mp.pi + 1 - len(g)) / T ** 2
    out['gamma_1700'] = mp.nstr(T, 20)
    out['N_0(gamma_1700) - 1700'] = mp.nstr(mp.siegeltheta(T) / mp.pi + 1 - len(g), 6)
    out['smooth tail above gamma_1700, in weight'] = mp.nstr(smooth_tail / m0, 10)
    out['sum n <= 1700 + tail + correction, minus 1'] = mp.nstr(s1700 + (smooth_tail + corr) / m0 - 1, 6)

out['seconds'] = round(time.time() - t0, 1)
json.dump(out, open('pakes_sizebias.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
