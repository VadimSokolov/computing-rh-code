# vdpair.py: the van Dantzig pair of Section ch:clock (ch/ch02.tex, Table tab:ch2:vd) and the remarks on it in Section
# ch:strip (ch/strip.tex). Added for the book in September 2026, from the discussion in the fourth version of the authors'
# note on [P26, Theorem A] (misc/reviews/nick-thm511/).
#  1. The two members Xi(theta)/Xi(0) = Re xi(1/2 + i theta)/xi(1/2) and xi(1/2)/xi(1/2 + theta), and their Hadamard
#     products prod (1 - theta^2/gamma^2) and prod (1 + theta^2/gamma^2)^(-1) over the 1700 zeros below 2197.26 in
#     zeros_1700.txt. The zeros above the last one enter through S2 = sum a_k^-1 = E T and S4 = sum a_k^-2 = Var T over all
#     zeros, read off the Taylor coefficients of log xi(1/2 + s), which need no hypothesis on the zeros: the tail factors are
#     exp(-theta^2 S2' - theta^4 S4'/2) and exp(-theta^2 S2' + theta^4 S4'/2), with S2' and S4' the sums over the zeros
#     above 2197.26, and the error is of order theta^6 sum gamma^-6 over those zeros.
#  2. The factor |1 + u^2/(gamma - i delta)^2|^(-2) that an off line quadruple 1/2 +- delta +- i gamma would contribute to
#     xi(1/2)/xi(1/2 + u): its inverse Fourier transform in closed form against quadrature at gamma = 1, delta = 0.4, and
#     the least value -(gamma^2 + delta^2) e^(-pi gamma/delta)/(4 gamma), at x = pi/delta.
#  3. Away from the centre: b_alpha = xi'(alpha)/xi(alpha) > 0, the corner of xi(alpha)/xi(alpha + |theta|) at 0, and
#     xi(alpha)/xi(alpha - theta) > 1 for small theta > 0, so that it is not a characteristic function.
# Run from any folder with zeros_1700.txt beside it: python3 vdpair.py; writes vdpair.json there.
import json

import mpmath as mp

mp.mp.dps = 40
pi, half = mp.pi, mp.mpf(1) / 2


def xi(s):   # xi(0) = xi(1) = 1/2, where the product below is 0 times the pole of zeta
    if s == 0 or s == 1:
        return half
    return s * (s - 1) / 2 * pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)


out = {}
xih = xi(half)
co = mp.taylor(lambda s: mp.log(xi(half + s)), 0, 4)
S2, S4 = co[2], -2 * co[4]
Z = [mp.mpf(x) for x in open('zeros_1700.txt').read().split()]
S2t, S4t = S2 - mp.fsum(g ** -2 for g in Z), S4 - mp.fsum(g ** -4 for g in Z)
out['sums'] = dict(S2=mp.nstr(S2, 12), S4=mp.nstr(S4, 12), S2_above_last=mp.nstr(S2t, 10), S4_above_last=mp.nstr(S4t, 10), zeros=len(Z), last=mp.nstr(Z[-1], 10))
rows = []
for th in [1, 2, 5, 10, 15, 20, 25, 30]:
    th = mp.mpf(th)
    phi = mp.re(xi(half + 1j * th)) / xih
    psi = xih / xi(half + th)
    below = sum(1 for g in Z if g < th)
    p1 = (-1) ** below * mp.exp(mp.fsum(mp.log(abs(1 - th * th / (g * g))) for g in Z) - th ** 2 * S2t - th ** 4 * S4t / 2)
    p2 = mp.exp(-mp.fsum(mp.log(1 + th * th / (g * g)) for g in Z) - th ** 2 * S2t + th ** 4 * S4t / 2)
    rows.append(dict(theta=int(th), zeros_below=below, phi=mp.nstr(phi, 10), psi=mp.nstr(psi, 10), phi_prod=mp.nstr(p1, 10), psi_prod=mp.nstr(p2, 10),
                     rel_diff=mp.nstr(max(abs(p1 / phi - 1), abs(p2 / psi - 1)), 2)))
    print(json.dumps(rows[-1]), flush=True)
out['table'] = rows


def dens(x, g, d):   # inverse Fourier transform of the off line factor
    x = abs(x)
    return (g * g + d * d) / (4 * g * d) * mp.exp(-g * x) * (g * mp.sin(d * x) + d * mp.cos(d * x))


def factor(u, g, d):
    return 1 / abs(1 + u * u / mp.mpc(g, -d) ** 2) ** 2


g, d = mp.mpf(1), mp.mpf('0.4')
chk = []
for x in [0, 2, 8, pi / d]:
    x = mp.mpf(x)
    q = mp.quad(lambda u: factor(u, g, d), [0, mp.inf]) / pi if x == 0 else mp.quadosc(lambda u: factor(u, g, d) * mp.cos(u * x), [0, mp.inf], omega=x) / pi
    chk.append(dict(x=mp.nstr(x, 6), closed_form=mp.nstr(dens(x, g, d), 10), quadrature=mp.nstr(q, 10)))
least = lambda g, d: -(g * g + d * d) * mp.exp(-pi * g / d) / (4 * g)
out['offline_factor'] = dict(check_at_gamma_one_delta_04=chk, least_at_gamma_one_delta_04=mp.nstr(least(g, d), 6),
                             least_at_first_zero_delta_04=mp.nstr(least(Z[0], d), 6), least_at_first_zero_delta_half=mp.nstr(least(Z[0], half), 6))
print(json.dumps(out['offline_factor']), flush=True)
out['away'] = [dict(alpha=float(a), b_alpha=mp.nstr(mp.diff(xi, a) / xi(a), 8), reciprocal_at_theta_025=mp.nstr(xi(a) / xi(a - mp.mpf('0.25')), 8))
               for a in [mp.mpf('0.75'), mp.mpf(1)]]
print(json.dumps(out['away']), flush=True)
json.dump(out, open('vdpair.json', 'w'), indent=1)
