# dbeta_window.py: the window count in the proof of Theorem thm:dbeta:minall(v) (Section sec:dbeta:comp, ch/basepoint.tex).
# Added for the book in October 2026, for the merge of the first author's note on the Dirichlet beta function in Thorin form
# (misc/reviews/nick-dbeta/), which leaves open the basepoints alpha in [1, 1/2 + 1/sqrt2) above the verified height.
# Claim: for f = xi(s, chi_4) and f = xi(s, chi_12), every alpha in [1, 1/2 + 1/sqrt2) and every theta > theta_0 =
# H_0 - (1+sqrt2)/2, with H_0 = 10^8/q the height to which Platt (Math. Comp. 85 (2016)) verified GRH for the modulus q,
# the closed window |gamma - theta| <= h holds n > V (alpha^2 + h^2) zeros, with V = (log f)''(1/2), h = 15 for chi_4 and
# h = 10 for chi_12; each of them contributes at least a/(alpha^2 + h^2), a = alpha - beta, to v_alpha(theta), and the a add
# up to n (alpha - 1/2), so v_alpha(theta) > (alpha - 1/2) V >= b_alpha.
# The count uses the explicit formula of Bennett, Martin, O'Bryant and Rechnitzer (Math. Comp. 90 (2021), Theorem 1.1):
# for a primitive character chi of conductor q > 1, T >= 5/7 and l = log(q(T+2)/(2 pi)) > 1.567,
# |N(T, chi) - (T/pi) log(qT/(2 pi e)) + chi(-1)/4| <= E(T) = 0.22737 l + 2 log(1+l) - 0.5, where N(T, chi) counts the
# zeros with |gamma| <= T. Since f has no real zeros and its zeros are symmetric about the real axis,
# n >= (N(theta+h) - N(theta-h))/2 >= (1/(2 pi)) int_{theta-h}^{theta+h} log(qt/(2 pi)) dt - E(theta+h)
#   >= low(theta) = (h/pi) log(q(theta-h)/(2 pi)) - E(theta+h),
# and low increases in theta, its derivative being at least (h/pi - 0.22737 - 2/(1+l(theta+h)))/(theta+h+2).
# Checked here in ball arithmetic (Arb through python-flint, 128 bits): the hypotheses of the formula at T = theta_0 - h
# (l increases with T, so they hold for every larger T); the positivity of h/pi - 0.22737 - 2/(1+l) at T = theta_0 + h, which
# gives the monotonicity for every theta >= theta_0; and low(theta_0) > V (alpha_top^2 + h^2), alpha_top = 1/2 + 1/sqrt2.
# V is computed from the Taylor series of log f at 1/2, as in dbeta_minimum.py.
# Run on Hopper (seconds): python3 dbeta_window.py; writes dbeta_window.json in the current folder.
import json

from flint import acb, acb_series, arb, ctx, dirichlet_char

ctx.prec = 128
PI = arb.pi()
SQRT2 = arb(2).sqrt()
CHI4 = dirichlet_char(4, 3)
CHI12 = dirichlet_char(12, 11)


def logser(name, x, n):
    s = acb_series([x, 1], prec=n)
    if name == 'chi4':
        return (s + 1) / 2 * acb((arb(4) / PI).log()) + ((s + 1) / 2).lgamma() + acb_series.dirichlet_l(s, CHI4).log()
    return s / 2 * acb((arb(12) / PI).log()) + (s / 2).lgamma() + acb_series.dirichlet_l(s, CHI12).log()


def ell(q, T):
    return (arb(q) * (T + 2) / (2 * PI)).log()


def E(q, T):
    lv = ell(q, T)
    return arb('0.22737') * lv + 2 * (1 + lv).log() - arb('0.5')


S = lambda x, n=10: x.str(n, radius=False)
out = {}
alpha_top = arb('0.5') + 1 / SQRT2
for name, q, h in [('chi4', 4, 15), ('chi12', 12, 10)]:
    V = 2 * logser(name, acb(0.5), 3).coeffs()[2].real
    H0 = arb(10) ** 8 / q
    theta0 = H0 - (1 + SQRT2) / 2
    Tlo, Thi = theta0 - h, theta0 + h
    hyp = bool(Tlo >= arb(5) / 7) and bool(ell(q, Tlo) > arb('1.567'))
    deriv = arb(h) / PI - arb('0.22737') - 2 / (1 + ell(q, Thi))
    mainterm = arb(h) / PI * (arb(q) * Tlo / (2 * PI)).log()
    low = mainterm - E(q, Thi)
    need = V * (alpha_top ** 2 + h * h)
    out[name] = {'q': q, 'h': h, 'V': S(V, 12), 'H_0 = 10^8/q': S(H0, 12), 'theta_0 = H_0 - (1+sqrt2)/2': S(theta0, 12),
                 'hypotheses T >= 5/7 and l > 1.567 at theta_0 - h': hyp, 'l at theta_0 - h': S(ell(q, Tlo), 10),
                 'log(q(theta_0-h)/(2 pi))': S((arb(q) * Tlo / (2 * PI)).log(), 10),
                 'E(theta_0+h)': S(E(q, Thi), 10), '(h/pi) log(q(theta_0-h)/(2 pi))': S(mainterm, 10),
                 'low(theta_0)': S(low, 10), 'alpha_top^2': S(alpha_top ** 2, 10),
                 'need V(alpha_top^2 + h^2)': S(need, 10), 'low(theta_0) > need': bool(low > need),
                 'h/pi - 0.22737 - 2/(1+l) at theta_0 + h': S(deriv, 8), 'monotone': bool(deriv > 0)}
    print(name, json.dumps(out[name], indent=1), flush=True)
    assert hyp and bool(low > need) and bool(deriv > 0)
json.dump(out, open('dbeta_window.json', 'w'), indent=1)
