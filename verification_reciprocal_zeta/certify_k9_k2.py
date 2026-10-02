"""Certified enclosures for Corollary rz:cor:Tfalse and Theorem rz:thm:class (K9 and K2 of the review).

K9: Re kappa(rho_213), with kappa(rho) = rho xi(3 - rho) / (2 pi xi'(rho)) (Proposition rz:prop:firstorder)
    and the Hardy function formula of Theorem rz:thm:kappaexact, both evaluated; Z'(gamma_213), so that
    rho_213 is simple; kappa at gamma_289 and gamma_379, and the sign of Re kappa at the first 460 zeros.
K2: an interval [T1, T2] below gamma_289 with N(T) = 288 and 1 + S(T) < 0, S(T) = N(T) - 1 - theta(T)/pi
    (the numerical remark in the proof of Theorem rz:thm:class), S(gamma_289^-), and the zeros among the
    first 460 with S(gamma^-) < -1.

Everything is python-flint 0.9 (Arb) ball arithmetic. The zeros come from acb.zeta_zero, which isolates each
zero of Hardy's Z function by a sign change and checks the count by Turing's method, so the n-th zero is
certified to lie on the critical line with the right index; N(T) comes from arb.zeta_nzeros. Writes
certify_k9_k2.json beside itself and prints the same JSON.
"""
import json
from pathlib import Path

from flint import arb, acb, ctx, arb_series, acb_series

ctx.prec = 256
PI = arb.pi()
out = {}


def s(x, d=25):
    return x.str(d, radius=True)


def theta(t):
    """Riemann Siegel theta via Arb's continuous log gamma."""
    return (acb(arb(1) / 4, t / 2).lgamma()).imag - t / 2 * PI.log()


def theta_series(t):
    return arb_series.riemann_siegel_theta(arb_series([t, 1], prec=ctx.prec))


def zprime(t):
    """Z'(t) as an enclosure over the ball t, from the Taylor series of Hardy's Z."""
    ser = arb_series.riemann_siegel_z(arb_series([t, 1], prec=ctx.prec))
    return ser.coeffs()[0], ser.coeffs()[1]


def xi_series(s0):
    """xi(s0 + x) as a power series to order x^1, over the ball s0."""
    S = acb_series([s0, 1], prec=ctx.prec)
    half = acb(1) / 2
    xs = half * S * (S - 1) * (-(S / 2) * acb(PI).log()).exp() * (S / 2).gamma() * S.zeta()
    return xs.coeffs()[0], xs.coeffs()[1]


def xi(s0):
    s0 = acb(s0)
    return s0 * (s0 - 1) / 2 * (-(s0 / 2) * acb(PI).log()).exp() * (s0 / 2).gamma() * s0.zeta()


def kappa_def(n):
    rho = acb.zeta_zero(n)
    xi0, xi1 = xi_series(rho)
    k = rho * xi(3 - rho) / (2 * PI * xi1)
    return rho, xi0, xi1, k


def kappa_hardy(g):
    th = theta(g)
    a = acb(0, th).exp() * acb(arb(5) / 2, g).zeta()
    z0, z1 = zprime(g)
    c = g * g - arb(15) / 4
    re = (c * a.imag - 4 * g * a.real) / (4 * PI ** 2 * z1)
    im = (4 * g * a.imag + c * a.real) / (4 * PI ** 2 * z1)
    return re, im, z0, z1, a, th


# ---------------- K9 ----------------
k9 = {}
rho, xi0, xi1, kd = kappa_def(213)
g = rho.imag
re_h, im_h, z0, z1, a, th = kappa_hardy(g)
k9["gamma_213"] = s(g, 40)
k9["gamma_213_radius"] = str(g.rad())
k9["Re_rho_213"] = s(rho.real, 10)
k9["neighbours"] = {"gamma_212": s(acb.zeta_zero(212).imag, 20), "gamma_214": s(acb.zeta_zero(214).imag, 20)}
k9["N(415)"] = s(arb(415).zeta_nzeros(), 5)
k9["N(415.4)"] = s(arb("415.4").zeta_nzeros(), 5)
k9["N(415.5)"] = s(arb("415.5").zeta_nzeros(), 5)
k9["Z(gamma_213) over the ball"] = s(z0, 5)
k9["Z'(gamma_213)"] = s(z1, 30)
k9["theta(gamma_213)"] = s(th, 30)
k9["theta via arb_series (check)"] = s(theta_series(g).coeffs()[0], 30)
k9["a(gamma_213) = e^{i theta} zeta(5/2 + i gamma)"] = [s(a.real, 25), s(a.imag, 25)]
k9["zeta(5/2+i gamma)"] = [s(acb(arb(5) / 2, g).zeta().real, 25), s(acb(arb(5) / 2, g).zeta().imag, 25)]
k9["xi(rho) over the ball"] = s(abs(xi0), 5)
k9["xi'(rho_213)"] = [s(xi1.real, 25), s(xi1.imag, 25)]
k9["kappa by definition"] = [s(kd.real, 30), s(kd.imag, 30)]
k9["kappa by Theorem rz:thm:kappaexact"] = [s(re_h, 30), s(im_h, 30)]
k9["|kappa|"] = s(abs(kd), 20)
k9["Re kappa lower end (definition)"] = str(kd.real.lower())
k9["Re kappa lower end (Hardy formula)"] = str(re_h.lower())
k9["Re kappa > 233.95 certified"] = bool(kd.real > arb("233.95")) and bool(re_h > arb("233.95"))
k9["Re kappa < 233.96 certified"] = bool(kd.real < arb("233.96"))
k9["Z'(gamma_213) > 1.9 certified (so the zero is simple)"] = bool(z1 > arb("1.9"))
out["K9"] = k9

# the other two zeros with Re kappa > 0 and the first four zeros, as printed in the book
others = {}
for n in (1, 2, 3, 4, 289, 379):
    r_, x0_, x1_, k_ = kappa_def(n)
    reh, imh, _, _, _, _ = kappa_hardy(r_.imag)
    others[str(n)] = {"gamma": s(r_.imag, 20), "kappa_def": [s(k_.real, 15), s(k_.imag, 15)],
                      "kappa_hardy": [s(reh, 15), s(imh, 15)]}
out["kappa_other_zeros"] = others

# certified sign scan of Re kappa over the first 460 zeros
pos, undecided, agree = [], [], True
zs = acb.zeta_zeros(1, 460)
for i, r_ in enumerate(zs, 1):
    reh, imh, z0_, z1_, a_, th_ = kappa_hardy(r_.imag)
    if reh > 0:
        pos.append(i)
    elif not (reh < 0):
        undecided.append(i)
    # sign of Im a / Z'
    q = a_.imag / z1_
    if not ((q > 0 and reh > 0) or (q < 0 and reh < 0)):
        agree = False
out["Re_kappa_scan_460"] = {"positive_at": pos, "undecided": undecided,
                            "sign agrees with sign of Im a / Z' at every zero": agree}

# ---------------- K2 ----------------
k2 = {}
g287 = acb.zeta_zero(287).imag
g288 = acb.zeta_zero(288).imag
g289 = acb.zeta_zero(289).imag
gram288 = arb.gram_point(288)
k2["gamma_287"] = s(g287, 25)
k2["gamma_288"] = s(g288, 25)
k2["gamma_289"] = s(g289, 25)
k2["gram_288 (theta = 288 pi)"] = s(gram288, 25)
k2["S(gamma_289^-) = 287 - theta(gamma_289)/pi"] = s(287 - theta(g289) / PI, 25)
k2["S(gamma_289^+) = 288 - theta(gamma_289)/pi"] = s(288 - theta(g289) / PI, 25)
T1, T2 = arb("527.71"), arb("527.90")
chk = {}
chk["gamma_288 < gram_288"] = bool(g288 < gram288)
chk["gram_288 < T1"] = bool(gram288 < T1)
chk["T1 < T2 < gamma_289"] = bool(T1 < T2) and bool(T2 < g289)
chk["gamma_288 < T1"] = bool(g288 < T1)
chk["N(T1) = 288 (zeta_nzeros)"] = s(T1.zeta_nzeros(), 5)
chk["N(T2) = 288 (zeta_nzeros)"] = s(T2.zeta_nzeros(), 5)
th1, th2 = theta(T1), theta(T2)
chk["theta(T1)/pi"] = s(th1 / PI, 25)
chk["theta(T2)/pi"] = s(th2 / PI, 25)
chk["theta(T1) > 288 pi"] = bool(th1 > 288 * PI)
# theta' > 0 on [T1, T2]: evaluate theta' over the whole interval as a ball
Tball = arb((T1 + T2) / 2, (T2 - T1) / 2)
thp = theta_series(Tball).coeffs()[1]
chk["theta' over [T1,T2]"] = s(thp, 10)
chk["theta' > 0 on [T1,T2]"] = bool(thp > 0)
chk["1+S(T1) = 288 - theta(T1)/pi"] = s(288 - th1 / PI, 25)
chk["1+S(T2) = 288 - theta(T2)/pi"] = s(288 - th2 / PI, 25)
chk["S on [T1,T2] lies in"] = [s(287 - th2 / PI, 12), s(287 - th1 / PI, 12)]
chk["backlund_s(T1)"] = s(T1.backlund_s(), 20)
chk["backlund_s(T2)"] = s(T2.backlund_s(), 20)
chk["backlund_s(T2) < -1"] = bool(T2.backlund_s() < -1)
k2["interval"] = ["527.71", "527.90"]
k2["checks"] = chk

# S(gamma_k^-) scan, k <= 460: jumps, the zeros with S(gamma^-) < -1, the most negative value
sm = []
for i, r_ in enumerate(zs, 1):
    t_ = r_.imag
    val = (i - 2) - theta(t_) / PI
    sm.append((i, val))
below = [(i, s(v, 8), s(v + arb(1) / 2, 6)) for i, v in sm if v < -1]
undec = [i for i, v in sm if not (v < -1) and not (v > -1)]
most = min(sm, key=lambda iv: iv[1].mid())
k2["zeros k<=460 with S(gamma_k^-) < -1 (k, S^-, midpoint)"] = below
k2["undecided against -1"] = undec
k2["most negative S(gamma_k^-) for k<=460"] = [most[0], s(most[1], 12)]
k2["S jumps at 127, 213, 289, 379 (S^-, S^+, midpoint)"] = {
    str(k): [s(sm[k - 1][1], 8), s(sm[k - 1][1] + 1, 8), s(sm[k - 1][1] + arb(1) / 2, 8)] for k in (127, 213, 289, 379)}
out["K2"] = k2

print(json.dumps(out, indent=1))
with open(Path(__file__).resolve().parent / "certify_k9_k2.json", "w") as f:
    json.dump(out, f, indent=1)
