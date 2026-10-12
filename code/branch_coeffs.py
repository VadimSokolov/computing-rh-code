# branch_coeffs.py: the candidate branching mechanism of the centre clock, Theorem pk:thm:nobranch in Section pk:sec:ladder
# (ch/pakes2.tex). Added for the book in October 2026, for the merge of the first author's chapter on Pakes, the Lambert
# law and the Thorin condition (misc/reviews/nick-pakes/, triage-C.md). If 1/G were the Laplace transform of the martingale
# limit of a supercritical continuous state branching process with finite Malthusian parameter mu, Pakes's relation
# mu theta omega'(theta) = C(omega(theta)) (Axioms 15 (2026), article 574, eq. (68)), with omega = log G up to a scale,
# would make the cumulant of the generating process a positive multiple of g = Psi_Z o (log G)^{-1}, where
# Psi_Z(u) = u (log G)'(u) is the Laplace exponent of the driving process, and would force (-1)^{k-1} g_k >= 0 for k >= 2
# in g(lambda) = sum_k g_k lambda^k. Certified here, in ball arithmetic (python-flint, Arb):
#  1. the Taylor coefficients of log G(u) = log(xi(1/2 + sqrt u)/xi(1/2)) = sum_n (-1)^{n+1} m_{n-1} u^n/n, with m_n the
#     Stieltjes moments of Section ch:heat, at 3200 bits from the power series of xi(1/2 + z) (series of zeta and gamma;
#     no zeros, no quadrature);
#  2. g_k for k <= 60, by reverting the series of log G and composing with that of Psi_Z, with every sign certified, and
#     by a second route, g = H/H' with H = (log G)^{-1}, which agrees for k <= 59 (it loses the last coefficient); the blocks
#     of k with (-1)^{k-1} g_k < 0;
#  3. the critical points of Xi in (gamma_1, gamma_2), (gamma_2, gamma_3) and (gamma_3, gamma_4), each by a certified sign
#     change of Xi', and the values of G = Xi/Xi(0) there; for the first, s* = -t*^2 and lambda* = log|G(s*)| + i pi, its
#     modulus and argument, 1/|lambda*| and the period 2 pi/(pi - arg lambda*) of the sign pattern.
# Outside the certificate: the ratios |g_{k+1}/g_k| and the root test |g_k|^{1/k}; the moments m_{k-1}, 2 <= k <= 12,
# against the first 1700 zeros (data/zeros_1700.txt) with a tail estimate; C'(s) = (log G)'(s) at points of (-gamma_1^2, 0);
# and a numerical continuation of (log G)^{-1} along the segment from 0 to lambda* in 3000 steps (Newton steps on recentred
# balls, so not a certificate), which approaches s* like the square root of the distance to lambda*.
# Run on Hopper (seconds), from a folder containing zeros_1700.txt: python3 branch_coeffs.py; writes branch_coeffs.json there.
import json
import time

from flint import arb, acb, arb_series, acb_series, ctx

T0 = time.time()
K = 60
ctx.prec = 3200
ctx.cap = 2 * K + 2
out, holds = {'K': K, 'prec_bits': 3200}, []


def s(x, n=15):
    return x.str(n, radius=False) if isinstance(x, arb) else x


def check(name, cond, **vals):
    holds.append(bool(cond))
    out[name] = dict(holds=bool(cond), **vals)
    print(name, bool(cond), round(time.time() - T0, 1), flush=True)


# 1. The series of xi(1/2 + z) and of log G.
z = arb_series([arb(1) / 2, 1])
xs = z * (z - 1) / 2 * (-(z / 2) * arb.pi().log()).exp() * (z / 2).gamma() * z.zeta()
c = xs.coeffs()
lc = (xs / c[0]).log().coeffs()
kap = [arb(0)] + [lc[2 * k] for k in range(1, K + 1)]  # log G(u) = sum_k kap_k u^k
mom = [(-1) ** k * (k + 1) * kap[k + 1] for k in range(K)]  # m_n = (-1)^n (n+1) kap_{n+1}
check('series', all(c[k].contains(0) for k in range(1, 2 * K + 2, 2)) and all(mom[n] > 0 for n in range(K)),
      xi_half=c[0].str(30), m_0_to_5=[mom[n].str(20) for n in range(6)],
      note='odd coefficients of xi(1/2 + z) contain 0 (xi is even about 1/2); every m_n, n < 60, is positive')

# 2. g = Psi_Z o (log G)^{-1}, two routes.
ctx.cap = K + 1
Cs = arb_series([arb(0)] + kap[1:])
Ps = arb_series([arb(0)] + [k * kap[k] for k in range(1, K + 1)])
H = Cs.reversion()
g1 = Ps(H).coeffs()
g1 += [arb(0)] * (K + 1 - len(g1))
g2 = (H / H.derivative()).coeffs()
g2 += [arb(0)] * (K + 1 - len(g2))
sign = {k: (1 if g1[k] > 0 else (-1 if g1[k] < 0 else 0)) for k in range(1, K + 1)}
wrong = [k for k in range(2, K + 1) if sign[k] != (-1) ** (k - 1)]
blocks = []
for k in wrong:
    if blocks and blocks[-1][1] == k - 1:
        blocks[-1][1] = k
    else:
        blocks.append([k, k])
check('coefficients', all(sign[k] != 0 for k in range(1, K + 1)) and all(g1[k].overlaps(g2[k]) for k in range(1, K))
      and g1[16] > 0 and g1[1].overlaps(arb(1)),
      g={k: g1[k].str(15) for k in range(1, K + 1)}, wrong_sign_blocks=blocks, g16=g1[16].str(15), g22=g1[22].str(15),
      routes_agree_for_k_le=max(k for k in range(1, K + 1) if all(g1[i].overlaps(g2[i]) for i in range(1, k + 1))))
out['ratios_abs_g_next_over_g'] = {k: round(abs(float((g1[k + 1] / g1[k]).mid())), 4) for k in range(1, K)}
out['root_test'] = {k: round(float(abs(g1[k]).mid()) ** (1.0 / k), 5) for k in (10, 20, 30, 40, 50, 60)}
rat = [abs(float((g1[k + 1] / g1[k]).mid())) for k in range(1, K)]
out['ratio_range_k_2_to_9'] = [round(min(rat[1:9]), 4), round(max(rat[1:9]), 4)]
out['ratio_range_k_10_to_59'] = [round(min(rat[9:]), 4), round(max(rat[9:]), 4)]

# 3. Critical points of Xi between the first four zeros.
ctx.prec = 256
zs = [arb(l.strip(), '1e-66') for l in open('zeros_1700.txt') if l.strip()]


def Xi_d(t):  # Xi(t) and Xi'(t), Xi(t) = xi(1/2 + it)
    ctx.cap = 2
    w = acb_series([acb(arb(1) / 2, t), acb(0, 1)])
    cc = (w * (w - 1) / 2 * (-(w / 2) * acb.pi().log()).exp() * (w / 2).gamma() * w.zeta()).coeffs()
    return cc[0].real, cc[1].real


Xi0 = Xi_d(arb(0))[0]
crit, ok = [], True
for i in range(3):
    lo, hi = arb(zs[i].mid()) + arb('1e-8'), arb(zs[i + 1].mid()) - arb('1e-8')
    dlo = Xi_d(lo)[1]
    for _ in range(170):
        mid = (lo + hi) / 2
        dm = Xi_d(mid)[1]
        if not (dm > 0 or dm < 0):
            break
        if (dm > 0) == (dlo > 0):
            lo, dlo = mid, dm
        else:
            hi = mid
    dhi = Xi_d(hi)[1]
    ok = ok and bool((dlo > 0 and dhi < 0) or (dlo < 0 and dhi > 0))
    tb = lo.union(hi)
    Gv = Xi_d(tb)[0] / Xi0
    crit.append(dict(t=tb.str(15), s=(-(tb * tb)).str(12), G=Gv.str(12), log_abs_G=abs(Gv).log().str(12)))
    if i == 0:
        tstar, Gstar = tb, Gv
lam = acb(abs(Gstar).log(), arb.pi())
check('critical_points', ok and Gstar < 0, critical_points=crit,
      lambda_star=dict(re=lam.real.str(12), im='pi', abs=abs(lam).str(10), arg=lam.arg().str(10),
                       inverse_abs=(1 / abs(lam)).str(8), period=(2 * arb.pi() / (arb.pi() - lam.arg())).str(8)))
out['all_hold'] = all(holds)

# Outside the certificate. (i) The moments against the zeros: m_{k-1} = sum_j a_j^{-k}; the 1700 zeros below height 2197.3
# lie on the line, and the zeros above add at most about (k/pi) G^{1-2k} (log G/(2k - 1) + 1/(2k - 1)^2), G = 2197.3.
ctx.prec = 300
Gam = arb('2197.3')
cz = {}
for k in range(2, 13):
    S = sum((x * x) ** (-k) for x in zs)
    tail = arb(k) / arb.pi() * Gam ** (1 - 2 * k) * (Gam.log() / (2 * k - 1) + arb(1) / (2 * k - 1) ** 2)
    cz[k] = dict(m_from_xi=mom[k - 1].str(20), sum_1700=S.str(20), tail_estimate=tail.str(4),
                 consistent=bool(mom[k - 1].overlaps(S + tail * arb(0.5, 0.5))))
out['moments_vs_zeros'] = cz

# (ii) C'(s) = -Xi'(t)/(2t Xi(t)) at s = -t^2 on (-gamma_1^2, 0).
ctx.prec = 160
out['C_prime_on_negative_axis'] = {sv: (-Xi_d(arb(sv[1:]).sqrt())[1] / (2 * arb(sv[1:]).sqrt() * Xi_d(arb(sv[1:]).sqrt())[0])
                                        ).str(10) for sv in ('-1', '-50', '-100', '-150', '-190', '-199', '-199.7')}


# (iii) Continuation of (log G)^{-1} from 0 along the segment to lambda*: Newton steps for log G(s) = tau lambda*, with
# log G followed continuously, log G(s_new) = log G(s_old) + log(G(s_new)/G(s_old)).
def G_and_Cp(sv):  # G(s) and C'(s) = (xi'/xi)(1/2 + sqrt s)/(2 sqrt s)
    w0 = acb(arb(1) / 2) + sv.sqrt()
    ctx.cap = 2
    w = acb_series([w0, 1])
    cc = (w * (w - 1) / 2 * (-(w / 2) * acb.pi().log()).exp() * (w / 2).gamma() * w.zeta()).coeffs()
    return cc[0] / XIH, cc[1] / cc[0] / (2 * sv.sqrt())


ctx.prec = 160
XIH = acb(c[0])
lam_s = acb(arb(lam.real.mid()), arb.pi())
sst = -(tstar * tstar)
N = 3000
sv = lam_s / (N * kap[1].mid())
Gv, Cp = G_and_Cp(sv)
Cv = Gv.log()
for _ in range(30):
    sv = (sv - (Cv - lam_s / N) / Cp).mid()
    Gv, Cp = G_and_Cp(sv)
    Cv = Gv.log()
path, minG = [], 1e9
for i in range(2, N):
    lam_i = lam_s * i / N
    s_old, G_old, C_old = sv, Gv, Cv
    sn = (s_old + (lam_i - C_old) / Cp).mid()
    for _ in range(12):
        Gn, Cpn = G_and_Cp(sn)
        sn = (sn - (C_old + (Gn / G_old).log() - lam_i) / Cpn).mid()
    Gv, Cp = G_and_Cp(sn)
    Cv = C_old + (Gv / G_old).log()
    sv = sn
    minG = min(minG, float(abs(Gv).mid()))
    if i % 300 == 0 or i >= N - 3:
        d = abs(sv - acb(sst))
        path.append(dict(tau=round(i / N, 5), s=[round(float(sv.real.mid()), 6), round(float(sv.imag.mid()), 6)],
                         residual=float(abs(Cv - lam_i).mid()), abs_C_prime=float(abs(Cp).mid()),
                         dist_to_s_star=float(d.mid()), dist_over_sqrt_1_minus_tau=float((d / arb(1 - i / N).sqrt()).mid())))
out['continuation_to_lambda_star'] = dict(steps=N, path=path, min_abs_G_on_path=minG)
out['seconds'] = round(time.time() - T0, 1)
json.dump(out, open('branch_coeffs.json', 'w'), indent=1)
print(json.dumps({k: out[k] for k in ('coefficients', 'critical_points', 'root_test', 'ratio_range_k_2_to_9',
                                      'ratio_range_k_10_to_59', 'all_hold', 'seconds')}, indent=1)[:6000])
print(json.dumps(out['continuation_to_lambda_star'], indent=1))
