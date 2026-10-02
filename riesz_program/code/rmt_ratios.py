# Added for the book (not part of the authors' package): the exact moments E|Z_N|^{2k} of the characteristic
# polynomial of an N x N Haar unitary matrix and the ratios E|Z_N|^{2k}/N^{k^2} for k = 1, 2, 3 and
# N = 10, 20, 10^2, 10^3, 10^4, which Chapter ch:almost quotes against the Barnes constants G(1+k)^2/G(1+2k) in the
# sentence after Table ag:tab:ks (Section ag:sec:moments).
# The text calls them "the exact values", so they are the field 'ratio' = exact/N^{k^2} of killing2.py carried to k = 3
# and to N = 10^2, 10^3, 10^4; no simulation is needed. Five independent evaluations, which must agree:
#   (1) the formula of killing2.py, prod_{j<=N} Gamma(j)Gamma(j+2k)/Gamma(j+k)^2, with mpmath at 30 digits as there;
#   (2) the same product in exact rational arithmetic (for integer k each factor is prod_{i<k} (j+k+i)/(j+i)),
#       and its telescoped closed form prod_{i<k} binom(N+k+i,k)/binom(k+i,k);
#   (3) the beta decomposition of Bourgade, Hughes, Nikeghbali and Yor: over a uniform angle w,
#       E|1+e^{iw}sqrt(B_{1,j-1})|^{2k} = sum_m binom(k,m)^2 E B^m with E B_{1,j-1}^m = m!/(j(j+1)...(j+m-1)),
#       multiplied over j = 1..N (exact);
#   (4) Heine's identity: E|Z_N|^{2k} is the N x N Toeplitz determinant of the symbol |1-e^{it}|^{2k}, whose Fourier
#       coefficients are (-1)^m binom(2k,k+m), |m| <= k; exact elimination of the band matrix gives every N <= 10^4;
#   (5) Arb balls of exp(sum_{j<=N} [lgamma(j)+lgamma(j+2k)-2 lgamma(j+k)])/N^{k^2}, a rigorous enclosure.
# The Barnes constants come exactly from G(n) = 0!1!...(n-2)!, and from mpmath and Arb. The closed form gives
# E|Z_N|^{2k}/N^{k^2} = G(1+k)^2/G(1+2k) (1 + k^3/N + O(N^{-2})), recorded as (ratio/Barnes - 1) N/k^3.
# Reads killing.json (if present) to compare its N = 10, 20 entries; writes rmt_ratios.json.
import json, os, socket, time
from decimal import Decimal, getcontext, ROUND_HALF_EVEN
from fractions import Fraction
from math import comb, factorial
import mpmath as mp
import flint
from flint import arb, acb, ctx

ctx.prec = 256
getcontext().prec = 80
KS = [1, 2, 3]
NS = [10, 20, 100, 1000, 10000]
NMAX = max(NS)
# the values quoted in ch/almost.tex (mantissa, power of ten)
QUOTED = {2: {10: ('0.1716', 0), 100: ('0.0902', 0), 1000: ('0.0840', 0), 10000: ('0.0834', 0)},
          3: {10: ('1.18', -3), 100: ('1.51', -4), 1000: ('1.189', -4), 10000: ('1.1605', -4)}}
QUOTED_BARNES = {2: ('1/12', None), 3: ('1.1574', -4)}
t_start = time.time()


def dec(x):
    return Decimal(x.numerator) / Decimal(x.denominator)


def rounded_like(x, q):
    """x (a Fraction) written with as many decimals as the quoted mantissa q[0] after scaling by 10^-q[1]."""
    m, e = q
    d = len(m.split('.')[1]) if '.' in m else 0
    return str((dec(x) / Decimal(10) ** e).quantize(Decimal(1).scaleb(-d), rounding=ROUND_HALF_EVEN))


def closed(N, k):
    r = Fraction(1)
    for i in range(k):
        r *= Fraction(comb(N + k + i, k), comb(k + i, k))
    return r


def ks_factor(j, k):
    """Gamma(j)Gamma(j+2k)/Gamma(j+k)^2 for integer k, exactly."""
    num = den = 1
    for i in range(k):
        num *= j + k + i
        den *= j + i
    return Fraction(num, den)


def bhny_factor(j, k):
    """E|1 + e^{iw} sqrt(B)|^{2k}, w uniform, B ~ Beta(1, j-1) (B = 1 for j = 1), exactly."""
    s = Fraction(0)
    for m in range(k + 1):
        rising = 1
        for i in range(m):
            rising *= j + i
        s += comb(k, m) ** 2 * Fraction(factorial(m), rising)
    return s


def toeplitz_dets(k, nmax):
    """D_n = det(c_{r-s})_{r,s<n} for n = 0..nmax, c_m = (-1)^m binom(2k,k+m): Gaussian elimination without pivoting
    (the matrix is symmetric positive definite, so every pivot D_{i+1}/D_i is positive and the band never fills)."""
    c = {m: Fraction((-1) ** abs(m) * comb(2 * k, k + m)) for m in range(-k, k + 1)}

    def row(r):
        return {s: c[r - s] for s in range(max(0, r - k), min(nmax, r + k + 1))}

    rows = {r: row(r) for r in range(min(nmax, k + 1))}
    D, dets, minpiv = Fraction(1), [Fraction(1)], None
    for i in range(nmax):
        piv = rows[i][i]
        assert piv > 0
        minpiv = piv if minpiv is None else min(minpiv, piv)
        D *= piv
        dets.append(D)
        for r in range(i + 1, min(nmax, i + k + 1)):
            if r not in rows:
                rows[r] = row(r)
            f = rows[r].pop(i, 0) / piv
            if f:
                for s, v in rows[i].items():
                    if s > i:
                        rows[r][s] = rows[r].get(s, 0) - f * v
        del rows[i]
    return dets, minpiv


out = dict(script='rmt_ratios.py', host=socket.gethostname(), slurm_job=os.environ.get('SLURM_JOB_ID'),
           versions=dict(python_flint=flint.__version__, mpmath=mp.__version__),
           definition='ratio = E|Z_N|^{2k} / N^{k^2}, E|Z_N|^{2k} = prod_{j=1}^N Gamma(j)Gamma(j+2k)/Gamma(j+k)^2',
           barnes={}, ratios=[], checks={})

# Barnes constants G(1+k)^2/G(1+2k)
def G_int(n):
    r = 1
    for m in range(n - 1):
        r *= factorial(m)
    return r


BAR = {}
for k in KS:
    ex = Fraction(G_int(1 + k) ** 2, G_int(1 + 2 * k))
    BAR[k] = ex
    with mp.workdps(30):
        bm = mp.barnesg(1 + k) ** 2 / mp.barnesg(1 + 2 * k)
        bm_rel = abs(bm / (mp.mpf(ex.numerator) / ex.denominator) - 1)
    ba = acb(1 + k).barnes_g() ** 2 / acb(1 + 2 * k).barnes_g()
    row = dict(k=k, exact=str(ex), value=float(ex), mpmath=float(bm), arb=str(ba.real),
               arb_contains_exact=bool(ba.real.overlaps(arb(ex.numerator) / arb(ex.denominator))),
               mpmath_rel_err=float(bm_rel))
    if k in QUOTED_BARNES:
        q = QUOTED_BARNES[k]
        row['quoted'] = q[0] if q[1] is None else '%se%d' % q
        row['quoted_ok'] = (q[0] == str(ex)) if q[1] is None else (rounded_like(ex, q) == q[0])
    out['barnes'][str(k)] = row
    print('k', k, 'G(1+k)^2/G(1+2k) =', ex, '=', float(ex), ' mpmath', float(bm), ' Arb', ba.real, flush=True)

allok = True
for k in KS:
    t0 = time.time()
    # (2) exact running product over j and (3) beta decomposition, compared with the closed form at every N <= NMAX
    E, EB, bad_prod, bad_bhny = Fraction(1), Fraction(1), 0, 0
    exact_at = {}
    for j in range(1, NMAX + 1):
        E *= ks_factor(j, k)
        EB *= bhny_factor(j, k)
        cf = closed(j, k)
        bad_prod += (E != cf)
        bad_bhny += (EB != cf)
        if j in NS:
            exact_at[j] = E
    # (4) Heine's identity, every N <= NMAX
    dets, minpiv = toeplitz_dets(k, NMAX)
    bad_toep = sum(dets[n] != closed(n, k) for n in range(1, NMAX + 1))
    # (5) Arb enclosure by lgamma sums
    s, arb_at = arb(0), {}
    for j in range(1, NMAX + 1):
        s += arb(j).lgamma() + arb(j + 2 * k).lgamma() - 2 * arb(j + k).lgamma()
        if j in NS:
            arb_at[j] = s.exp() / arb(j) ** (k * k)
    out['checks'][str(k)] = dict(N_checked=NMAX, product_vs_closed_form_mismatches=int(bad_prod),
                                 beta_decomposition_vs_closed_form_mismatches=int(bad_bhny),
                                 toeplitz_vs_closed_form_mismatches=int(bad_toep),
                                 toeplitz_smallest_pivot=float(minpiv), seconds=round(time.time() - t0, 2))
    allok &= bad_prod == 0 and bad_bhny == 0 and bad_toep == 0
    print('k', k, 'mismatches with the closed form for N <= %d: product %d, beta decomposition %d, Toeplitz %d'
          % (NMAX, bad_prod, bad_bhny, bad_toep), flush=True)
    for N in NS:
        Ex = exact_at[N]
        assert Ex == closed(N, k) == dets[N]
        x = Ex / Fraction(N) ** (k * k)
        # (1) the formula of killing2.py (mp.dps = 30 there)
        with mp.workdps(30):
            ex_mp = mp.fprod([mp.gamma(j) * mp.gamma(j + 2 * k) / mp.gamma(j + k) ** 2 for j in range(1, N + 1)])
            ratio_mp = ex_mp / N ** (k * k)
        with mp.workdps(50):
            xm = mp.mpf(x.numerator) / x.denominator
            rel_mp = abs(ratio_mp - xm) / xm
            excess = (xm / (mp.mpf(BAR[k].numerator) / BAR[k].denominator) - 1) * N / k ** 3
        ball = arb_at[N]
        in_ball = bool(ball.overlaps(arb(x.numerator) / arb(x.denominator)))
        allok &= in_ball and rel_mp < 1e-20
        row = dict(N=N, k=k, exact_moment=str(Ex), ratio=float(x), ratio_30_digits=mp.nstr(xm, 30),
                   killing2_formula_ratio=float(ratio_mp), killing2_formula_rel_err=float(rel_mp),
                   arb_ratio=str(ball), arb_contains_exact=in_ball,
                   barnes=float(BAR[k]), ratio_over_barnes=float(x / BAR[k]),
                   excess_times_N_over_k3=float(excess))
        if k in QUOTED and N in QUOTED[k]:
            q = QUOTED[k][N]
            row['quoted'] = q[0] if q[1] == 0 else '%se%d' % q
            row['exact_rounded_like_quoted'] = rounded_like(x, q)
            row['quoted_ok'] = row['exact_rounded_like_quoted'] == q[0]
        out['ratios'].append(row)
        print('N %6d k %d  E|Z_N|^{2k} = %s  ratio %s  (killing2 formula %.12g, Arb %s)  quoted %s ok %s  (ratio/Barnes-1)N/k^3 %.6f'
              % (N, k, Ex if Ex < 10 ** 30 else '%.6e' % float(Ex), mp.nstr(xm, 12), float(ratio_mp), ball.str(12, radius=False),
                 row.get('quoted'), row.get('quoted_ok'), float(excess)), flush=True)

# comparison with the archived output of killing2.py (N = 10, 20; k = 1, 2)
if os.path.exists('killing.json'):
    arch = json.load(open('killing.json')).get('ks', [])
    cmp = []
    for r in arch:
        N, k = r['N'], r['k']
        ex = closed(N, k)
        cmp.append(dict(N=N, k=k, archived_exact=r['exact'], exact=str(ex), archived_ratio=r['ratio'],
                        ratio=float(ex / Fraction(N) ** (k * k)),
                        agree=abs(r['exact'] - float(ex)) <= 1e-9 * float(ex)
                        and abs(r['ratio'] - float(ex / Fraction(N) ** (k * k))) <= 1e-12))
    out['killing_json_comparison'] = cmp
    allok &= all(c['agree'] for c in cmp)
    print('killing.json entries agree:', all(c['agree'] for c in cmp), flush=True)

quoted_rows = [r for r in out['ratios'] if 'quoted' in r] + [b for b in out['barnes'].values() if 'quoted' in b]
out['all_quoted_values_correct'] = all(r['quoted_ok'] for r in quoted_rows)
out['all_checks_pass'] = bool(allok)
out['seconds'] = round(time.time() - t_start, 1)
print('all quoted values correct:', out['all_quoted_values_correct'], ' all checks pass:', out['all_checks_pass'],
      ' seconds', out['seconds'])
json.dump(out, open('rmt_ratios.json', 'w'), indent=1)
