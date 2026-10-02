# saddle.py: Table tab:ch12:saddle of Section ch:heat (ch/ch12.tex), added for the book in September 2026 from the check
# script of the GPT review (misc/gpt-review/triage/numerics/N1/b_saddle_table.py); no earlier script wrote the table.
# For y = 50, 200, 1000 it encloses in Arb ball arithmetic (python-flint, at 128 and at 256 bits)
#   col1 = log(Xi(iy)/Xi(0)) = log xi(1/2 + y) - log xi(1/2)   (Xi(iy) = xi(1/2 - y) = xi(1/2 + y) by the functional equation),
#   col2 = (y/2) log(y/(2 pi e)),   and the ratio col1/col2,
# with xi(s) = s (s - 1) pi^(-s/2) Gamma(s/2) zeta(s) / 2, and the difference col1 - col2 - (7/4) log y, which Stirling's
# formula sends to C = (1/4) log(pi/2) - log xi(1/2) = 0.81182...  It checks that every printed entry is its enclosure
# rounded to three decimals and that mpmath at 50 digits falls in the enclosures of col1 and col2, prints the table and
# writes saddle.json.
# Run on Hopper (one core, under a second), from the folder data/, which receives saddle.json:
#   cd data && bash ../misc/tools/hopper_run.sh -c 1 -m 2G -t 10 -g saddle.json ../code/saddle.py
import json, time
from flint import arb, ctx
import mpmath as mp

BOOK = {50: ('34.511', '26.854', '1.285'), 200: ('256.128', '246.044', '1.041'), 1000: ('2047.839', '2034.939', '1.006')}
t0 = time.time()

def logxi(s):
    # log xi(s) for real s > 0; at s = 1/2 both s (s - 1) and zeta(1/2) are negative
    return (abs(s * (s - 1)) / 2).log() - s / 2 * arb.pi().log() + (s / 2).lgamma() + abs(s.zeta()).log()

out = {'book': {str(y): v for y, v in BOOK.items()}}
checks = {}
for prec in (128, 256):
    ctx.prec = prec
    half = arb(1) / 2
    L0 = logxi(half)
    rows = {'xi_half': L0.exp().str(20), 'C': ((arb.pi() / 2).log() / 4 - L0).str(20)}
    ball = {}
    for y in BOOK:
        Y = arb(y)
        c1 = logxi(half + Y) - L0
        c2 = Y / 2 * (Y / (2 * arb.pi() * arb.const_e())).log()
        ball[y] = (c1, c2, c1 / c2)
        rows[str(y)] = dict(col1=c1.str(25), col2=c2.str(25), ratio=(c1 / c2).str(20), diff=(c1 - c2).str(15),
                            diff_minus_7_4_logy=(c1 - c2 - arb(7) / 4 * Y.log()).str(15))
        for name, x, b in zip(('col1', 'col2', 'ratio'), ball[y], BOOK[y]):
            checks['%d_%s_%d' % (y, name, prec)] = bool(abs(x - arb(b)) < arb(5) / 10000)
    out['arb_%d' % prec] = rows

mp.mp.dps = 50
def mlogxi(s):
    return mp.log(abs(s * (s - 1) * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s) / 2))
M0 = mlogxi(mp.mpf(1) / 2)
out['mpmath_50'] = {}
for y in BOOK:
    m1 = mlogxi(mp.mpf(1) / 2 + y) - M0
    m2 = mp.mpf(y) / 2 * mp.log(y / (2 * mp.pi * mp.e))
    out['mpmath_50'][str(y)] = dict(col1=mp.nstr(m1, 25), col2=mp.nstr(m2, 25))
    for name, m, x in (('col1', m1, ball[y][0]), ('col2', m2, ball[y][1])):
        checks['%d_%s_mpmath' % (y, name)] = bool(abs(arb(mp.nstr(m, 45)) - x) < arb(10) ** -40)

out['checks'] = checks
out['all_checks_true'] = all(checks.values())
out['elapsed_s'] = round(time.time() - t0, 2)
json.dump(out, open('saddle.json', 'w'), indent=1)
for y in BOOK:
    r = out['arb_256'][str(y)]
    print(y, r['col1'][:14], r['col2'][:14], r['ratio'][:10], '  book', *BOOK[y])
print('C =', out['arb_256']['C'], '  checks true:', sum(checks.values()), 'of', len(checks))
if not out['all_checks_true']:
    raise SystemExit('a check failed: ' + ', '.join(k for k, v in checks.items() if not v))
