# CHECK SCRIPT for the reconstruction in misc/repro/reconstruct/hankel. It is neither an archive script nor
# the authors'; it tests the reconstructed hankel.py and the archived data/hankel.json in four ways.
# 1. It encloses all 18 numbers in Arb ball arithmetic (python-flint), from the first 1700 zeros computed by
#    Arb itself (acb.zeta_zeros) and t = 1/100, 1/20 exactly. The balls are rigorous for the quantities as the
#    book defines them from 1700 zeros; the precision is doubled until every relative radius is below 1e-30.
# 2. It repeats the enclosure with the float64 ordinates of g_*.npy taken as exact, and with the first 80, 100
#    and 120 zeros only (the round 1 recomputation used 80, Exercise 1 of Chapter ch:heat asks for 100, and
#    Numerical observation 7.2 of the September 2026 draft of P17 used 120).
# 3. It reruns the mpmath method of hankel.py at 20 to 100 digits, to show the precision the stored strings need,
#    and once with t as an exact decimal instead of a float.
# 4. It compares the regenerated hankel.json with the archived file number by number (the walk of
#    misc/repro/compare.py), tests whether each archived string is the correctly rounded value of its enclosure,
#    and tests Table tab:ch12:hankel (ch/ch12.tex lines 55 to 58), the orders of magnitude quoted in the text,
#    the round 1 recomputation (misc/reviews/round1/U09-review.md line 275) and Observation 7.2 the same way.
# Run after hankel.py: python3 hankel_check.py ARCHIVED_JSON   (writes hankel_check.json)
import sys, json, time
import numpy as np, mpmath as mp
from flint import arb, acb, arb_mat, ctx

ARCH = sys.argv[1] if len(sys.argv) > 1 else '/scratch/vsokolov/rh_book_repro/book_data/hankel.json'
TS = ['0.01', '0.05']
TEXACT = {'0.01': (1, 100), '0.05': (1, 20)}
mp.mp.dps = 60
T0 = time.time()
report = {}

def say(*a):
    print(*a, flush=True)

# ---------- Arb enclosures ----------
def quantities_arb(gam, t, nz):
    """c_0 c_2 / c_1^2 and the normalised even and odd minors, N = 2..5, from the first nz ordinates."""
    a = [x*x for x in gam[:nz]]
    p = [(-x*t).exp() for x in a]
    c = []
    for k in range(10):
        s = arb(0)
        for y in p:
            s += y
        c.append(s)
        p = [y*x for x, y in zip(a, p)]
    ev = []; od = []
    for N in range(2, 6):
        de = arb_mat([[c[i+j] for j in range(N)] for i in range(N)]).det()
        do = arb_mat([[c[i+j+1] for j in range(N)] for i in range(N)]).det()
        pe = arb(1); po = arb(1)
        for i in range(N):
            pe *= c[2*i]; po *= c[2*i+1]
        ev.append(de/pe); od.append(do/po)
    return dict(ratio=c[0]*c[2]/(c[1]*c[1]), even=ev, odd=od)

def flat(q):
    out = []
    for t in TS:
        out.append((f'/{t}/ratio', q[t]['ratio'], 10))
        for key in ['even', 'odd']:
            for i, v in enumerate(q[t][key]):
                out.append((f'/{t}/{key}[{i}]', v, 4))
    return out

def mpv(x, n=50):
    return mp.mpf(x.str(n, radius=False))

def lohi(x):
    return mpv(x.lower()), mpv(x.upper())

def mid(x):
    return mpv(x.mid())

def relrad(x):
    return mp.mpf(x.rad().str(5, radius=False))/abs(mid(x))

def rounds_to(x, s, d):
    """True if every point of the ball (arb) or the number (mpf) x rounds to the decimal string s at d digits."""
    lo, hi = lohi(x) if isinstance(x, arb) else (x, x)
    return mp.mpf(mp.nstr(lo, d)) == mp.mpf(s) and mp.mpf(mp.nstr(hi, d)) == mp.mpf(s)

def sig(s):
    m = s.strip().lstrip('+-').lower().split('e')[0].replace('.', '').lstrip('0')
    return len(m)

g = np.concatenate([np.load(f) for f in ['g_1_400.npy', 'g_401_1000.npy', 'g_1001_1700.npy']])
say('float64 zeros', len(g), g[0], g[-1])

prec = 320
while True:
    ctx.prec = prec
    t1 = time.time()
    zs = acb.zeta_zeros(1, 1700)
    gam = [z.imag for z in zs]
    tz = time.time() - t1
    tt = {t: arb(TEXACT[t][0])/TEXACT[t][1] for t in TS}
    A = {t: quantities_arb(gam, tt[t], 1700) for t in TS}
    worst = max(relrad(v) for _, v, _ in flat(A))
    say(f'Arb at {prec} bits: zeros in {tz:.1f} s, worst relative radius {mp.nstr(worst, 3)}')
    if worst < mp.mpf('1e-30'):
        break
    prec *= 2
report['arb_prec_bits'] = prec
g64 = [arb(float(x)) for x in g]
maxzd = max(abs(mid(x - y)) for x, y in zip(gam, g64))
say('largest |Arb zero minus float64 zero| over 1700 zeros:', mp.nstr(maxzd, 3))
report['max_abs_zero_diff_float64'] = mp.nstr(maxzd, 3)
A64 = {t: quantities_arb(g64, tt[t], 1700) for t in TS}
AN = {nz: {t: quantities_arb(gam, tt[t], nz) for t in TS} for nz in [80, 100, 120]}

fa = flat(A); f64 = flat(A64)
say('\nArb enclosures from 1700 Arb zeros (midpoint to 20 digits, relative radius)')
report['arb'] = {}
for (p, v, _), (_, w, _) in zip(fa, f64):
    r64 = abs(mid(w)/mid(v) - 1)
    report['arb'][p] = dict(mid=mp.nstr(mid(v), 20), relrad=mp.nstr(relrad(v), 3), rel_effect_float64_zeros=mp.nstr(r64, 3))
    say(f'  {p:14s} {mp.nstr(mid(v), 20):>28s}  relrad {mp.nstr(relrad(v), 3):>9s}  float64 zeros change it by {mp.nstr(r64, 3)}')
report['truncation'] = {}
for nz in [80, 100, 120]:
    worst = mp.mpf(0)
    for (p, v, _), (_, w, _) in zip(fa, flat(AN[nz])):
        d = abs((w - v)/v).upper()
        worst = max(worst, mpv(d))
    report['truncation'][str(nz)] = mp.nstr(worst, 3)
    say(f'first {nz} zeros instead of 1700: largest relative change of the 18 numbers at most {mp.nstr(worst, 3)}')

# ---------- mpmath, as in hankel.py ----------
def quantities_mp(dps, tdec=False, nz=1700):
    with mp.workdps(dps):
        a = [mp.mpf(float(x))**2 for x in g[:nz]]
        q = {}
        for t in TS:
            T = mp.mpf(t) if tdec else mp.mpf(float(t))
            e = [mp.exp(-x*T) for x in a]
            c = [mp.fsum(x**k*y for x, y in zip(a, e)) for k in range(10)]
            ev = []; od = []
            for N in range(2, 6):
                He = mp.matrix([[c[i+j] for j in range(N)] for i in range(N)])
                Ho = mp.matrix([[c[i+j+1] for j in range(N)] for i in range(N)])
                ev.append(mp.det(He)/mp.fprod(c[2*i] for i in range(N)))
                od.append(mp.det(Ho)/mp.fprod(c[2*i+1] for i in range(N)))
            q[t] = dict(ratio=c[0]*c[2]/c[1]**2, even=ev, odd=od)
    return q

arch = json.load(open(ARCH))
archflat = {f'/{t}/ratio': arch[t]['ratio'] for t in TS}
for t in TS:
    for key in ['even', 'odd']:
        for i, s in enumerate(arch[t][key]):
            archflat[f'/{t}/{key}[{i}]'] = s
say('\nmpmath method of hankel.py (float64 zeros, t as float) at several precisions')
report['precision_scan'] = {}
for dps in [20, 25, 30, 35, 40, 45, 50, 55, 60, 80, 100]:
    fq = flat(quantities_mp(dps))
    same = sum(mp.nstr(v, d) == archflat[p] for p, v, d in fq)
    err = max(abs(v/mid(w) - 1) for (p, v, d), (_, w, _) in zip(fq, f64))
    report['precision_scan'][str(dps)] = dict(strings_equal_to_archive=same, max_rel_err_vs_arb_float64_zeros=mp.nstr(err, 3))
    say(f'  dps {dps:3d}: {same:2d} of 18 strings equal the archive, largest relative error {mp.nstr(err, 3)}')
q60 = flat(quantities_mp(60)); q60d = flat(quantities_mp(60, tdec=True))
dt = max(abs(v/w - 1) for (_, v, _), (_, w, _) in zip(q60, q60d))
report['t_float_vs_decimal_max_rel'] = mp.nstr(dt, 3)
say('t as float instead of exact decimal changes the 18 numbers by at most', mp.nstr(dt, 3))
e60 = max(abs(v/mid(w) - 1) for (_, v, _), (_, w, _) in zip(q60, fa))
report['mp60_vs_arb_exact_zeros_max_rel'] = mp.nstr(e60, 3)
say('mpmath at 60 digits from float64 zeros against Arb from exact zeros: largest relative difference', mp.nstr(e60, 3))

# ---------- regenerated file against the archived file ----------
def num(x):
    if isinstance(x, bool): return None
    if isinstance(x, (int, float)): return mp.mpf(x)
    if isinstance(x, str):
        try: return mp.mpf(x.strip().replace('[', '').split('+/-')[0])
        except Exception: return None
    return None
def walk(a, b, path, diffs, stats):
    if isinstance(a, dict) and isinstance(b, dict):
        for k in a:
            if k not in b: diffs.append((path + '/' + str(k), 'missing in rerun', '')); continue
            walk(a[k], b[k], path + '/' + str(k), diffs, stats)
        for k in b:
            if k not in a: diffs.append((path + '/' + str(k), 'extra in rerun', ''))
        return
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b): diffs.append((path, f'length {len(a)} vs {len(b)}', ''))
        for i, (x, y) in enumerate(zip(a, b)): walk(x, y, f'{path}[{i}]', diffs, stats)
        return
    na, nb = num(a), num(b)
    if na is not None and nb is not None:
        stats['n'] += 1
        if na == nb: return
        d = abs(na - nb) / max(abs(na), abs(nb), mp.mpf(10) ** -300)
        stats['maxrel'] = max(stats['maxrel'], d)
        if d > mp.mpf('1e-9'): diffs.append((path, mp.nstr(na, 15), mp.nstr(nb, 15) + f'  (rel {mp.nstr(d, 3)})'))
        return
    if a != b: diffs.append((path, str(a)[:60], str(b)[:60]))
regen = json.load(open('hankel.json'))
diffs = []; stats = {'n': 0, 'maxrel': mp.mpf(0)}
walk(arch, regen, '', diffs, stats)
same_bytes = open('hankel.json', 'rb').read() == open(ARCH, 'rb').read()
report['regenerated_vs_archived'] = dict(numbers=stats['n'], max_rel=mp.nstr(stats['maxrel'], 3), differences=diffs, byte_identical=same_bytes)
say(f'\nregenerated hankel.json against {ARCH}: {stats["n"]} numbers, {len(diffs)} differences beyond 1e-9, max rel diff {mp.nstr(stats["maxrel"], 3)}, byte identical {same_bytes}')

# ---------- archived strings, book table, round 1 values and Observation 7.2 against the enclosures ----------
def test(name, values, ref):
    """values: path -> decimal string; ref: path -> arb. Correct rounding at the printed digits, and deviation."""
    rows = {}; ok = 0; worst = mp.mpf(0)
    for p, s in values.items():
        v = ref[p]; d = sig(s)
        good = rounds_to(v, s, d)
        rel = abs(mp.mpf(s)/mid(v) - 1)
        worst = max(worst, rel); ok += good
        rows[p] = dict(printed=s, digits=d, correctly_rounded=good, rel_dev=mp.nstr(rel, 3), true=mp.nstr(mid(v), 12))
    report[name] = dict(n=len(values), correctly_rounded=ok, max_rel_dev=mp.nstr(worst, 3), rows=rows)
    say(f'{name}: {ok} of {len(values)} printed values are the correctly rounded true values; largest relative deviation {mp.nstr(worst, 3)}')
    for p, r in rows.items():
        if not r['correctly_rounded']:
            say(f'    NOT correctly rounded: {p} printed {r["printed"]} true {r["true"]}')
ref = {p: v for p, v, _ in fa}
ref120 = {p: v for p, v, _ in flat(AN[120])}
test('archived_strings', archflat, ref)
BOOK = {'/0.01/ratio': '1.1346548', '/0.01/even[0]': '0.119', '/0.01/even[1]': '2.45e-3', '/0.01/even[2]': '1.05e-5', '/0.01/even[3]': '1.04e-8',
        '/0.01/odd[0]': '0.182', '/0.01/odd[1]': '5.87e-3', '/0.01/odd[2]': '3.30e-5', '/0.01/odd[3]': '4.57e-8',
        '/0.05/ratio': '1.0000081', '/0.05/even[0]': '8.11e-6', '/0.05/even[1]': '1.77e-14', '/0.05/even[2]': '5.35e-28', '/0.05/even[3]': '2.20e-44',
        '/0.05/odd[0]': '1.79e-5', '/0.05/odd[1]': '1.22e-13', '/0.05/odd[2]': '1.71e-26', '/0.05/odd[3]': '3.81e-42'}
test('book_table_tab_ch12_hankel', BOOK, ref)
REVIEW = {'/0.01/ratio': '1.13465478943', '/0.05/ratio': '1.00000811308',
          '/0.01/even[0]': '0.118675', '/0.01/even[1]': '2.45125e-3', '/0.01/even[2]': '1.0458e-5', '/0.01/even[3]': '1.04143e-8',
          '/0.01/odd[0]': '0.181646', '/0.01/odd[1]': '5.86623e-3', '/0.01/odd[2]': '3.29553e-5', '/0.01/odd[3]': '4.56716e-8',
          '/0.05/even[0]': '8.11302e-6', '/0.05/even[1]': '1.76901e-14', '/0.05/even[2]': '5.34791e-28', '/0.05/even[3]': '2.19922e-44',
          '/0.05/odd[0]': '1.79472e-5', '/0.05/odd[1]': '1.22491e-13', '/0.05/odd[2]': '1.71497e-26', '/0.05/odd[3]': '3.81297e-42'}
test('round1_U09_recomputation', REVIEW, ref)
OBS = {'/0.01/ratio': '1.1347', '/0.05/ratio': '1.0000081',
       '/0.01/even[0]': '0.119', '/0.01/even[1]': '2.45e-3', '/0.01/even[2]': '1.05e-5', '/0.01/even[3]': '1.04e-8',
       '/0.01/odd[0]': '0.182', '/0.01/odd[1]': '5.87e-3', '/0.01/odd[2]': '3.30e-5', '/0.01/odd[3]': '4.57e-8'}
test('P17_draft_observation_7_2_vs_120_zeros', OBS, ref120)
# Could the table have been rounded from the stored 4 digit strings? Round each stored string to 3 digits
# through a binary double, as a %.2e format would.
viaf = sum(mp.mpf('%.2e' % float(archflat[p])) == mp.mpf(s) for p, s in BOOK.items() if 'ratio' not in p)
report['table_minors_equal_to_%.2e_of_stored_strings'] = viaf
say(f'table minors equal to %.2e of the stored 4 digit strings: {viaf} of 16')
mags = {p: int(mp.nint(mp.log10(mid(ref[p])))) for p in ref if p.startswith('/0.05/')}
report['orders_of_magnitude_t_0.05'] = mags
say('nearest powers of ten at t = 0.05:', mags)
report['seconds'] = round(time.time() - T0, 1)
json.dump(report, open('hankel_check.json', 'w'), indent=1)
say('done in', report['seconds'], 's')
