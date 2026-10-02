# Reconstruction helper (not an authors' script): compares the rerun riesz_explicit.json (argument 1) with the archived copies
# basepoint_one/wiener/riesz_explicit.json and riesz_program/data/riesz_explicit.json (arguments 2 and 3), number by number,
# with the relative tolerance 1e-9 of misc/repro/compare.py, and also byte by byte.
import json, sys, mpmath as mp
mp.mp.dps = 60
new, arch = sys.argv[1], sys.argv[2:]
def flat(a, path=''):
    if isinstance(a, dict):
        for k in a: yield from flat(a[k], path + '/' + k)
    elif isinstance(a, list):
        for i, v in enumerate(a): yield from flat(v, '%s[%d]' % (path, i))
    else: yield path, a
b = json.load(open(new)); bb = open(new, 'rb').read()
for f in arch:
    a = json.load(open(f)); ab = open(f, 'rb').read()
    print('== %s against %s: byte identical %s (%d and %d bytes)' % (new, f, ab == bb, len(bb), len(ab)))
    fa, fb = dict(flat(a)), dict(flat(b))
    print('   keys equal:', list(fa) == list(fb))
    n = agree = exact = 0; worst = (mp.mpf(0), '')
    for k in fa:
        x, y = mp.mpf(fa[k]), mp.mpf(fb[k]); n += 1
        d = abs(x - y) / max(abs(x), abs(y), mp.mpf(10) ** -300)
        exact += (fa[k] == fb[k]); agree += (d <= mp.mpf('1e-9'))
        if d >= worst[0]: worst = (d, k)
        print('   %-10s archived %-24r rerun %-24r rel diff %s' % (k, fa[k], fb[k], mp.nstr(d, 3)))
    print('   %d numbers, %d agree to 1e-9, %d bit identical, worst rel diff %s at %s' % (n, agree, exact, mp.nstr(worst[0], 3), worst[1]))
if len(arch) == 2:
    print('== the two archived copies are byte identical:', open(arch[0], 'rb').read() == open(arch[1], 'rb').read())
