# Compare the output of the reconstructed spacing.py with the archived data/spacing.json, number by number, with
# the method of misc/repro/compare.py (60 digit arithmetic, relative tolerance 1e-9).
# Run: python3 compare_spacing.py ARCHIVED.json NEW.json [NEW2.json ...]   (writes compare_spacing.txt)
import json, sys
import mpmath as mp
mp.mp.dps = 60
TOL = mp.mpf('1e-9')


def flat(a, path=''):
    if isinstance(a, dict):
        for k in a: yield from flat(a[k], path + '/' + str(k))
    elif isinstance(a, list):
        for i, v in enumerate(a): yield from flat(v, '%s[%d]' % (path, i))
    else:
        yield path, a


def compare(fa, fb):
    A, B = json.load(open(fa)), json.load(open(fb))
    a, b = dict(flat(A)), dict(flat(B))
    lines = ['%s against %s' % (fb, fa)]
    if list(A) != list(B): lines.append('  keys differ or are in another order: %s | %s' % (list(A), list(B)))
    miss = [p for p in a if p not in b] + [p for p in b if p not in a]
    if miss: lines.append('  paths present in one file only: %s' % miss)
    rows = []
    for p in a:
        if p not in b: continue
        x, y = mp.mpf(a[p]), mp.mpf(b[p])
        d = mp.mpf(0) if x == y else abs(x - y) / max(abs(x), abs(y))
        rows.append((d, p, a[p], b[p]))
    n, same = len(rows), sum(1 for r in rows if r[0] == 0)
    bad = [r for r in rows if r[0] > TOL]
    worst = sorted(rows, key=lambda r: -r[0])[:6]
    lines.append('  %d numbers, %d identical in float64, %d beyond relative 1e-9, max relative difference %s'
                 % (n, same, len(bad), mp.nstr(worst[0][0], 3)))
    for d, p, x, y in worst:
        if d > 0: lines.append('    %-12s archived %-24r new %-24r rel %s' % (p, x, y, mp.nstr(d, 3)))
    return lines


if __name__ == '__main__':
    out = []
    for f in sys.argv[2:]: out += compare(sys.argv[1], f)
    for i in range(2, len(sys.argv)):
        for j in range(i + 1, len(sys.argv)): out += compare(sys.argv[i], sys.argv[j])
    open('compare_spacing.txt', 'w').write('\n'.join(out) + '\n'); print('\n'.join(out))
