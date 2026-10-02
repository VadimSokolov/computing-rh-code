# Check helper of unit R2-23 (not an archive script): compare two JSON files number by number.
# Usage: python3 cmpjson.py ARCHIVED.json NEW.json [TOL]   prints the count of numbers, how many are bit identical,
# the largest relative difference and the largest difference in units in the last place (floats only).
import json, sys, math, hashlib


def flat(a, path=''):
    if isinstance(a, dict):
        for k in a:
            yield from flat(a[k], path + '/' + str(k))
    elif isinstance(a, list):
        for i, v in enumerate(a):
            yield from flat(v, '%s[%d]' % (path, i))
    else:
        yield path, a


fa, fb = sys.argv[1], sys.argv[2]
tol = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-9
ra, rb = open(fa, 'rb').read(), open(fb, 'rb').read()
print('sha256 %s  %s' % (hashlib.sha256(ra).hexdigest(), fa))
print('sha256 %s  %s' % (hashlib.sha256(rb).hexdigest(), fb))
print('byte identical:', ra == rb)
A, B = dict(flat(json.loads(ra))), dict(flat(json.loads(rb)))
print('same paths in the same order:', list(A) == list(B))
n = same = over = 0
worst = (0.0, None, None, None)
maxulp = 0
for p in A:
    x, y = A[p], B.get(p)
    if isinstance(x, bool) or not isinstance(x, (int, float)):
        same += (x == y); n += 1
        continue
    n += 1
    if x == y and type(x) == type(y):
        same += 1
        continue
    d = abs(x - y) / max(abs(x), abs(y))
    if d > tol:
        over += 1
    if isinstance(x, float) and isinstance(y, float):
        u = abs(x - y) / math.ulp(max(abs(x), abs(y)))
        maxulp = max(maxulp, u)
    if d > worst[0]:
        worst = (d, p, x, y)
print('%d values, %d identical, %d beyond relative %g, max relative difference %.3g%s, max ulp %.0f'
      % (n, same, over, tol, worst[0], '' if worst[1] is None else ' at %s (archived %r, new %r)' % worst[1:], maxulp))
