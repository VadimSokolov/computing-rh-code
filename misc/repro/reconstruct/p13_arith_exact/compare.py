# New comparison script (item arith_exact): the patched arith.json (60 digits) against the archived data/arith.json
# (here arith_archived.json), entry by entry, and the patched runs at 40 and 80 digits against the one at 60 digits;
# then the pixels of the figure drawn by heattrace.py from the new and from the archived json.
# Usage: python3 compare.py   (writes compare.json)
import json, mpmath as mp, numpy as np
from PIL import Image
mp.mp.dps = 50
new = json.load(open('arith.json')); old = json.load(open('arith_archived.json'))
d40 = json.load(open('arith_dps40.json')); d80 = json.load(open('arith_dps80.json'))
names = ['t', 'Pi', 'A', 'P', 'Pi+A-P', 'W zeros', 'digits lost']
res = []
print('new (60 digits) against archived data/arith.json')
for rn, ro in zip(new, old):
    assert rn[0] == ro[0]
    for j in range(1, 7):
        a, b = mp.mpf(rn[j]), mp.mpf(ro[j])
        same = str(rn[j]) == str(ro[j])
        rel = abs(a - b)/abs(a)
        res.append(dict(t=rn[0], entry=names[j], archived=ro[j], new=rn[j], identical=same, abs_diff=mp.nstr(abs(a - b), 4), rel_diff=mp.nstr(rel, 4)))
        if not same:
            print('  t=%-6s %-12s archived %-24s new %-24s |diff| %-10s rel %s' % (rn[0], names[j], ro[j], rn[j], mp.nstr(abs(a - b), 4), mp.nstr(rel, 4)))
nid = sum(r['identical'] for r in res)
print('  %d of %d entries identical as stored' % (nid, len(res)))
print('runs at 40 and 80 digits against 60 digits (entries that differ as stored)')
ctrl = []
for lab, D in [('40', d40), ('80', d80)]:
    k = 0
    for rn, rd in zip(new, D):
        for j in range(1, 7):
            if str(rn[j]) != str(rd[j]):
                k += 1
                a, b = mp.mpf(rn[j]), mp.mpf(rd[j])
                ctrl.append(dict(dps=lab, t=rn[0], entry=names[j], dps60=rn[j], other=rd[j], rel_diff=mp.nstr(abs(a - b)/abs(a), 4)))
                print('  dps %s t=%-6s %-12s 60: %-24s %s: %-24s rel %s' % (lab, rn[0], names[j], rn[j], lab, rd[j], mp.nstr(abs(a - b)/abs(a), 4)))
    print('  dps %s: %d entries differ from dps 60' % (lab, k))
pix = {}
for a, b in [('heattrace.png', 'heattrace_archived_json.png')]:
    x = np.asarray(Image.open(a).convert('RGB'), dtype=int); y = np.asarray(Image.open(b).convert('RGB'), dtype=int)
    pix[a + ' vs ' + b] = dict(shape=list(x.shape), same_shape=x.shape == y.shape, differing_pixels=int((np.abs(x - y).sum(axis=2) > 0).sum()) if x.shape == y.shape else None)
print('pixels:', pix)
json.dump(dict(new_vs_archived=res, controls=ctrl, pixels=pix), open('compare.json', 'w'), indent=1)
