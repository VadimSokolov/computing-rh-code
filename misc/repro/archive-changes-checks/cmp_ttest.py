# Check helper of unit R2-23 (not an archive script): compare a new run of ttest2.py and ttest.py with the archived
# reconstruction run (misc/repro/reconstruct/ttest2/out/archived_1241681.log), section by section, ignoring timings.
# Usage: python3 cmp_ttest.py ARCHIVED.log NEW.log
import sys


def sections(fname):
    out, cur = {}, None
    for line in open(fname):
        line = line.rstrip('\n')
        if line.startswith('== '):
            cur = line[3:].split(' (')[0]
            out[cur] = []
        elif cur is not None and line.strip():
            out[cur].append(line)
    return out


def key(line):
    t = line.split()
    if t[0] == 'nodes':
        return t
    if t[0] == 'setup':
        return [t[0]] + t[2:]            # drop the setup time
    if t[1] == 'err':
        return t[:5]                     # height, err, relrad; drop 'sec' and the time
    return t[:-1]                        # ttest.py: height, error, ball radius; drop the time


A, B = sections(sys.argv[1]), sections(sys.argv[2])
nall = nok = 0
for s in A:
    if s == 'done':
        continue
    a, b = A[s], B.get(s, [])
    print('==', s)
    for i, la in enumerate(a):
        lb = b[i] if i < len(b) else ''
        ok = bool(lb) and key(la) == key(lb)
        nall += 1; nok += ok
        print('  %-4s archived: %-60s new: %s' % ('same' if ok else 'DIFF', la, lb))
print('%d of %d lines agree apart from the timings' % (nok, nall))
