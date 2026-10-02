# Check helper of unit R2-23 (not an archive script): compare the joined block outputs of the edited
# code_bosch/verify_bosch_xi_v7.py with (1) the joined output of the reconstruction's fixed copy (job 1242945, task 0)
# and (2) the archived code_bosch/verify_bosch_xi_v7_output.txt, line by line; timings "(123s)" are masked.
# Usage: python3 cmp_bosch.py NEW_JOINED.txt RECONSTRUCT_FIXED_JOINED.txt ARCHIVED_OUTPUT.txt
import re, sys

mask = lambda s: re.sub(r'\(\d+s\)', '(Ns)', s)
new = open(sys.argv[1]).read().splitlines()
rec = open(sys.argv[2]).read().splitlines()
arc = open(sys.argv[3]).read().splitlines()
print('lines: new %d, reconstruction fixed copy %d, archive %d' % (len(new), len(rec), len(arc)))
print('new == fixed copy, byte for byte:', new == rec)
d = [i for i in range(max(len(new), len(rec))) if i >= len(new) or i >= len(rec) or new[i] != rec[i]]
dm = [i for i in d if i >= len(new) or i >= len(rec) or mask(new[i]) != mask(rec[i])]
print('lines that differ from the fixed copy: %d, of which %d apart from timings' % (len(d), len(dm)))
for i in d:
    print('  line %d\n    new:   %s\n    fixed: %s' % (i + 1, new[i] if i < len(new) else '', rec[i] if i < len(rec) else ''))
# alignment with the archive: walk through the archive in order, matching each new line to the next equal archive line
j, matched, exact, unmatched = 0, [], 0, []
for i, s in enumerate(new):
    k = j
    while k < len(arc) and mask(arc[k]) != mask(s):
        k += 1
    if k < len(arc):
        matched.append((i + 1, k + 1)); exact += (arc[k] == s); j = k + 1
    else:
        unmatched.append(i + 1)
skipped = sorted(set(range(1, len(arc) + 1)) - {k for _, k in matched})
print('new lines found in order in the archive: %d of %d (%d byte for byte, %d only apart from timings); '
      'new lines not found: %s' % (len(matched), len(new), exact, len(matched) - exact, unmatched))
print('archive lines with no counterpart in the new run: %s' % skipped)
for i, k in matched:
    if arc[k - 1] != new[i - 1]:
        print('  timing differs: new line %d / archive line %d\n    new:     %s\n    archive: %s' % (i, k, new[i - 1], arc[k - 1]))
