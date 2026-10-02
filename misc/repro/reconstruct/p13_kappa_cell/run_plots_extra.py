# Driver of this folder, not a book script: copies plotcurves.npz and results_prime.json of the earlier rerun
# (/scratch/vsokolov/rh_book_repro/run, read only) beside plots_extra_patched.py, runs it, and compares the
# extra.json it writes with the archived data/extra.json (/scratch/vsokolov/rh_book_repro/book_data) and with the
# certified values of kappa_cell.json (Arb, this folder). Output: extra_patched.json, compare_extra.txt.
import json, shutil, subprocess, sys
RUN = '/scratch/vsokolov/rh_book_repro/run/'
BOOK = '/scratch/vsokolov/rh_book_repro/book_data/'
for f in ['plotcurves.npz', 'results_prime.json']:
    shutil.copy(RUN + f, f)
subprocess.run([sys.executable, '-u', 'plots_extra_patched.py'], check=True)
shutil.copy('extra.json', 'extra_patched.json')
new, old = json.load(open('extra_patched.json')), json.load(open(BOOK + 'extra.json'))
ref = json.load(open('kappa_cell.json'))['A_exact_cells']
lines = ['keys archived: %s; keys patched: %s' % (list(old), list(new))]
same = [k for k in old if new.get(k) == old[k]]
lines.append('archived keys with identical values: %s of %d (%s)' % (len(same), len(old), ', '.join(same)))
lines.append('kappa_cell patched %.16g; Arb at the exact cells %s (difference %.2e, from the float64 cells)' % (
    new['kappa_cell'], ref['kappa_cell'], new['kappa_cell'] - ref['kappa_cell_float']))
d = max(abs(x - y) for x, y in zip(new['cell_leaks'], ref['cell_errors']))
lines.append('cell_leaks: %d values, max |patched - Arb| = %.2e, sum %.16g, min %.6f (cell %d), max %.6f (cell %d)' % (
    len(new['cell_leaks']), d, sum(new['cell_leaks']), min(new['cell_leaks']),
    new['cell_leaks'].index(min(new['cell_leaks'])) + 1, max(new['cell_leaks']),
    new['cell_leaks'].index(max(new['cell_leaks'])) + 1))
open('compare_extra.txt', 'w').write('\n'.join(lines) + '\n')
print('\n'.join(lines))
