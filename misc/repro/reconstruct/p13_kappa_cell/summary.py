# Collector of this folder, not a book script: reads kappa_cell.json, interp_check.json and extra_patched.json
# (this folder) and the archived data/spacing.json and data/analysis.json (/scratch/vsokolov/rh_book_repro/book_data,
# read only), and prints every number that REPORT.md quotes, rounded as the book would print it. Output: summary.txt.
import json
import numpy as np
BOOK = '/scratch/vsokolov/rh_book_repro/book_data/'
R = json.load(open('kappa_cell.json')); I = json.load(open('interp_check.json')); X = json.load(open('extra_patched.json'))
SP = json.load(open(BOOK + 'spacing.json')); AN = json.load(open(BOOK + 'analysis.json'))
A = R['A_exact_cells']; H = R['C_truncated_hp']; D = R['D_finite_eps']['rows']
out = []
p = lambda *a: out.append(' '.join(str(x) for x in a))
k = A['kappa_cell_float']
p('kappa_cell (Arb, exact cells):', A['kappa_cell'])
p('  to 6 significant digits %.6g, to 8 %.8g; float64 cells of analyse.py %s; Polya route block of plots.py %.16g' % (
    k, k, R['A2_float64_cells']['kappa_cell'], X['kappa_cell']))
p('truncated (1700 zeros): spacing.json coef %.16g, high precision %s; kappa - kappa_1700 = %s = %.4f%% of kappa, %.4f%% of kappa_1700' % (
    SP['coef'], H['kappa_1700'], H['kappa_minus_kappa_1700'], H['percent_of_kappa'], H['percent_of_kappa_1700']))
p('kappa_1700 + smooth tail beyond T* = %s: %s (kappa minus it %s)' % (H['Tstar'], H['kappa_1700_plus_tail'],
                                                                      H['kappa_minus_kappa_1700_plus_tail']))
ex, tr = np.array(A['cell_errors']), np.array(SP['leaks'])
b = np.array([float(x) for x in A['boundaries']]); w = np.diff(b)
dd = ex - tr
p('exact minus truncated cell errors: min %.3e, max %.3e, mean %.3e; (exact - truncated)/|C_k| from %.4e to %.4e' % (
    dd.min(), dd.max(), dd.mean(), (dd / w).min(), (dd / w).max()))
p('  (2/pi) sum_{j>1700} gamma_j^-2 by the smooth density: %.4e' % (
    2 / np.pi * (np.log(2197.88 / (2 * np.pi)) + 1) / (2 * np.pi * 2197.88)))
p('sum of cell errors: exact %s = -l(200)/pi %s; truncated %.7f; extra.json pred_coeff %.17g' % (
    A['sum_cell_errors'], A['minus_l200_over_pi'], float(np.sum(tr)), X['pred_coeff']))
p('range: exact %.5f (cell %d) to %.5f (cell %d); truncated %.5f (cell %d) to %.5f (cell %d)' % (
    A['min'], A['argmin'], A['max'], A['argmax'], tr.min(), tr.argmin() + 1, tr.max(), tr.argmax() + 1))
p('largest positive, exact:', A['largest_five'][:3], '; truncated:', R['C_truncated_float64']['largest_five'][:3])
p('eps^2 coefficient of d_cell/eps: %s' % A['c2_eps2_coefficient'])
for r in D:
    s = '  eps %-6g exact d_cell %.10f slope %.7f slope-kappa %+.3e cells [%.5f, %.5f] mass %.5f' % (
        r['eps'], r['d_cell'], r['slope'], r['slope_minus_kappa'], r['cell_min'], r['cell_max'], r['mass'])
    if 'archived_tv_cells' in r:
        s += ' | archived %.10f slope %.7f cells [%.5f, %.5f]' % (r['archived_tv_cells'], r['archived_slope'],
                                                                 r['archived_cell_min'], r['archived_cell_max'])
    p(s)
p('Richardson of the archived slopes at 0.01 and 0.005: %.9f' % R['D_finite_eps']['richardson_archived'])
printed = {1.0: '4.995', 0.5: '2.656', 0.25: '1.351', 0.1: '0.5431', 0.05: '0.2718', 0.02: '0.1087', 0.01: '0.05437',
           0.005: '0.02718'}
p('tab:tv cell TV: printed | archived analysis.json | exact | Hermite interpolation of the same grid')
for r in I['rows']:
    e = r['eps']; nd = len(printed[e].split('.')[1])
    fmt = '%%.%df' % nd
    p('  eps %-6g %-8s | %s | %s | %s   (grid step %.5f, linear interpolation error at the cells up to %.1e)' % (
        e, printed[e], fmt % r['tv_cells_archived'], fmt % r['tv_cells_exact'], fmt % r['tv_cells_hermite'],
        r['grid_step'], r['max_abs_interp_error_at_cells']))
open('summary.txt', 'w').write('\n'.join(out) + '\n')
print('\n'.join(out))
