# Comparison tool for the item ag_figures (not part of the authors' package): compares the old fig/ag_almost.pdf and
# fig/ag_cancel.pdf with the ones drawn by riesz_program/code/figs_almost.py. For both versions it reads the vector paths
# (pdfpaths.py), maps them to data coordinates through the grid lines (mapdata.py), and measures every vertex of every
# curve against the exact values from the archived data (the formula of almost.py with g_*.npy and almost.json, more.json,
# high.json); it also compares the axis limits, colours, line widths and dash patterns panel by panel. The PNG renderings
# for the comparison by eye are made by render_png.py.
# Usage (run folder holding old_*.pdf, the new ag_*.pdf, the JSON and npy files): python3 compare_figs.py
import json
import numpy as np
import pdfpaths, mapdata

A = json.load(open('almost.json')); M = json.load(open('more.json')); H = json.load(open('high.json'))
g = np.concatenate([np.load(f) for f in ['g_1_400.npy', 'g_401_1000.npy', 'g_1001_1700.npy']])
tj = A['synthetic']['tj']; off = [(0.75, t) for t in tj]
P = lambda a, x: a / (np.pi * (a * a + x * x))


def rho(alpha, off, th):
    r = np.zeros_like(th)
    for gg in g:
        r += P(alpha - 0.5, th - gg) + P(alpha - 0.5, th + gg)
    for (b, t) in off:
        r += P(alpha - b, th - t) + P(alpha - (1 - b), th - t) + P(alpha - b, th + t) + P(alpha - (1 - b), th + t)
    return r


def exact(fig, ip, c):
    """Exact values at the vertices of curve c of panel ip, and the distance of the vertices to the data grid."""
    x = np.array(c['x']); col, lw = c['color'], c['lw']
    if fig == 'ag_almost' and ip == 0:
        f = {'#d68910': lambda s: 3 * (1 - s) / (2 - s), '#c0392b': lambda s: 15 * (1 - s) / (3 + 5 * s),
             '#808080': lambda s: 2 * (1 - s), '#000000': lambda s: np.ones_like(s)}[col]
        return f(x), 0.0, {'#d68910': 'Ingham', '#c0392b': 'Guth and Maynard', '#808080': 'density hypothesis', '#000000': 'level 1'}[col]
    if fig == 'ag_almost':
        alpha = 0.6 if ip == 1 else 0.7
        if col == '#000000':
            return (np.zeros_like(x), 0.0, 'axis 0') if lw > 0.55 else (None, float(np.max(np.abs(x - tj[4]))), 'dotted line at tj[4]')
        k = np.round((x - 284.0) / 0.005); xs = np.linspace(284.0, 300.0, 3201)[k.astype(int)]
        return rho(alpha, [] if col == '#808080' else off, xs), float(np.max(np.abs(x - xs))), ('line only' if col == '#808080' else 'with hypothetical zeros')
    if fig == 'ag_cancel' and ip == 0:
        if c['dash'][0] and len(x) == 2:
            return None, 0.0, 'reference line'
        beta = {'#c0392b': 0.6, '#d68910': 0.75, '#1b4f72': 0.9}[col]
        rows = [r for r in M['curve'] if r['beta'] == beta]
        return np.interp(x, [r['alpha'] for r in rows], [r['neg'] for r in rows]), 0.0, 'beta %s' % beta
    key = '300.0' if ip == 1 else '100000.0'
    h = H[key]; gam0 = float(key); tt = np.array(h['th']); bg = np.array(h['bg'])
    if col == '#000000':
        return np.zeros_like(x), 0.0, 'axis 0'
    k = np.round((x + gam0 - tt[0]) / (tt[1] - tt[0])).astype(int)
    grid = tt[k] - gam0
    if col == '#808080':
        y = bg[k]; name = 'background'
    else:
        beta = 0.7 if col == '#c0392b' else 0.9
        y = bg[k] + P(-(beta - 0.6), tt[k] - gam0) + P(0.6 - (1 - beta), tt[k] - gam0); name = 'beta %s' % beta
    return y, float(np.max(np.abs(x - grid))), name


report = {}
for fig in ('ag_almost', 'ag_cancel'):
    D = {}
    for ver, fn in (('old', 'old_%s.pdf' % fig), ('new', '%s.pdf' % fig)):
        pdfpaths.main(fn)
        mapdata.run(fn + '.paths.json')
        D[ver] = json.load(open(fn + '.data.json'))
    report[fig] = []
    for ip in range(3):
        po, pn = D['old']['panels'][ip], D['new']['panels'][ip]
        rec = dict(panel=ip, xlim_old=po['xlim'], xlim_new=pn['xlim'], ylim_old=po['ylim'], ylim_new=pn['ylim'],
                   styles_old=sorted('%s lw %.2f dash %s' % (c['color'], c['lw'], c['dash']) for c in po['curves']),
                   styles_new=sorted('%s lw %.2f dash %s' % (c['color'], c['lw'], c['dash']) for c in pn['curves']), curves=[])
        for ver, pan in (('old', po), ('new', pn)):
            for c in pan['curves']:
                y, dgrid, name = exact(fig, ip, c)
                err = None if y is None else float(np.max(np.abs(np.array(c['y']) - y)))
                rec['curves'].append(dict(version=ver, name=name, color=c['color'], vertices=len(c['x']),
                                          x_range=[min(c['x']), max(c['x'])], y_range=[min(c['y']), max(c['y'])],
                                          max_distance_to_data_grid=dgrid, max_abs_error_against_data=err))
            for col, pts in pan['markers'].items():
                beta = {'#c0392b': 0.6, '#d68910': 0.75, '#1b4f72': 0.9}[col]
                rows = {round(r['alpha'], 6): r['neg'] for r in M['curve'] if r['beta'] == beta}
                errs = [abs(p[1] - rows[round(p[0], 2)]) for p in pts]
                rec['curves'].append(dict(version=ver, name='markers beta %s' % beta, color=col, vertices=len(pts),
                                          max_abs_error_against_data=float(max(errs))))
        report[fig].append(rec)
        print('== %s panel %d: xlim old [%.4f, %.4f] new [%.4f, %.4f]; ylim old [%.4f, %.4f] new [%.4f, %.4f]; styles %s'
              % (fig, ip, *po['xlim'], *pn['xlim'], *po['ylim'], *pn['ylim'], 'equal' if rec['styles_old'] == rec['styles_new'] else 'DIFFER'))
        if rec['styles_old'] != rec['styles_new']:
            print('   old', rec['styles_old']); print('   new', rec['styles_new'])
        for c in rec['curves']:
            print('   %s %-26s %s %4d vertices, error against the data %s, distance to data grid %s'
                  % (c['version'], c['name'], c['color'], c['vertices'],
                     'n/a' if c['max_abs_error_against_data'] is None else '%.1e' % c['max_abs_error_against_data'],
                     '%.1e' % c.get('max_distance_to_data_grid', 0.0)))
json.dump(report, open('compare_figs.json', 'w'), indent=1)
