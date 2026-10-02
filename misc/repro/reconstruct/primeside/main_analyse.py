"""RECONSTRUCTION AID, NOT THE AUTHORS' CODE.

Analysis stage of the attempt to regenerate data/results_main.json (see main_compute.py).
For each field of the archived file it evaluates candidate definitions on the arrays saved
by main_compute.py and reports, per candidate, the largest relative difference from the
archived values over the six resolutions. A candidate that agrees to 1e-9 or better at every
resolution is taken as the definition the lost script used.
Usage: python main_analyse.py BOOK_results_main.json   (reads main_EPS.npz; writes main_analyse.json
       and main_reconstructed.json, the file rebuilt from the best candidates)
"""
import json, sys
import numpy as np
from scipy import integrate

book = json.load(open(sys.argv[1]))
EPS = [r['eps'] for r in book['rows']]
D = {e: np.load(f'main_{e}.npz') for e in EPS}
trap = integrate.trapezoid


def rel(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(np.max(np.abs(a - b) / np.maximum(np.maximum(np.abs(a), np.abs(b)), 1e-300)))


report = {}
d0 = D[EPS[0]]
gw = d0['gw']
top = dict(N_ok=(book['N'] == 1700), T_rel=rel(book['T'], float(d0['g1700'])), T_bitwise=(book['T'] == float(d0['g1700'])),
           B_ok=(book['B'] == 58.0), gamma_window_rel=rel(book['gamma_window'], gw),
           gamma_window_bitwise=sum(1 for x, y in zip(book['gamma_window'], gw) if x == float(y)), n_window=len(gw))
report['top'] = top
print('top level:', top)

cands = {}  # field -> name -> list over eps


def add(field, name, e, val):
    cands.setdefault(field, {}).setdefault(name, {})[e] = val


for e in EPS:
    d = D[e]; th = d['th']; F = d['F']; rho = F / np.pi; arg = d['arg']
    add('npts', 'len', e, len(th))
    add('minrho', 'min(F)/pi', e, float(F.min() / np.pi))
    add('minrho', 'min(F/pi)', e, float(rho.min()))
    add('mass', 'trapezoid', e, float(trap(rho, th)))
    add('mass', 'simpson', e, float(integrate.simpson(rho, x=th)))
    add('mass', 'argument', e, float(arg[-1] / np.pi))
    for rn in ['Rld', 'R64']:
        R = d[rn]
        add('tv_raw', f'trapezoid_{rn}', e, float(trap(np.abs(rho - R), th)))
        add('tv_raw', f'simpson_{rn}', e, float(integrate.simpson(np.abs(rho - R), x=th)))
        for tn in [k for k in d.files if k.startswith('tail_')]:
            diff = rho - R - d[tn]
            add('tv', f'trapezoid_{rn}_{tn}', e, float(trap(np.abs(diff), th)))
            add('sup', f'max_{rn}_{tn}', e, float(np.abs(diff).max()))
    mids = [(gw[k] + gw[k + 1]) / 2 for k in range(len(gw) - 1)]
    for lastname, last in [('B', 58.0), ('mid13', (gw[-1] + float(d['gnext'])) / 2)]:
        edges = np.array([0.0] + mids + [last])
        C = integrate.cumulative_trapezoid(rho, th, initial=0.0)
        M = arg / np.pi
        for how in ['mask_lt', 'mask_le', 'mask_gt', 'interp_trapz', 'interp_arg']:
            am = []
            for a, b in zip(edges[:-1], edges[1:]):
                if how == 'mask_lt':
                    m = (th >= a) & (th < b); am.append(float(trap(rho[m], th[m])))
                elif how == 'mask_le':
                    m = (th >= a) & (th <= b); am.append(float(trap(rho[m], th[m])))
                elif how == 'mask_gt':
                    m = (th > a) & (th <= b); am.append(float(trap(rho[m], th[m])))
            if how == 'interp_trapz':
                am = list(np.diff(np.interp(edges, th, C)))
            if how == 'interp_arg':
                am = list(np.diff(np.interp(edges, th, M)))
            add('atom_masses', f'{how}_last{lastname}', e, [float(x) for x in am])

rows = {r['eps']: r for r in book['rows']}
best = {}
print()
for field, cs in cands.items():
    res = []
    for name, vals in cs.items():
        r = max(rel(vals[e], rows[e][field]) for e in EPS)
        per = {e: rel(vals[e], rows[e][field]) for e in EPS}
        nb = sum(1 for e in EPS if np.array_equal(np.asarray(vals[e], float), np.asarray(rows[e][field], float)))
        res.append((r, name, per, nb))
    res.sort(key=lambda x: x[0])
    best[field] = res[0][1]
    report[field] = [dict(candidate=n, max_rel=r, per_eps=p, bitwise_rows=nb) for r, n, p, nb in res]
    print(f'field {field}: {len(res)} candidates; best:')
    for r, n, p, nb in res[:4]:
        print(f'   {n:40s} max rel {r:.3e}  ' + ' '.join(f'{p[e]:.1e}' for e in EPS) + f'  ({nb}/{len(EPS)} rows bitwise)')
# the file rebuilt from the best candidate of each field (secs: our Polya timings, not comparable)
out = dict(N=1700, T=float(d0['g1700']), B=58.0, gamma_window=[float(x) for x in gw], rows=[])
for e in EPS:
    row = dict(eps=e)
    for field in ['npts', 'tv_raw', 'tv', 'sup', 'mass', 'minrho', 'atom_masses']:
        row[field] = cands[field][best[field]][e]
    row['secs'] = float(D[e]['secs_polya'])
    out['rows'].append(row)
json.dump(out, open('main_reconstructed.json', 'w'), indent=1)
report['best'] = best
json.dump(report, open('main_analyse.json', 'w'), indent=1)
print('\nbest candidates:', best)
