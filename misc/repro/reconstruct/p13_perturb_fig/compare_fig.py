# Helper of the reconstruction in misc/repro/reconstruct/p13_perturb_fig/ (not the authors' code): compare the curve
# of the archived fig/perturb.pdf (kept here as old_perturb.pdf) with the curves that perturb.py computes. Each PDF is
# read in its own data units: its seven vertical grid lines are t = 1e-8, ..., 1e-2 and its seven horizontal grid
# lines are log10 = -3000, ..., 0. Every vertex that matplotlib kept in the stroked curve is mapped back to the index
# i of the grid t = logspace(-8, -2, 400) and to its plotted value, which is compared with perturb.json.
# Run: python3 compare_fig.py old_perturb.pdf perturb.json [candidate.pdf ...]   (writes compare_fig.json)
import json, sys
import numpy as np
from pdf_paths import streams, paths


def curve(fname):
    P = paths(max(streams(fname), key=len))
    ic = max(range(len(P)), key=lambda k: P[k]['n'])
    vert = sorted(p['pts'][0][0] for p in P[:ic] if p['n'] == 2 and p['pts'][0][0] == p['pts'][1][0])
    horz = sorted(p['pts'][0][1] for p in P[:ic] if p['n'] == 2 and p['pts'][0][1] == p['pts'][1][1])
    assert len(vert) == 7 and len(horz) == 7, (vert, horz)
    X0, X6, Y0, Y6 = vert[0], vert[-1], horz[0], horz[-1]
    pts = P[ic]['pts']
    lt = [-8 + 6 * (x - X0) / (X6 - X0) for x, y in pts]
    val = [-3000 + 3000 * (y - Y0) / (Y6 - Y0) for x, y in pts]
    return dict(n=len(pts), lw=P[ic]['lw'], rgb=P[ic]['rgb'], grid_x=vert, grid_y=horz, pts=pts, log10t=lt, value=val,
                axes_x=(min(p['pts'][0][0] for p in P if p['n'] == 2), max(p['pts'][1][0] for p in P if p['n'] == 2)))


if __name__ == '__main__':
    old = curve(sys.argv[1]); R = json.load(open(sys.argv[2]))
    Y = {k: np.array(v) for k, v in R['log10_rel'].items()}
    out = dict(archived=dict(n=old['n'], lw=old['lw'], rgb=old['rgb'], grid_x=old['grid_x'], grid_y=old['grid_y']))
    rows = []
    for j, (lt, v, p) in enumerate(zip(old['log10t'], old['value'], old['pts'])):
        fi = (lt + 8) * 399 / 6; i = int(round(fi))
        row = dict(vertex=j, pdf=p, index=fi, value=v)
        if abs(fi - i) < 1e-3 and p[1] > 0:
            row.update(t=R['t'][i], **{'diff_' + k: v - Y[k][i] for k in Y})
        rows.append(row)
    out['archived_vertices'] = rows
    kept = [r for r in rows if 'diff_w1700' in r]
    out['summary'] = {k: dict(max_abs_diff=max(abs(r['diff_' + k]) for r in kept),
                              first_vertex_diff=kept[0]['diff_' + k]) for k in Y}
    out['summary']['vertices_on_grid'] = len(kept)
    out['summary']['vertices_total'] = len(rows)
    print('archived curve: %d vertices, %d on the grid; lw %s rgb %s' % (len(rows), len(kept), old['lw'], old['rgb']))
    for k in Y:
        print('  W = %-9s max |archived - computed| = %.3g (log10 units); at t = 1e-8: %.6f' % (
            k, out['summary'][k]['max_abs_diff'], out['summary'][k]['first_vertex_diff']))
    for r in rows:
        print('  vertex %2d  pdf (%.6f, %.6f)  i = %9.4f  value %12.5f  ' % (r['vertex'], r['pdf'][0], r['pdf'][1],
              r['index'], r['value']) + ' '.join('%s %+.2e' % (k, r['diff_' + k]) for k in Y if 'diff_' + k in r))
    old_idx = [int(round(r['index'])) for r in kept]
    old_stream = max(streams(sys.argv[1]), key=len)
    for f in sys.argv[3:]:
        c = curve(f)
        idx = [(lt + 8) * 399 / 6 for lt in c['log10t']]
        on = [k for k, x in enumerate(idx) if abs(x - round(x)) < 1e-3 and c['pts'][k][1] > 0]
        new_idx = [int(round(idx[k])) for k in on]
        common = sorted(set(new_idx) & set(old_idx))
        vo = {int(round(r['index'])): r['value'] for r in kept}; vn = {int(round(idx[k])): c['value'][k] for k in on}
        new_stream = max(streams(f), key=len)
        out[f] = dict(n=c['n'], lw=c['lw'], rgb=c['rgb'], grid_x=c['grid_x'], grid_y=c['grid_y'],
                      same_grid_lines=(c['grid_x'] == old['grid_x'] and c['grid_y'] == old['grid_y']),
                      vertex_indices=new_idx, same_vertex_indices=(new_idx == old_idx),
                      only_in_archived=sorted(set(old_idx) - set(new_idx)), only_in_candidate=sorted(set(new_idx) - set(old_idx)),
                      max_value_diff_common=max(abs(vn[i] - vo[i]) for i in common),
                      max_pdf_coordinate_diff=(max(max(abs(a[0] - b[0]), abs(a[1] - b[1])) for a, b in zip(c['pts'], old['pts']))
                                               if c['n'] == old['n'] else None),
                      clipped_end=c['pts'][-1], clipped_end_index=idx[-1],
                      content_stream_identical=(new_stream == old_stream),
                      content_stream_lines_differing=sum(a != b for a, b in zip(new_stream.splitlines(), old_stream.splitlines()))
                      + abs(len(new_stream.splitlines()) - len(old_stream.splitlines())))
        print('==', f, json.dumps(out[f]))
    json.dump(out, open('compare_fig.json', 'w'), indent=1)
