# vertex_check.py (round 2 item kernel_fig, not part of the authors' package). A detail of the comparison made by
# pdf_compare.py: for each curve of the right panel it lists, in data coordinates, the path vertices of the archived
# fig/kernel.pdf (old_kernel.pdf) and of the new kernel.pdf that lie near the bottom edge 1e-12 of the axes or off the
# grid of 601 points, to explain why the number of compared vertices of one curve differs by one between the two files.
import json, numpy as np
from pdf_compare import page, boxes

J = json.load(open('kernel_curves.json'))
x = np.array(J['x']); ax1 = J['axes'][1]; ly = np.log10(ax1['ylim'])
for f in ['old_kernel.pdf', 'kernel.pdf']:
    P = page(f); b = boxes(P)[1]
    print(f, 'right box', b)
    for col, (d, c) in zip(['#1b4f72', '#c0392b', '#27ae60', '#d68910'], J['curves'].items()):
        sub = [s['pts'] for s in P['strokes'] if s['lw'] == 1.4 and s['col'] == col and len(s['pts']) > 3]
        pts = np.array([p for q in sub for p in q])
        X = ax1['xlim'][0] + (pts[:, 0] - b[0]) / (b[2] - b[0]) * (ax1['xlim'][1] - ax1['xlim'][0])
        LY = ly[0] + (pts[:, 1] - b[1]) / (b[3] - b[1]) * (ly[1] - ly[0])
        k = np.clip(np.rint((X - x[0]) / (x[1] - x[0])).astype(int), 0, len(x) - 1)
        ref = np.log10(np.maximum(np.array(c['v']), 1e-30))
        rows = [(round(float(X[i]), 6), round(float(LY[i]), 4), round(float(pts[i, 1] - b[1]), 4), round(float(ref[k[i]]), 4), bool(abs(X[i] - x[k[i]]) < 1e-6))
                for i in range(len(X)) if LY[i] < -11.9 or abs(X[i] - x[k[i]]) >= 1e-6]
        print('  d=%s: %d vertices; (x, log10 y from the PDF, height above the box bottom in pt, log10 of the computed value, on grid):' % (d, len(X)))
        for r in rows: print('    ', r)
