# Helper of the reconstruction in misc/repro/reconstruct/spacing/ (not the authors' code): compare the histogram of
# the book's fig/spacing.pdf with the candidate drawn by spacing.py and with the exact histogram of the reconstructed
# spacings (spacing_hist.json). Each figure is read in its own data units: x = 0 and x = 3 at the ends of the two
# curves, y = 0 at the start of the surmise (0 at s = 0) and y = 1 at the start of the Poisson curve (e^0 = 1).
# Run: python3 fig_check.py BOOK.pdf CANDIDATE.pdf spacing_hist.json   (writes fig_check.json)
import json, sys
import numpy as np
from pdf_bars import parse


def bars(fname):
    r = parse(fname)
    curves = sorted([p for p in r['paths'] if p['op'] == 'S'], key=lambda p: -p['n'])[:2]
    gue = min(curves, key=lambda p: p['first'][1]); poi = max(curves, key=lambda p: p['first'][1])
    x0, x3 = gue['first'][0], gue['last'][0]
    y0, y1 = gue['first'][1], poi['first'][1]
    fx = lambda X: (X - x0) / (x3 - x0) * 3; fy = lambda Y: (Y - y0) / (y1 - y0)
    B = sorted([q for q in r['rects'] if q[4] == 'B' and abs(q[1] - y0) < 1e-6], key=lambda q: q[0])
    return dict(left=[fx(q[0]) for q in B], right=[fx(q[2]) for q in B], height=[fy(q[3]) for q in B])


if __name__ == '__main__':
    book, cand = bars(sys.argv[1]), bars(sys.argv[2])
    H = json.load(open(sys.argv[3])); e = np.array(H['edges']); d = np.array(H['density'])
    out = {}
    for name, f in [('book', book), ('candidate', cand)]:
        n = len(f['left'])
        pos = np.abs(d) > 0
        out[name] = dict(bars=n, bins_in_hist=len(d), nonzero_bins=int(pos.sum()))
        if n == int(pos.sum()):
            out[name].update(
                max_edge_diff=float(max(np.max(np.abs(np.array(f['left']) - e[:-1][pos])),
                                        np.max(np.abs(np.array(f['right']) - e[1:][pos])))),
                max_height_diff=float(np.max(np.abs(np.array(f['height']) - d[pos]))),
                max_rel_height_diff=float(np.max(np.abs(np.array(f['height']) - d[pos]) / d[pos])))
    out['first_edges_book'] = book['left'][:3]; out['first_edges_hist'] = e[:3].tolist()
    out['heights_book_first5'] = book['height'][:5]; out['density_first5'] = d[:5].tolist()
    json.dump(out, open('fig_check.json', 'w'), indent=1)
    print(json.dumps(out, indent=1))
