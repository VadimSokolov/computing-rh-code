# Comparison tool for the item bp1_figures (not part of the authors' package): compares two matplotlib PDFs panel by panel,
# grid lines, curves and markers, after mapping each panel of the old file onto the same panel of the new file through the
# panels' clip boxes (the axes rectangles), so that fig/bp1_wiener.pdf and fig/bp1_salem_riesz.pdf can be compared
# numerically with the ones drawn by basepoint_one/wiener/figs_wiener.py. When the grid lines of the two panels coincide
# the two panels have the same axis limits, and the distances printed for the curves and markers are then distances between
# the plotted data, in points of the new figure (the PDF path simplification of matplotlib moves a curve by at most about
# 0.1 point, and a line is 0.8 to 1.5 points wide).
# It reads the PDFs with pdfpaths.py of the item ag_figures (misc/repro/reconstruct/part6_ag_figures/pdfpaths.py), which must
# sit next to it. Usage: python3 compare_paths.py old1.pdf new1.pdf [old2.pdf new2.pdf ...]
import sys, numpy as np
import pdfpaths as P


def load(fn):
    data = open(fn, 'rb').read()
    ss = P.pdf_streams(data)
    page = max((raw for d, raw in ss.values() if b'/Subtype' not in d), key=lambda r: r.count(b' l\n') + r.count(b' re'))
    return P.parse(page, P.ext_alpha(data), P.xobject_names(data))


def hexcol(c):
    return '#%02x%02x%02x' % tuple(int(round(255 * q)) for q in c)


def key(clip):
    return tuple(round(q, 3) for q in clip)


def segments(s):
    """Line segments of a stroked path (a Bezier piece is replaced by its chord; the data curves have none)."""
    out, it, cur, start = [], iter(s['verts']), None, None
    for c in s['codes']:
        if c == 'h':
            if cur is not None and start is not None:
                out.append((cur, start))
                cur = start
            continue
        p = next(it)
        if c == 'm':
            cur = start = p
        else:
            if cur is not None:
                out.append((cur, p))
            cur = p
    return out


def dist_to_segments(pts, segs):
    """Distance from each point to the nearest segment."""
    pts = np.asarray(pts, float)
    if len(segs) == 0:
        return np.full(len(pts), np.inf)
    A = np.array([s[0] for s in segs], float)
    B = np.array([s[1] for s in segs], float)
    D = B - A
    L2 = np.maximum((D ** 2).sum(1), 1e-300)
    best = np.empty(len(pts))
    for i in range(0, len(pts), 400):
        p = pts[i:i + 400, None, :]
        t = np.clip(((p - A) * D).sum(2) / L2, 0, 1)
        q = A + t[..., None] * D
        best[i:i + 400] = np.sqrt(((p - q) ** 2).sum(2)).min(1)
    return best


def panels(res):
    """Clip boxes that hold grid lines, i.e. the axes rectangles, from left to right."""
    boxes = {}
    for s in res['strokes']:
        if s['clip'] and abs(s['alpha'] - 0.25) < 1e-9:
            boxes[key(s['clip'])] = s['clip']
    return [boxes[k] for k in sorted(boxes, key=lambda k: (k[0], -k[1]))]


def compare(old_fn, new_fn):
    O, N = load(old_fn), load(new_fn)
    po, pn = panels(O), panels(N)
    print('== %s (old) against %s (new): %d and %d panels' % (old_fn, new_fn, len(po), len(pn)))
    worst = 0.0
    for k, (bo, bn) in enumerate(zip(po, pn)):
        sx, sy = (bn[2] - bn[0]) / (bo[2] - bo[0]), (bn[3] - bn[1]) / (bo[3] - bo[1])
        tr = lambda p: [bn[0] + (p[0] - bo[0]) * sx, bn[1] + (p[1] - bo[1]) * sy]
        print(' panel %d: old axes %s, new axes %s (points)' % (k + 1, [round(q, 2) for q in bo], [round(q, 2) for q in bn]))
        so = [s for s in O['strokes'] if s['clip'] and key(s['clip']) == key(bo)]
        sn = [s for s in N['strokes'] if s['clip'] and key(s['clip']) == key(bn)]
        # grid lines: position of each vertical and horizontal line relative to the axes, old mapped onto new
        for name, axis in (('vertical', 0), ('horizontal', 1)):
            go = sorted(tr(s['verts'][0])[axis] for s in so if abs(s['alpha'] - 0.25) < 1e-9 and abs(s['verts'][0][axis] - s['verts'][-1][axis]) < 1e-6)
            gn = sorted(s['verts'][0][axis] for s in sn if abs(s['alpha'] - 0.25) < 1e-9 and abs(s['verts'][0][axis] - s['verts'][-1][axis]) < 1e-6)
            d = max(abs(a - b) for a, b in zip(go, gn)) if len(go) == len(gn) and go else float('nan')
            print('   %-10s grid lines: old %d, new %d, largest shift %.2e points' % (name, len(go), len(gn), d))
            worst = max(worst, d if d == d else np.inf)
        # curves, grouped by colour
        cols = sorted({tuple(round(q, 4) for q in s['color']) for s in so + sn if abs(s['alpha'] - 0.25) > 1e-9})
        for c in cols:
            ao = [s for s in so if tuple(round(q, 4) for q in s['color']) == c and abs(s['alpha'] - 0.25) > 1e-9]
            an = [s for s in sn if tuple(round(q, 4) for q in s['color']) == c and abs(s['alpha'] - 0.25) > 1e-9]
            vo = [tr(p) for s in ao for p in s['verts']]
            vn = [p for s in an for p in s['verts']]
            go_ = [(tr(a), tr(b)) for s in ao for a, b in segments(s)]
            gn_ = [sg for s in an for sg in segments(s)]
            d1 = dist_to_segments(vo, gn_).max() if vo else float('nan')
            d2 = dist_to_segments(vn, go_).max() if vn else float('nan')
            lw = sorted({(round(s['lw'], 3), str(s['dash'][0])) for s in ao}), sorted({(round(s['lw'], 3), str(s['dash'][0])) for s in an})
            print('   curve %s: old %d paths %d vertices, new %d paths %d vertices; old to new %.3f, new to old %.3f points; width, dash old %s new %s'
                  % (hexcol(c), len(ao), len(vo), len(an), len(vn), d1, d2, lw[0], lw[1]))
            worst = max(worst, d1 if d1 == d1 else np.inf, d2 if d2 == d2 else np.inf)
        # markers (XObjects placed inside the axes), grouped by fill colour
        mo = [d for d in O['dos'] if d['clip'] and key(d['clip']) == key(bo)]
        mn = [d for d in N['dos'] if d['clip'] and key(d['clip']) == key(bn)]
        for c in sorted({tuple(round(q, 4) for q in d['fill']) for d in mo + mn}):
            xo = np.array([tr((d['x'], d['y'])) for d in mo if tuple(round(q, 4) for q in d['fill']) == c]).reshape(-1, 2)
            xn = np.array([(d['x'], d['y']) for d in mn if tuple(round(q, 4) for q in d['fill']) == c]).reshape(-1, 2)
            if len(xo) and len(xn):
                dd = np.sqrt(((xo[:, None, :] - xn[None, :, :]) ** 2).sum(2))
                d1, d2 = dd.min(1).max(), dd.min(0).max()
            else:
                d1 = d2 = np.inf
            print('   markers %s: old %d, new %d; old to new %.3f, new to old %.3f points' % (hexcol(c), len(xo), len(xn), d1, d2))
            worst = max(worst, d1, d2)
    print(' largest distance over all grid lines, curves and markers: %.3f points' % worst)
    return worst


if __name__ == '__main__':
    a = sys.argv[1:]
    for i in range(0, len(a), 2):
        compare(a[i], a[i + 1])
