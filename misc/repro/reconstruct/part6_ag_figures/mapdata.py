# Comparison tool for the item ag_figures (not part of the authors' package): maps the curves that pdfpaths.py read
# from fig/ag_almost.pdf and fig/ag_cancel.pdf (old or regenerated) to data coordinates, through the grid lines of each
# panel and the tick values printed on the figure, and writes <name>.data.json with the axis limits and every curve.
# Usage: python3 mapdata.py name.paths.json [...]; the name must contain ag_almost or ag_cancel.
import sys, json
import numpy as np

TICKS = {
    'ag_almost': [dict(x=[0.5, 0.6, 0.7, 0.8, 0.9, 1.0], y=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0]),
                  dict(x=[285, 287.5, 290, 292.5, 295, 297.5, 300], y=[-1, 0, 1, 2, 3]),
                  dict(x=[285, 287.5, 290, 292.5, 295, 297.5, 300], y=[-5, -4, -3, -2, -1, 0, 1, 2])],
    'ag_cancel': [dict(x=[0.5, 0.6, 0.7, 0.8, 0.9], y=[0, 2, 4, 6, 8, 10, 12, 14]),
                  dict(x=[-4, -2, 0, 2, 4], y=[-2, -1, 0, 1, 2, 3, 4]),
                  dict(x=[-4, -2, 0, 2, 4], y=[-2, -1, 0, 1, 2, 3, 4])],
}
GRID = (0.6902, 0.6902, 0.6902)


def hexcol(c):
    return '#%02x%02x%02x' % tuple(int(round(255 * v)) for v in c)


def run(fn):
    key = 'ag_almost' if 'ag_almost' in fn else 'ag_cancel'
    R = json.load(open(fn))
    clips = sorted({tuple(round(c, 2) for c in s['clip']) for s in R['strokes'] if s['clip']}, key=lambda c: c[0])
    out = dict(panels=[])
    for ip, cb in enumerate(clips):
        inpanel = [s for s in R['strokes'] if s['clip'] and tuple(round(c, 2) for c in s['clip']) == cb]
        grid = [s for s in inpanel if tuple(round(c, 4) for c in s['color']) == GRID and s['alpha'] == 0.25]
        gx = sorted(s['verts'][0][0] for s in grid if abs(s['verts'][0][0] - s['verts'][1][0]) < 1e-6)
        gy = sorted(s['verts'][0][1] for s in grid if abs(s['verts'][0][1] - s['verts'][1][1]) < 1e-6)
        tk = TICKS[key][ip]
        assert len(gx) == len(tk['x']) and len(gy) == len(tk['y']), (fn, ip, len(gx), len(gy))
        ax_, bx_ = np.polyfit(gx, tk['x'], 1)          # data x = ax_ * X + bx_
        ay_, by_ = np.polyfit(gy, tk['y'], 1)
        resx = np.max(np.abs(ax_ * np.array(gx) + bx_ - np.array(tk['x'])))
        resy = np.max(np.abs(ay_ * np.array(gy) + by_ - np.array(tk['y'])))
        X = lambda v: ax_ * v + bx_
        Y = lambda v: ay_ * v + by_
        pan = dict(xlim=[X(cb[0]), X(cb[2])], ylim=[Y(cb[1]), Y(cb[3])], curves=[], markers={})
        print('== %s panel %d: xlim [%.4f, %.4f] ylim [%.4f, %.4f] (tick fit residuals %.1e, %.1e; %.5f data units per point in y)'
              % (fn, ip, pan['xlim'][0], pan['xlim'][1], pan['ylim'][0], pan['ylim'][1], resx, resy, ay_))
        for s in inpanel:
            if s in grid:
                continue
            v = np.array(s['verts'])
            x, y = X(v[:, 0]), Y(v[:, 1])
            c = dict(color=hexcol(s['color']), lw=s['lw'], dash=s['dash'], x=x.tolist(), y=y.tolist())
            pan['curves'].append(c)
            dx = np.diff(x)
            i = int(np.argmin(y))
            print('   %s lw %.2f dash %s: %d vertices, x [%.4f, %.4f], y [%.4f, %.4f], min at x = %.4f, median step %.4f'
                  % (c['color'], s['lw'], s['dash'], len(x), x.min(), x.max(), y.min(), y.max(), x[i], np.median(dx) if len(dx) else 0))
        for d in R['dos']:
            if d['name'].startswith('M') and d['clip'] and tuple(round(c, 2) for c in d['clip']) == cb:
                pan['markers'].setdefault(hexcol(d['fill']), []).append([X(d['x']), Y(d['y'])])
        for col, pts in pan['markers'].items():
            print('   markers %s: %s' % (col, ' '.join('(%.3f, %.3f)' % tuple(p) for p in pts)))
        out['panels'].append(pan)
    json.dump(out, open(fn.replace('.paths.json', '.data.json'), 'w'))


if __name__ == '__main__':
    for fn in sys.argv[1:]:
        run(fn)
