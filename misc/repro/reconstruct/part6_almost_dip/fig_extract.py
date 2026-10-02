# Check for the item part6_almost_dip: reads the vector paths of the authors' figure riesz_program/notes/fig_almost.pdf
# (the same file as fig/ag_almost.pdf before this review) and converts the curves of its middle and right panels to data
# coordinates through the grid lines, to find the window, the sampling step and the lowest plotted point of each curve.
# Pure Python (re, zlib); run on Hopper: bash misc/tools/hopper_run.sh -c 1 -t 10 fig_extract.py fig_almost.pdf
import re, sys, zlib, json
from collections import defaultdict

pdf = open(sys.argv[1] if len(sys.argv) > 1 else 'fig_almost.pdf', 'rb').read()
streams = []
for m in re.finditer(rb'(?<!end)stream\r?\n', pdf):
    start = m.end()
    end = pdf.find(b'endstream', start)
    raw = pdf[start:end]
    try:
        streams.append(zlib.decompressobj().decompress(raw))
    except Exception:
        streams.append(raw)

tok_re = re.compile(rb'/[^\s/\[\]()<>{}%]+|\[|\]|<<|>>|<[0-9A-Fa-f\s]*>|\((?:\\.|[^\\)])*\)|[-+]?(?:\d+\.\d*|\.\d+|\d+)|[A-Za-z\'"*]+')

def mul(a, b):
    # PDF matrices [a b c d e f]; returns a x b
    return [a[0]*b[0] + a[1]*b[2], a[0]*b[1] + a[1]*b[3], a[2]*b[0] + a[3]*b[2], a[2]*b[1] + a[3]*b[3],
            a[4]*b[0] + a[5]*b[2] + b[4], a[4]*b[1] + a[5]*b[3] + b[5]]

def app(M, x, y):
    return (M[0]*x + M[2]*y + M[4], M[1]*x + M[3]*y + M[5])

paths = []
for s in streams:
    if b' m' not in s and b'\nm' not in s:
        continue
    ops, stack = [], []
    st = dict(ctm=[1, 0, 0, 1, 0, 0], rgb=(0, 0, 0), lw=1.0, dash='', gs='')
    saved = []
    cur, curve = [], False
    for t in tok_re.findall(s):
        if re.fullmatch(rb'[-+]?(?:\d+\.\d*|\.\d+|\d+)', t):
            stack.append(float(t)); continue
        if t.startswith(b'/') or t in (b'[', b']') or t.startswith(b'<') or t.startswith(b'('):
            stack.append(t.decode('latin1')); continue
        op = t.decode('latin1')
        if op == 'q':
            saved.append(dict(st)); st['ctm'] = list(st['ctm'])
        elif op == 'Q':
            if saved: st = saved.pop()
        elif op == 'cm':
            M = [float(v) for v in stack[-6:]]; st['ctm'] = mul(M, st['ctm'])
        elif op == 'RG':
            st['rgb'] = tuple(round(float(v), 4) for v in stack[-3:])
        elif op == 'G':
            v = round(float(stack[-1]), 4); st['rgb'] = (v, v, v)
        elif op == 'w':
            st['lw'] = float(stack[-1])
        elif op == 'd':
            st['dash'] = ' '.join(str(v) for v in stack)
        elif op == 'gs':
            st['gs'] = stack[-1]
        elif op == 'm':
            cur.append([app(st['ctm'], stack[-2], stack[-1])])
        elif op == 'l':
            if cur: cur[-1].append(app(st['ctm'], stack[-2], stack[-1]))
        elif op == 'c':
            if cur: cur[-1].append(app(st['ctm'], stack[-2], stack[-1])); curve = True
        elif op == 're':
            x, y, w, h = stack[-4:]
            cur.append([app(st['ctm'], x, y), app(st['ctm'], x + w, y), app(st['ctm'], x + w, y + h), app(st['ctm'], x, y + h)])
        elif op in ('S', 's', 'f', 'F', 'f*', 'B', 'B*', 'b', 'b*', 'n'):
            for sub in cur:
                paths.append(dict(pts=sub, rgb=st['rgb'], lw=st['lw'], dash=st['dash'], gs=st['gs'], op=op, curve=curve))
            cur, curve = [], False
        stack = []

print('streams', len(streams), 'painted subpaths', len(paths))
long = [p for p in paths if len(p['pts']) >= 50 and p['op'] in ('S', 's')]
print('long stroked paths', len(long))
for p in long:
    xs = [q[0] for q in p['pts']]; ys = [q[1] for q in p['pts']]
    print('  n=%d rgb=%s lw=%.2f dash=%r gs=%s x=[%.2f,%.2f] y=[%.2f,%.2f]' % (len(xs), p['rgb'], p['lw'], p['dash'], p['gs'], min(xs), max(xs), min(ys), max(ys)))

# grid lines: stroked two point axis aligned segments with the matplotlib grid colour #b0b0b0
grid = [p for p in paths if len(p['pts']) == 2 and p['op'] in ('S', 's') and abs(p['rgb'][0] - 0.6902) < 0.01]
vert = sorted(set((round(p['pts'][0][0], 3), round(min(p['pts'][0][1], p['pts'][1][1]), 2), round(max(p['pts'][0][1], p['pts'][1][1]), 2)) for p in grid if abs(p['pts'][0][0] - p['pts'][1][0]) < 1e-6))
hori = sorted(set((round(p['pts'][0][1], 3), round(min(p['pts'][0][0], p['pts'][1][0]), 2), round(max(p['pts'][0][0], p['pts'][1][0]), 2)) for p in grid if abs(p['pts'][0][1] - p['pts'][1][1]) < 1e-6))
print('vertical grid lines', len(vert)); print('horizontal grid lines', len(hori))
# panels: group the horizontal grid lines by their x extent
panels = defaultdict(lambda: dict(v=[], h=[]))
for y, x0, x1 in hori:
    panels[(x0, x1)]['h'].append(y)
for x, y0, y1 in vert:
    for (x0, x1), P in panels.items():
        if x0 - 1e-6 <= x <= x1 + 1e-6:
            P['v'].append(x)
keys = sorted(panels)
for k in keys:
    print('panel x extent', k, 'vertical', [round(v, 2) for v in panels[k]['v']], 'horizontal', [round(v, 2) for v in panels[k]['h']])

def fit(pix, val):
    n = len(pix); mx = sum(pix) / n; mv = sum(val) / n
    a = sum((p - mx) * (v - mv) for p, v in zip(pix, val)) / sum((p - mx) ** 2 for p in pix)
    b = mv - a * mx
    res = max(abs(a * p + b - v) for p, v in zip(pix, val))
    return a, b, res

xt = [285 + 2.5 * i for i in range(7)]
ticks_y = {1: [-1, 0, 1, 2, 3], 2: [-5, -4, -3, -2, -1, 0, 1, 2]}
out = {}
for idx in (1, 2):
    k = keys[idx]; P = panels[k]
    vx = sorted(P['v']); hy = sorted(P['h'])
    assert len(vx) == 7 and len(hy) == len(ticks_y[idx]), (len(vx), len(hy))
    ax_, bx_, rx = fit(vx, xt); ay_, by_, ry = fit(hy, ticks_y[idx])
    print('panel', idx, 'x map residual %.2e, y map residual %.2e' % (rx, ry))
    curves = []
    for p in long:
        xs = [q[0] for q in p['pts']]
        if not (k[0] - 1 <= min(xs) and max(xs) <= k[1] + 1):
            continue
        X = [ax_ * q[0] + bx_ for q in p['pts']]; Y = [ay_ * q[1] + by_ for q in p['pts']]
        i = min(range(len(Y)), key=lambda j: Y[j])
        dx = sorted(X[j + 1] - X[j] for j in range(len(X) - 1))
        c = dict(rgb=p['rgb'], n=len(X), x_first=X[0], x_last=X[-1], y_min=Y[i], x_at_min=X[i],
                 dx_min=dx[0], dx_median=dx[len(dx) // 2], dx_max=dx[-1],
                 neighbours=[(round(X[j], 4), round(Y[j], 4)) for j in range(max(0, i - 3), min(len(X), i + 4))],
                 vertices=[[X[j], Y[j]] for j in range(len(X))])
        curves.append(c)
        print('  curve rgb=%s n=%d x=[%.4f, %.4f] min %.4f at %.4f; dx min %.2e median %.2e max %.2e' % (p['rgb'], len(X), X[0], X[-1], Y[i], X[i], dx[0], dx[len(dx) // 2], dx[-1]))
        print('     points around the minimum', c['neighbours'])
    # dotted vertical line (axvline at the hypothetical ordinate): two point vertical, dashed, not grid colour
    dots = [p for p in paths if len(p['pts']) == 2 and p['op'] in ('S', 's') and p['dash'] and abs(p['pts'][0][0] - p['pts'][1][0]) < 1e-6 and k[0] <= p['pts'][0][0] <= k[1]]
    dl = [ax_ * p['pts'][0][0] + bx_ for p in dots]
    print('  dotted vertical lines at', [round(v, 4) for v in dl])
    out['panel%d' % idx] = dict(curves=curves, dotted=dl, xmap=[ax_, bx_], ymap=[ay_, by_])
json.dump(out, open('fig_extract.json', 'w'), indent=1)
