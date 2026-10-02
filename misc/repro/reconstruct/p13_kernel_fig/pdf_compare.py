# pdf_compare.py (round 2 item kernel_fig, not part of the authors' package). Compares the archived fig/kernel.pdf
# (old_kernel.pdf) with the two PDFs drawn on Hopper by kernel_eta.py: kernel_oldlabels.pdf (archived labels, a control
# for the change of matplotlib version) and kernel.pdf (labels of notation decision 15). It reads the vector content of
# each page (stroked and filled paths with colour and line width, glyph strings of every text block, placed XObjects),
# compares the three drawings element by element, maps the curves of the archived PDF back to data coordinates and
# checks them against the values computed by kernel_eta.py, and makes a pixel difference of the three PNG renders.
# Usage: python3 pdf_compare.py   (in a folder holding the three PDFs, the two *_curves.json files and the three
# pdftoppm renders at 150 dpi (pdftoppm -png -r 150 -singlefile) old_kernel.png, kernel_oldlabels.png, kernel.png); writes pdf_compare.json and diff_*.png.
import re, zlib, json, sys
import numpy as np

def objects(data):
    out = {}
    for m in re.finditer(rb'(\d+) 0 obj\s*(.*?)endobj', data, re.S):
        body = m.group(2)
        k = body.find(b'stream')
        if k >= 0 and body[:k].rstrip().endswith(b'>>'):
            head = body[:k]
            raw = body[k + 6:]
            raw = raw[1:] if raw[:1] == b'\n' else raw[2:] if raw[:2] == b'\r\n' else raw
            raw = raw[:raw.rfind(b'endstream')].rstrip(b'\r\n')
            if b'FlateDecode' in head:
                raw = zlib.decompressobj().decompress(raw)
            out[int(m.group(1))] = (head, raw)
        else:
            out[int(m.group(1))] = (body, None)
    return out

def tokens(s):
    i, n = 0, len(s)
    ws, delim = b' \t\r\n\x00\x0c', b'()<>[]{}/%'
    while i < n:
        c = s[i:i + 1]
        if c in (b' ', b'\t', b'\r', b'\n', b'\x00', b'\x0c'):
            i += 1; continue
        if c == b'(':
            depth, j, buf = 1, i + 1, bytearray()
            while j < n:
                ch = s[j:j + 1]
                if ch == b'\\':
                    nx = s[j + 1:j + 2]
                    if nx.isdigit():
                        k = j + 1
                        while k < j + 4 and s[k:k + 1].isdigit(): k += 1
                        buf.append(int(s[j + 1:k], 8)); j = k; continue
                    buf += {b'n': b'\n', b'r': b'\r', b't': b'\t', b'b': b'\b', b'f': b'\f'}.get(nx, nx); j += 2; continue
                if ch == b'(': depth += 1
                elif ch == b')':
                    depth -= 1
                    if depth == 0: j += 1; break
                buf += ch; j += 1
            yield ('str', bytes(buf)); i = j; continue
        if s[i:i + 2] in (b'<<', b'>>'):
            yield ('dict', s[i:i + 2]); i += 2; continue
        if c == b'<':
            j = s.index(b'>', i); yield ('str', bytes.fromhex(s[i + 1:j].decode())); i = j + 1; continue
        if c in (b'[', b']'):
            yield ('arr', c); i += 1; continue
        if c == b'/':
            j = i + 1
            while j < n and s[j:j + 1] not in ws and s[j:j + 1] not in delim: j += 1
            yield ('name', s[i + 1:j].decode('latin-1')); i = j; continue
        j = i
        while j < n and s[j:j + 1] not in ws and s[j:j + 1] not in delim: j += 1
        w = s[i:j]
        try: yield ('num', float(w))
        except ValueError: yield ('op', w.decode('latin-1'))
        i = j

def mul(a, b):
    return [a[0] * b[0] + a[1] * b[2], a[0] * b[1] + a[1] * b[3], a[2] * b[0] + a[3] * b[2], a[2] * b[1] + a[3] * b[3],
            a[4] * b[0] + a[5] * b[2] + b[4], a[4] * b[1] + a[5] * b[3] + b[5]]

def app(m, x, y):
    return (m[0] * x + m[2] * y + m[4], m[1] * x + m[3] * y + m[5])

def hexcol(c):
    return '#%02x%02x%02x' % tuple(int(round(255 * v)) for v in c)

def resources(objs):
    """font resource name -> {code: glyph name}; ExtGState name -> (CA, ca)."""
    fonts, alphas = {}, {}
    alltext = b''.join(h for h, r in objs.values())
    for m in re.finditer(rb'/(F\d+)\s+(\d+)\s+0\s+R', alltext):
        head = objs.get(int(m.group(2)), (b'', None))[0]
        enc = {}
        dm = re.search(rb'/Differences\s*\[(.*?)\]', head, re.S)
        if dm:
            code = 0
            for t in re.findall(rb'/[^\s/\[\]]+|\d+', dm.group(1)):
                if t.startswith(b'/'): enc[code] = t[1:].decode('latin-1'); code += 1
                else: code = int(t)
        bf = re.search(rb'/BaseFont\s*/([^\s/]+)', head)
        fonts[m.group(1).decode()] = dict(base=bf.group(1).decode() if bf else '?', enc=enc)
    for m in re.finditer(rb'/(A\d+)\s*<<(.*?)>>', alltext, re.S):
        ca = re.search(rb'/CA\s+([\d.]+)', m.group(2)); cf = re.search(rb'/ca\s+([\d.]+)', m.group(2))
        alphas[m.group(1).decode()] = (float(ca.group(1)) if ca else 1.0, float(cf.group(1)) if cf else 1.0)
    return fonts, alphas

def page(fn):
    data = open(fn, 'rb').read()
    objs = objects(data)
    pg = [h for h, r in objs.values() if re.search(rb'/Type\s*/Page\b', h)][0]
    content = objs[int(re.search(rb'/Contents\s+(\d+)\s+0\s+R', pg).group(1))][1]
    fonts, alphas = resources(objs)
    st = dict(ctm=[1, 0, 0, 1, 0, 0], RG=(0, 0, 0), rg=(0, 0, 0), lw=1.0, CA=1.0, ca=1.0, dash='')
    stack, ops, path, arr = [], [], [], None
    strokes, fills, dos, texts = [], [], [], []
    tb = None
    for kind, v in tokens(content):
        if kind == 'arr':
            if v == b'[': arr = []
            else: ops.append(arr); arr = None
            continue
        if arr is not None: arr.append(v); continue
        if kind != 'op': ops.append(v); continue
        op, o, ops = v, ops, []
        if op == 'q': stack.append(dict(st))
        elif op == 'Q': st = stack.pop()
        elif op == 'cm': st['ctm'] = mul(o[-6:], st['ctm'])
        elif op == 'w': st['lw'] = o[-1]
        elif op == 'd': st['dash'] = str(o)
        elif op == 'RG': st['RG'] = tuple(o[-3:])
        elif op == 'rg': st['rg'] = tuple(o[-3:])
        elif op == 'G': st['RG'] = (o[-1],) * 3
        elif op == 'g': st['rg'] = (o[-1],) * 3
        elif op == 'gs':
            st['CA'], st['ca'] = alphas.get(o[-1], (1.0, 1.0))
        elif op == 'm': path.append([app(st['ctm'], o[-2], o[-1])])
        elif op == 'l': path[-1].append(app(st['ctm'], o[-2], o[-1]))
        elif op == 'c': path[-1].append(app(st['ctm'], o[-2], o[-1]))
        elif op == 'h': pass
        elif op == 're':
            x, y, w, h = o[-4:]
            path.append([app(st['ctm'], x, y), app(st['ctm'], x + w, y), app(st['ctm'], x + w, y + h), app(st['ctm'], x, y + h)])
        elif op in ('S', 's', 'f', 'f*', 'F', 'B', 'B*', 'b', 'b*', 'n', 'W', 'W*'):
            if op in ('W', 'W*'): continue
            if op in ('S', 's', 'B', 'B*', 'b', 'b*'):
                for sub in path: strokes.append(dict(col=hexcol(st['RG']), lw=round(st['lw'], 4), a=st['CA'], dash=st['dash'], pts=sub))
            if op in ('f', 'f*', 'F', 'B', 'B*', 'b', 'b*'):
                for sub in path: fills.append(dict(col=hexcol(st['rg']), a=st['ca'], pts=sub))
            path = []
        elif op == 'Do':
            dos.append(dict(name=o[-1], at=app(st['ctm'], 0, 0), scale=st['ctm'][0]))
        elif op == 'BT':
            tb = dict(font=None, size=None, glyphs=[], at=None, tm=[1, 0, 0, 1, 0, 0], tl=[1, 0, 0, 1, 0, 0])
        elif op == 'ET':
            texts.append(tb); tb = None
        elif op == 'Tf': tb['font'], tb['size'] = o[-2], o[-1]
        elif op == 'Td':
            tb['tl'] = mul([1, 0, 0, 1, o[-2], o[-1]], tb['tl']); tb['tm'] = list(tb['tl'])
        elif op == 'Tm':
            tb['tl'] = list(o[-6:]); tb['tm'] = list(o[-6:])
        elif op in ('Tj', 'TJ'):
            parts = o[-1] if op == 'TJ' else [o[-1]]
            enc = fonts.get(tb['font'], {}).get('enc', {})
            pos = app(mul(tb['tm'], st['ctm']), 0, 0)
            if tb['at'] is None: tb['at'] = pos
            for p in parts:
                if isinstance(p, bytes):
                    tb['glyphs'] += [enc.get(b, 'code%d' % b) for b in p]
    return dict(strokes=strokes, fills=fills, dos=dos, texts=texts, fonts={k: v['base'] for k, v in fonts.items()})

def pdiff(a, b):
    """max abs difference of two point lists of the same length."""
    return float(np.max(np.abs(np.array(a) - np.array(b)))) if len(a) == len(b) and len(a) else None

def compare(A, B, tag):
    r = dict(pair=tag)
    # strokes: same number, same colours and widths in the same order, and the largest vertex difference
    sa, sb = A['strokes'], B['strokes']
    r['n_strokes'] = [len(sa), len(sb)]
    same_style = len(sa) == len(sb) and all((p['col'], p['lw'], len(p['pts'])) == (q['col'], q['lw'], len(q['pts'])) for p, q in zip(sa, sb))
    r['strokes_same_style_and_length'] = same_style
    if same_style:
        d = [pdiff(p['pts'], q['pts']) for p, q in zip(sa, sb)]
        r['strokes_max_vertex_diff_pt'] = max(d)
        bad = [(i, sa[i]['col'], sa[i]['lw'], len(sa[i]['pts']), d[i]) for i in range(len(d)) if d[i] > 1e-3]
        r['strokes_differing'] = bad
    else:
        r['strokes_a'] = sorted({(p['col'], p['lw']) for p in sa}); r['strokes_b'] = sorted({(p['col'], p['lw']) for p in sb})
        # match by colour and width and compare those with equal length
        info = []
        for key in sorted({(p['col'], p['lw']) for p in sa} | {(p['col'], p['lw']) for p in sb}):
            pa = [p for p in sa if (p['col'], p['lw']) == key]; pb = [p for p in sb if (p['col'], p['lw']) == key]
            info.append(dict(key=key, na=len(pa), nb=len(pb), lens_a=[len(p['pts']) for p in pa], lens_b=[len(p['pts']) for p in pb],
                             diffs=[pdiff(p['pts'], q['pts']) for p, q in zip(pa, pb)]))
        r['strokes_by_style'] = info
    fa, fb = A['fills'], B['fills']
    r['n_fills'] = [len(fa), len(fb)]
    r['fills_max_vertex_diff_pt'] = max([pdiff(p['pts'], q['pts']) or 0 for p, q in zip(fa, fb)] or [0]) if len(fa) == len(fb) and all(len(p['pts']) == len(q['pts']) for p, q in zip(fa, fb)) else 'differ'
    if r['fills_max_vertex_diff_pt'] == 'differ':
        r['fills_a'] = [(p['col'], p['a'], len(p['pts'])) for p in fa]; r['fills_b'] = [(p['col'], p['a'], len(p['pts'])) for p in fb]
    # placed XObjects (tick marks, and glyphs outside the 8 bit range such as Greek letters of mathtext)
    da = [(x['name'], round(x['at'][0], 3), round(x['at'][1], 3)) for x in A['dos']]
    db = [(x['name'], round(x['at'][0], 3), round(x['at'][1], 3)) for x in B['dos']]
    r['n_xobjects'] = [len(da), len(db)]
    r['xobjects_only_a'] = sorted(set(da) - set(db)); r['xobjects_only_b'] = sorted(set(db) - set(da))
    # text blocks: glyph strings and positions
    ta = [(' '.join(t['glyphs']), round(t['size'] or 0, 2), round(t['at'][0], 2), round(t['at'][1], 2)) for t in A['texts'] if t['at']]
    tb = [(' '.join(t['glyphs']), round(t['size'] or 0, 2), round(t['at'][0], 2), round(t['at'][1], 2)) for t in B['texts'] if t['at']]
    r['n_texts'] = [len(ta), len(tb)]
    r['texts_only_a'] = sorted(set(ta) - set(tb)); r['texts_only_b'] = sorted(set(tb) - set(ta))
    return r

def boxes(P):
    """axes boxes from the black spines (lw 0.8, longer than 50 pt; the tick marks are 3.5 pt strokes of the same style):
    left spines are vertical, bottom spines horizontal."""
    sp = [s for s in P['strokes'] if s['col'] == '#000000' and abs(s['lw'] - 0.8) < 1e-6 and len(s['pts']) == 2
          and np.hypot(s['pts'][1][0] - s['pts'][0][0], s['pts'][1][1] - s['pts'][0][1]) > 50]
    vert = sorted([s['pts'] for s in sp if abs(s['pts'][0][0] - s['pts'][1][0]) < 1e-9], key=lambda p: p[0][0])
    hori = sorted([s['pts'] for s in sp if abs(s['pts'][0][1] - s['pts'][1][1]) < 1e-9], key=lambda p: min(p[0][0], p[1][0]))
    out = []
    for v, h in zip(vert, hori):
        out.append([v[0][0], min(v[0][1], v[1][1]), max(h[0][0], h[1][0]), max(v[0][1], v[1][1])])
    return out

def styles(P):
    """summary of the strokes: (colour, width) -> number of subpaths and their vertex counts."""
    out = {}
    for s in P['strokes']:
        k = '%s lw=%s' % (s['col'], s['lw'])
        out.setdefault(k, []).append(len(s['pts']))
    return {k: dict(n=len(v), vertices=v if len(v) <= 12 else '%d paths, %d vertices in all' % (len(v), sum(v))) for k, v in out.items()}

def to_data(P, J):
    """curves of a PDF in data coordinates, checked against the values of a *_curves.json."""
    bx = boxes(P)
    res = dict(boxes_pt=bx)
    ax0, ax1 = J['axes']
    u, Phi, x = np.array(J['u']), np.array(J['Phi']), np.array(J['x'])
    def mapx(px, b, lim): return lim[0] + (px - b[0]) / (b[2] - b[0]) * (lim[1] - lim[0])
    # left panel: the one stroke of width 2
    if len(bx) < 2:
        return res
    left = [s for s in P['strokes'] if s['lw'] == 2.0 and len(s['pts']) > 10]
    b = bx[0]
    pts = np.array([p for s in left for p in s['pts']])
    res['left_subpaths'] = len(left)
    X = mapx(pts[:, 0], b, ax0['xlim']); Y = ax0['ylim'][0] + (pts[:, 1] - b[1]) / (b[3] - b[1]) * (ax0['ylim'][1] - ax0['ylim'][0])
    res['left'] = dict(n=len(pts), max_abs_err=float(np.max(np.abs(Y - np.interp(X, u, Phi)))), peak=float(Y.max()))
    # right panel: the strokes of width 1.4 with more than 2 vertices (legend handles have 2), one colour per curve
    colours = ['#1b4f72', '#c0392b', '#27ae60', '#d68910']
    b = bx[1]; ly = np.log10(ax1['ylim'])
    res['right'] = []
    for col, (d, c) in zip(colours, J['curves'].items()):
        sub = [s['pts'] for s in P['strokes'] if s['lw'] == 1.4 and s['col'] == col and len(s['pts']) > 2]
        if not sub:
            res['right'].append(dict(d=d, colour=col, n=0)); continue
        s = dict(col=col, pts=[p for q in sub for p in q])
        pts = np.array(s['pts'])
        X = mapx(pts[:, 0], b, ax1['xlim']); LY = ly[0] + (pts[:, 1] - b[1]) / (b[3] - b[1]) * (ly[1] - ly[0])
        inside = (pts[:, 1] >= b[1] - 1e-6) & (pts[:, 1] <= b[3] + 1e-6)
        ref = np.log10(np.maximum(np.array(c['v']), 1e-30))
        # a simplified path keeps a subset of the 601 vertices: compare each kept vertex with the nearest grid point
        k = np.clip(np.rint((X - x[0]) / (x[1] - x[0])).astype(int), 0, len(x) - 1)
        on_grid = np.abs(X - x[k]) < 1e-6
        m = inside & on_grid
        res['right'].append(dict(d=d, colour=s['col'], subpaths=len(sub), n=len(pts), n_compared=int(m.sum()), n_off_grid_inside=int((inside & ~on_grid).sum()),
                                 max_abs_err_log10=float(np.max(np.abs(LY[m] - ref[k[m]]))) if m.any() else None,
                                 peak_log10=float(LY[inside].max()) if inside.any() else None))
    return res

def pixels(fa, fb, out):
    from PIL import Image
    from scipy import ndimage
    A = np.asarray(Image.open(fa).convert('RGB')).astype(int); B = np.asarray(Image.open(fb).convert('RGB')).astype(int)
    if A.shape != B.shape: return dict(shape=[A.shape, B.shape])
    diff = np.abs(A - B).max(axis=2)
    mask = diff > 24
    lab, n = ndimage.label(ndimage.binary_dilation(mask, iterations=6))
    regions = []
    for sl in ndimage.find_objects(lab):
        ys, xs = sl
        regions.append(dict(x=[xs.start, xs.stop], y=[ys.start, ys.stop], pixels=int(mask[sl].sum())))
    img = (0.3 * B.mean(axis=2) + 0.7 * 255).astype(np.uint8)
    rgb = np.stack([img] * 3, axis=2)
    rgb[mask] = [220, 0, 0]
    Image.fromarray(rgb).save(out)
    return dict(shape=list(A.shape), n_pixels_differ=int(mask.sum()), n_pixels_any_diff=int((diff > 0).sum()), max_channel_diff=int(diff.max()), regions=regions)

if __name__ == '__main__':
    P = {f: page(f) for f in ['old_kernel.pdf', 'kernel_oldlabels.pdf', 'kernel.pdf']}
    J = {k: json.load(open(k + '_curves.json')) for k in ['kernel_oldlabels', 'kernel']}
    R = dict(fonts={f: P[f]['fonts'] for f in P}, strokes={f: styles(P[f]) for f in P}, fills={f: [(p['col'], p['a'], len(p['pts'])) for p in P[f]['fills']] for f in P})
    R['archived_vs_control'] = compare(P['old_kernel.pdf'], P['kernel_oldlabels.pdf'], 'old_kernel.pdf vs kernel_oldlabels.pdf')
    R['control_vs_new'] = compare(P['kernel_oldlabels.pdf'], P['kernel.pdf'], 'kernel_oldlabels.pdf vs kernel.pdf')
    R['archived_vs_new'] = compare(P['old_kernel.pdf'], P['kernel.pdf'], 'old_kernel.pdf vs kernel.pdf')
    R['data_check'] = {f: to_data(P[f], J['kernel']) for f in P}
    R['control_json_equals_new_json'] = (J['kernel']['Phi'] == J['kernel_oldlabels']['Phi'] and all(J['kernel']['curves'][d]['v'] == J['kernel_oldlabels']['curves'][d]['v'] for d in J['kernel']['curves']))
    R['axes_geometry'] = {k: J[k]['axes'] for k in J}
    R['matplotlib'] = {k: J[k]['matplotlib'] for k in J}
    import os
    if os.path.exists('kernel_curves_mpl3112.json'):
        # the same computation drawn by matplotlib 3.11.2 (shared venv): the plotted values must be identical
        K = json.load(open('kernel_curves_mpl3112.json'))
        R['values_equal_across_matplotlib_versions'] = (K['Phi'] == J['kernel']['Phi'] and all(K['curves'][d]['v'] == J['kernel']['curves'][d]['v'] for d in K['curves']))
    try:
        R['pixels_archived_vs_new'] = pixels('old_kernel.png', 'kernel.png', 'diff_old_new.png')
        R['pixels_archived_vs_control'] = pixels('old_kernel.png', 'kernel_oldlabels.png', 'diff_old_control.png')
        R['pixels_control_vs_new'] = pixels('kernel_oldlabels.png', 'kernel.png', 'diff_control_new.png')
    except FileNotFoundError as e:
        R['pixels'] = 'missing %s' % e
    json.dump(R, open('pdf_compare.json', 'w'), indent=1, default=str)
    print(json.dumps(R, indent=1, default=str))
