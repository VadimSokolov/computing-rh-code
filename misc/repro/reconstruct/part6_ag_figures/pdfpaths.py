# Comparison tool for the item ag_figures (not part of the authors' package): reads the vector paths of a matplotlib PDF
# (stroked and filled paths with colour, line width, dash and clip box, marker and glyph XObject positions, text font sizes)
# and maps the curves of every panel to data coordinates through its grid lines, so that the old figures
# fig/ag_almost.pdf and fig/ag_cancel.pdf can be compared numerically with the ones drawn by figs_almost.py.
# Usage: python3 pdfpaths.py file.pdf [file.pdf ...]; writes <file>.paths.json and prints a summary.
import sys, re, zlib, json


def pdf_streams(data):
    out = {}
    for m in re.finditer(rb'(\d+) 0 obj\s*<<(.*?)>>\s*stream\r?\n', data, re.S):
        start = m.end()
        end = data.find(b'endstream', start)
        raw = data[start:end]
        if raw.endswith(b'\r\n'):
            raw = raw[:-2]
        elif raw.endswith(b'\n'):
            raw = raw[:-1]
        if b'FlateDecode' in m.group(2):
            raw = zlib.decompressobj().decompress(raw)
        out[int(m.group(1))] = (m.group(2), raw)
    return out


def ext_alpha(data):
    """ExtGState names -> (stroke alpha, fill alpha)."""
    res = {}
    for m in re.finditer(rb'/(A\d+)\s*<<([^<>]*)>>', data):
        d = m.group(2)
        ca = re.search(rb'/CA\s+([\d.]+)', d)
        cf = re.search(rb'/ca\s+([\d.]+)', d)
        res[m.group(1).decode()] = (float(ca.group(1)) if ca else 1.0, float(cf.group(1)) if cf else 1.0)
    return res


def xobject_names(data):
    """XObject resource names -> object numbers (for glyph and marker XObjects)."""
    res = {}
    for m in re.finditer(rb'/([A-Za-z0-9_.\-]+)\s+(\d+)\s+0\s+R', data):
        res.setdefault(m.group(1).decode('latin-1'), int(m.group(2)))
    return res


def tokens(s):
    i, n = 0, len(s)
    ws = b' \t\r\n\x00\x0c'
    delim = b'()<>[]{}/%'
    while i < n:
        c = s[i:i + 1]
        if c in (b' ', b'\t', b'\r', b'\n', b'\x00', b'\x0c'):
            i += 1
            continue
        if c == b'%':
            while i < n and s[i:i + 1] not in (b'\r', b'\n'):
                i += 1
            continue
        if c == b'(':
            depth, j, buf = 1, i + 1, bytearray()
            while j < n and depth:
                ch = s[j:j + 1]
                if ch == b'\\':
                    buf += s[j:j + 2]
                    j += 2
                    continue
                if ch == b'(':
                    depth += 1
                elif ch == b')':
                    depth -= 1
                    if depth == 0:
                        j += 1
                        break
                buf += ch
                j += 1
            yield ('str', bytes(buf))
            i = j
            continue
        if s[i:i + 2] in (b'<<', b'>>'):
            yield ('dict', s[i:i + 2])
            i += 2
            continue
        if c == b'<':
            j = s.index(b'>', i)
            yield ('hex', s[i + 1:j])
            i = j + 1
            continue
        if c in (b'[', b']'):
            yield ('arr', c)
            i += 1
            continue
        if c == b'/':
            j = i + 1
            while j < n and s[j:j + 1] not in ws and s[j:j + 1] not in delim:
                j += 1
            yield ('name', s[i + 1:j].decode('latin-1'))
            i = j
            continue
        j = i
        while j < n and s[j:j + 1] not in ws and s[j:j + 1] not in delim:
            j += 1
        w = s[i:j]
        try:
            yield ('num', float(w))
        except ValueError:
            yield ('op', w.decode('latin-1'))
        i = j


def mul(a, b):
    """PDF matrices [a b c d e f]: returns a*b (apply a first, then b)."""
    return [a[0] * b[0] + a[1] * b[2], a[0] * b[1] + a[1] * b[3],
            a[2] * b[0] + a[3] * b[2], a[2] * b[1] + a[3] * b[3],
            a[4] * b[0] + a[5] * b[2] + b[4], a[4] * b[1] + a[5] * b[3] + b[5]]


def apply(m, x, y):
    return (m[0] * x + m[2] * y + m[4], m[1] * x + m[3] * y + m[5])


def parse(content, alphas, xobj):
    st = dict(ctm=[1, 0, 0, 1, 0, 0], RG=(0, 0, 0), rg=(0, 0, 0), lw=1.0, dash=([], 0), gs=None, clip=None)
    stack, operands, path, pending_clip = [], [], [], False
    strokes, fills, dos, texts = [], [], [], []
    font = None
    arr = None
    for kind, v in tokens(content):
        if kind == 'arr':
            if v == b'[':
                arr = []
            else:
                operands.append(arr)
                arr = None
            continue
        if arr is not None:
            arr.append(v)
            continue
        if kind != 'op':
            operands.append(v)
            continue
        op, o = v, operands
        operands = []
        if op == 'q':
            stack.append(dict(st))
        elif op == 'Q':
            st = stack.pop()
        elif op == 'cm':
            st['ctm'] = mul(o[-6:], st['ctm'])
        elif op == 'w':
            st['lw'] = o[-1]
        elif op == 'd':
            st['dash'] = (o[-2], o[-1])
        elif op == 'RG':
            st['RG'] = tuple(o[-3:])
        elif op == 'rg':
            st['rg'] = tuple(o[-3:])
        elif op == 'G':
            st['RG'] = (o[-1],) * 3
        elif op == 'g':
            st['rg'] = (o[-1],) * 3
        elif op == 'gs':
            st['gs'] = o[-1]
        elif op == 'm':
            path.append(['m', apply(st['ctm'], o[-2], o[-1])])
        elif op == 'l':
            path.append(['l', apply(st['ctm'], o[-2], o[-1])])
        elif op == 'c':
            path.append(['c', apply(st['ctm'], o[-2], o[-1])])
        elif op in ('v', 'y'):
            path.append(['c', apply(st['ctm'], o[-2], o[-1])])
        elif op == 'h':
            path.append(['h', None])
        elif op == 're':
            x, y, w, hh = o[-4:]
            for k, (px, py) in enumerate([(x, y), (x + w, y), (x + w, y + hh), (x, y + hh)]):
                path.append(['m' if k == 0 else 'l', apply(st['ctm'], px, py)])
            path.append(['h', None])
        elif op in ('W', 'W*'):
            pending_clip = True
        elif op in ('S', 's', 'f', 'F', 'f*', 'B', 'B*', 'b', 'b*', 'n'):
            if pending_clip:
                pts = [p for c, p in path if p is not None]
                xs, ys = [p[0] for p in pts], [p[1] for p in pts]
                st['clip'] = (min(xs), min(ys), max(xs), max(ys))
                pending_clip = False
            if op != 'n' and path:
                a = alphas.get(st['gs'], (1.0, 1.0)) if st['gs'] else (1.0, 1.0)
                rec = dict(verts=[[round(p[0], 6), round(p[1], 6)] for c, p in path if p is not None],
                           codes=''.join(c for c, p in path), clip=st['clip'])
                if op in ('S', 's', 'B', 'B*', 'b', 'b*'):
                    strokes.append(dict(rec, color=st['RG'], lw=st['lw'], dash=st['dash'], alpha=a[0]))
                if op in ('f', 'F', 'f*', 'B', 'B*', 'b', 'b*'):
                    fills.append(dict(rec, color=st['rg'], alpha=a[1]))
            path = []
        elif op == 'Do':
            x, y = apply(st['ctm'], 0, 0)
            dos.append(dict(name=o[-1], x=round(x, 6), y=round(y, 6), clip=st['clip'], color=st['RG'], fill=st['rg']))
        elif op == 'Tf':
            font = (o[-2], o[-1])
        elif op in ('Tj', 'TJ'):
            texts.append(dict(font=font, pos=apply(st['ctm'], 0, 0), s=repr(o[-1])[:60]))
    return dict(strokes=strokes, fills=fills, dos=dos, texts=texts)


def main(fn):
    data = open(fn, 'rb').read()
    ss = pdf_streams(data)
    alphas = ext_alpha(data)
    xobj = xobject_names(data)
    page = max((raw for d, raw in ss.values() if b'/Subtype' not in d), key=lambda r: r.count(b' l\n') + r.count(b' re'))
    res = parse(page, alphas, xobj)
    res['alphas'] = alphas
    json.dump(res, open(fn.rsplit('/', 1)[-1] + '.paths.json', 'w'))
    # summary: clip boxes (panels), then per panel the strokes by colour
    clips = sorted({tuple(round(c, 2) for c in s['clip']) for s in res['strokes'] if s['clip']}, key=lambda c: c[0])
    print('==', fn, 'strokes', len(res['strokes']), 'fills', len(res['fills']), 'Do', len(res['dos']), 'texts', len(res['texts']))
    print('ExtGState alphas', alphas)
    for c in clips:
        print('clip box', c)
        grp = {}
        for s in res['strokes']:
            if s['clip'] and tuple(round(q, 2) for q in s['clip']) == c:
                key = (tuple(round(q, 4) for q in s['color']), round(s['lw'], 3), str(s['dash']), s['alpha'])
                grp.setdefault(key, []).append(len(s['verts']))
        for k, v in grp.items():
            print('   stroke colour %s lw %s dash %s alpha %s: %d paths, vertices %s' % (k[0], k[1], k[2], k[3], len(v), v[:12]))
    fonts = {}
    for t in res['texts']:
        fonts[str(t['font'])] = fonts.get(str(t['font']), 0) + 1
    print('fonts', fonts)
    names = {}
    for d in res['dos']:
        names[d['name']] = names.get(d['name'], 0) + 1
    print('XObjects used', names)


if __name__ == '__main__':
    for fn in sys.argv[1:]:
        main(fn)
