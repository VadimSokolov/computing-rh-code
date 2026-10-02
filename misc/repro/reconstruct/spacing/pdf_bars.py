# Helper of the reconstruction in misc/repro/reconstruct/spacing/ (not the authors' code): read the bars and the
# curves of matplotlib PDF figures (the book's fig/spacing.pdf, and the candidate drawn by spacing.py) from their
# content streams, to recover the binning of the histogram of Figure fig:ch15:spacing. Prints the rectangles and
# the end points of long paths in PDF units.
# Run: python3 pdf_bars.py FIGURE.pdf [FIGURE2.pdf ...]   (writes pdf_bars.json)
import json, re, sys, zlib

FILL = ('f', 'F', 'f*', 'S', 's', 'B', 'b', 'n')


def parse(fname):
    raw = open(fname, 'rb').read()
    streams = []
    for m in re.finditer(rb'<<(.*?)>>\s*stream\r?\n', raw, re.S):
        start = m.end(); end = raw.find(b'endstream', start)
        data = raw[start:end]
        if b'FlateDecode' in m.group(1):
            try: data = zlib.decompress(data)
            except zlib.error: continue
        streams.append(data.decode('latin-1'))
    rects, paths = [], []
    for s in streams:
        toks = s.split(); cur = []
        for i, t in enumerate(toks):
            try:
                if t in ('m', 'l') and i >= 2:
                    if t == 'm': cur = []
                    cur.append((float(toks[i - 2]), float(toks[i - 1])))
                elif t == 'c' and i >= 6:
                    cur.append((float(toks[i - 2]), float(toks[i - 1])))
                elif t == 're' and i >= 4:
                    x, y, w, h = map(float, toks[i - 4:i]); rects.append((x, y, x + w, y + h, 're'))
                elif t in FILL and cur:
                    xs = sorted(set(p[0] for p in cur)); ys = sorted(set(p[1] for p in cur))
                    if len(cur) in (4, 5) and len(xs) == 2 and len(ys) == 2:
                        rects.append((xs[0], ys[0], xs[1], ys[1], t))
                    elif len(cur) > 20:
                        paths.append(dict(op=t, n=len(cur), first=cur[0], last=cur[-1],
                                          ymin=min(p[1] for p in cur), ymax=max(p[1] for p in cur)))
                    cur = []
            except ValueError:
                cur = []
    big = max(streams, key=len)
    open(fname.replace('/', '_') + '.content.txt', 'w').write(big)
    alpha = sorted(set(re.findall(r'/c[aA] [\d.]+', raw.decode('latin-1'))))
    return dict(rects=rects, paths=paths, nstreams=len(streams), alpha=alpha, head=big[:1200])


if __name__ == '__main__':
    out = {}
    for f in sys.argv[1:]:
        r = out[f] = parse(f)
        print('==', f, ': streams', r['nstreams'], 'rectangles', len(r['rects']), 'long paths', len(r['paths']),
              'alpha settings', r['alpha'])
        for x in r['rects']: print('rect', x)
        for p in r['paths']: print('path', p)
        print('--- head of the largest stream ---'); print(r['head'])
    json.dump(out, open('pdf_bars.json', 'w'), indent=0)
