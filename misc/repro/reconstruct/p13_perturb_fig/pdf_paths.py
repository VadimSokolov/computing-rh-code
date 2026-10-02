# Helper of the reconstruction in misc/repro/reconstruct/p13_perturb_fig/ (not the authors' code): read the stroked
# paths of a matplotlib PDF figure (the book's fig/perturb.pdf, kept here as old_perturb.pdf, or a candidate) from its
# content stream. Writes NAME.content.txt (the decompressed page stream) and NAME.paths.json (every stroked path with
# all of its vertices in PDF points, the colour and line width in force, and the text of the tick labels).
# Run: python3 pdf_paths.py FIGURE.pdf [FIGURE2.pdf ...]
import json, re, sys, zlib


def streams(fname):
    raw = open(fname, 'rb').read()
    out = []
    for m in re.finditer(rb'<<(.*?)>>\s*stream\r?\n', raw, re.S):
        start = m.end(); end = raw.find(b'endstream', start)
        data = raw[start:end]
        if b'FlateDecode' in m.group(1):
            try: data = zlib.decompress(data)
            except zlib.error: continue
        out.append(data.decode('latin-1'))
    return out


def paths(content):
    toks = content.split()
    res, cur, lw, rgb, cm = [], [], None, None, None
    for i, t in enumerate(toks):
        try:
            if t == 'm':
                cur = [(float(toks[i - 2]), float(toks[i - 1]))]
            elif t == 'l':
                cur.append((float(toks[i - 2]), float(toks[i - 1])))
            elif t == 'c':
                cur.append((float(toks[i - 2]), float(toks[i - 1])))
            elif t == 'w':
                lw = float(toks[i - 1])
            elif t == 'RG':
                rgb = tuple(float(x) for x in toks[i - 3:i])
            elif t == 'cm':
                cm = tuple(float(x) for x in toks[i - 6:i])
            elif t in ('S', 's'):
                if cur: res.append(dict(n=len(cur), lw=lw, rgb=rgb, cm=cm, pts=cur))
                cur = []
            elif t in ('f', 'F', 'f*', 'B', 'b', 'n'):
                cur = []
        except ValueError:
            pass
    return res


if __name__ == '__main__':
    for f in sys.argv[1:]:
        ss = streams(f)
        big = max(ss, key=len)
        base = f.rsplit('/', 1)[-1].rsplit('.', 1)[0]
        open(base + '.content.txt', 'w').write(big)
        P = paths(big)
        json.dump(P, open(base + '.paths.json', 'w'))
        print('==', f, 'streams', len(ss), 'content length', len(big), 'stroked paths', len(P))
        for k, p in enumerate(P):
            print(k, 'n', p['n'], 'lw', p['lw'], 'rgb', p['rgb'], 'first', p['pts'][0], 'last', p['pts'][-1])
