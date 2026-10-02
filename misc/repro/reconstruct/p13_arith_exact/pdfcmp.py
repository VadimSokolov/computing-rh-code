# New check script (item arith_exact): compares the drawing of two matplotlib PDFs. It inflates every FlateDecode
# stream (page content, fonts, glyph procedures) and reports whether the streams agree, then whether the files
# agree byte for byte. Used to compare the redrawn heattrace.pdf with the archived fig/heattrace.pdf (old_heattrace.pdf).
# Usage: python3 pdfcmp.py a.pdf b.pdf
import sys, re, zlib

def streams(path):
    data = open(path, 'rb').read()
    out = []
    for m in re.finditer(rb'stream\r?\n(.*?)\r?\nendstream', data, re.S):
        s = m.group(1)
        try:
            out.append(zlib.decompress(s))
        except zlib.error:
            out.append(s)
    return data, out

a, sa = streams(sys.argv[1]); b, sb = streams(sys.argv[2])
same = [x == y for x, y in zip(sa, sb)]
print('%s: %d streams, %d bytes; %s: %d streams, %d bytes' % (sys.argv[1], len(sa), len(a), sys.argv[2], len(sb), len(b)))
print('streams identical: %s (%d of %d pairs equal)' % (len(sa) == len(sb) and all(same), sum(same), min(len(sa), len(sb))))
for k, (x, y) in enumerate(zip(sa, sb)):
    if x != y:
        i = next(j for j in range(min(len(x), len(y))) if x[j] != y[j]) if any(p != q for p, q in zip(x, y)) else min(len(x), len(y))
        print('  stream %d differs (lengths %d, %d) from byte %d: %r | %r' % (k, len(x), len(y), i, x[max(0, i - 40):i + 40], y[max(0, i - 40):i + 40]))
print('files byte identical:', a == b)
