# Comparison tool for the item ag_figures: prints the grid of the old density curves and the marker XObjects.
import json, re, zlib
import numpy as np
D = json.load(open('old_ag_almost.pdf.data.json'))
for ip in (1, 2):
    for c in D['panels'][ip]['curves']:
        x = np.array(c['x'])
        if len(x) < 100:
            continue
        for step in (0.05, 0.01, 0.005, 0.0025, 0.002, 0.001):
            k = (x - 284.0) / step
            print('panel', ip, c['color'], 'n', len(x), 'step %.4f: max distance to grid %.2e' % (step, np.max(np.abs(k - np.round(k))) * step))
        print('   first x', np.round(x[:12], 5).tolist())
        i = int(np.argmin(c['y'])); print('   min', x[i], c['y'][i], 'neighbours', np.round(x[i-3:i+4], 5).tolist(), np.round(c['y'][i-3:i+4], 4).tolist())
L = D['panels'][0]['curves']
for c in L:
    print('panel 0', c['color'], 'x', np.round(c['x'][:5], 6).tolist(), '...', np.round(c['x'][-3:], 6).tolist())
data = open('old_ag_cancel.pdf', 'rb').read()
for m in re.finditer(rb'(\d+) 0 obj\s*<<(.*?)>>\s*stream\r?\n', data, re.S):
    if b'/BBox' in m.group(2) and b'/Form' in m.group(2):
        start = m.end(); end = data.find(b'endstream', start); raw = data[start:end]
        if b'FlateDecode' in m.group(2): raw = zlib.decompressobj().decompress(raw)
        if len(raw) < 600 and b' c' in raw and b'Glyph' not in m.group(2):
            print('XObject', m.group(1), m.group(2)[:120], raw[:300])
