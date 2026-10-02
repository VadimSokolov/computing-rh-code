# Reconstructed collector, not the authors' script. It reads the outputs of search_task.py (the archived
# code/ttest2.py run for each candidate step h) and compares them with Table tab:tilt of ch/ch04.tex.
# An error entry matches when the printed error rounds to the two significant digits of the table; an
# entry "<b" matches when the printed error is below b. Writes search_results.json and prints a summary.
import json, glob, re
from decimal import Decimal as D

TABLE = {'A': {'nodes': 2001, 200: ('=', '1.6e-45'), 500: ('=', '9.8e-42'), 1000: ('=', '1.8e-20')},
         'B': {'nodes': 8401, 1000: ('<', '5e-45'), 2000: ('<', '4e-45'), 5000: ('=', '2.9e-25')}}

def two_digits(x):
    return float('%.1e' % x)

def parse(stdout):
    nodes = None; rows = {}
    for line in stdout.splitlines():
        m = re.match(r'nodes (\d+)', line)
        if m: nodes = int(m.group(1)); continue
        m = re.match(r'(\d+) err (\S+) relrad (\S+) sec (\S+)', line)
        if m: rows[int(m.group(1))] = dict(err=float(m.group(2)), rad=float(m.group(3)), sec=float(m.group(4)))
    return nodes, rows

res = []; raw = []
for f in sorted(glob.glob('out/search_*.json'), key=lambda s: int(re.findall(r'\d+', s)[-1])):
    r = json.load(open(f)); raw.append(r)
    nodes, rows = parse(r['stdout'])
    tab = TABLE[r['config']]
    ok = {}
    for th, (op, v) in ((k, v) for k, v in tab.items() if k != 'nodes'):
        e = rows.get(th, {}).get('err')
        if e is None: ok[th] = False
        elif op == '=': ok[th] = two_digits(e) == float(v)
        else: ok[th] = e < float(v)
    res.append(dict(task=r['task'], config=r['config'], h=r['h'], xmax=r['xmax'], rc=r['rc'], nodes=nodes,
                    nodes_ok=(nodes == tab['nodes']), rows=rows, match=ok, all_match=all(ok.values()) and nodes == tab['nodes']))
json.dump(res, open('search_results.json', 'w'), indent=1)
json.dump(raw, open('search_raw.json', 'w'), indent=1)   # the stdout of every archived ttest2.py run

for c in 'AB':
    R = [x for x in res if x['config'] == c]
    ths = [k for k in TABLE[c] if k != 'nodes']
    print(f'== config {c}: {len(R)} runs, table', {k: TABLE[c][k] for k in ths})
    print('   h          xmax       nodes  ' + '  '.join(f'err({t})' for t in ths) + '   matches')
    for x in sorted(R, key=lambda x: (D(x['h']), D(x['xmax']))):
        errs = '  '.join('%9.3g' % x['rows'].get(t, {}).get('err', float('nan')) for t in ths)
        flags = ''.join('Y' if x['match'][t] else '.' for t in ths)
        print(f"   {x['h']:<10} {x['xmax']:<10} {x['nodes']}  {errs}   {flags}{'  ALL' if x['all_match'] else ''}")
    for t in ths:
        print(f'   entry theta={t}: matched by h =', [x['h'] for x in R if x['match'][t]])
    print('   all entries matched by (h, xmax) =', [(x['h'], x['xmax']) for x in R if x['all_match']])
    bad = [x['task'] for x in R if x['rc'] != 0 or not x['nodes_ok']]
    print('   runs with rc != 0 or wrong node count:', bad)
