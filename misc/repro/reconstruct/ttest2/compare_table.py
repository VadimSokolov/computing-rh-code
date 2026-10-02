# Reconstructed comparison, not the authors' script. It compares Table tab:tilt of ch/ch04.tex, entry by
# entry, with the output of the archived code/ttest2.py run with the recovered command lines (the log of
# run_archived.slurm, given as the first argument) and with floor_check.json. Errors are quoted in the table
# to two significant digits, so an error entry agrees when the printed error rounds to it; "<b" agrees when
# the error is below b. Nodes and the last column (0.341 theta digits) are exact; times are compared as a ratio.
import sys, re, json, math

TABLE = [  # delta, bits, nodes, theta, error, ms, untilted digits  (ch/ch04.tex, Table tab:tilt)
    (0.10, 200, 2001, 200, '1.6e-45', 7, 68), (0.10, 200, 2001, 500, '9.8e-42', 7, 171),
    (0.10, 200, 2001, 1000, '1.8e-20', 7, 341), (0.02, 220, 8401, 1000, '<5e-45', 28, 341),
    (0.02, 220, 8401, 2000, '<4e-45', 28, 682), (0.02, 220, 8401, 5000, '2.9e-25', 28, 1705)]

log = open(sys.argv[1]).read()
blocks = re.split(r'== ttest2.py ', log)[1:3]
run = {}
for b in blocks:
    d = float(b.split()[0]); nodes = int(re.search(r'nodes (\d+)', b).group(1))
    for m in re.finditer(r'^(\d+) err (\S+) relrad (\S+) sec (\S+)', b, re.M):
        run[(d, int(m.group(1)))] = dict(nodes=nodes, err=float(m.group(2)), sec=float(m.group(4)))
fc = {(r['delta'], r['theta']): r for r in json.load(open('floor_check.json'))['rows']}

rows = []; n_ok = 0; n = 0
print('delta bits theta | nodes book/run | error book / run (3 digits) / full / rel diff | ms book / run median | digits book / 0.341094 theta')
for d, bits, nodes, th, err, ms, dig in TABLE:
    r = run[(d, th)]; f = fc[(d, th)]
    full = float(f['err'])
    if err.startswith('<'):
        ok_err = full < float(err[1:]); rel = None
    else:
        ok_err = float('%.1e' % full) == float(err); rel = abs(full - float(err)) / float(err)
    ok_nodes = r['nodes'] == nodes
    dig_run = round(math.pi / (4 * math.log(10)) * th)
    ok_dig = dig_run == dig
    ms_run = f['ms_median']
    n += 3; n_ok += ok_err + ok_nodes + ok_dig
    rows.append(dict(delta=d, bits=bits, theta=th, nodes_book=nodes, nodes_run=r['nodes'], nodes_ok=ok_nodes,
                     err_book=err, err_run_printed=r['err'], err_run_full=f['err'], err_rel_diff=rel, err_ok=ok_err,
                     err_re_run=f['err_re'], err_true=f['err_true'], ms_book=ms, ms_run_single=round(1000*r['sec'], 1),
                     ms_run_median=ms_run, ms_ratio=round(ms_run/ms, 2), digits_book=dig, digits_formula=dig_run, digits_ok=ok_dig))
    print(f"{d:5} {bits} {th:5} | {nodes}/{r['nodes']} | {err:>8} / {r['err']:.3g} / {f['err']} / "
          f"{'-' if rel is None else '%.3f' % rel} {'ok' if ok_err else 'DIFFERS'} | {ms} / {ms_run} | {dig} / {dig_run}")
json.dump(dict(rows=rows, checked=n, agree=n_ok), open('compare.json', 'w'), indent=1)
print(f'{n_ok} of {n} computed entries agree (nodes, error, last column); times are machine dependent')
