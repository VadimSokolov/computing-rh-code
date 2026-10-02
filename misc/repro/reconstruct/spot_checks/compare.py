#!/usr/bin/env python3
# NOT an author script. Comparison helper for the reproducibility audit of
# verification_reciprocal_zeta/spot_checks_report.md.
#
# Inputs, all in the working folder:
#   archived_spot_checks_report.md  the archived report (verification_reciprocal_zeta/spot_checks_report.md)
#   spot_checks_report.md           the report regenerated on Hopper by make_report.py
#   run_all_spot_checks.out         the output of the same script in the full rerun (logs/spot_checks.out)
#   book_quotes.tsv                 the numbers the book quotes from the report, with the exact TeX they sit in
#   rz1.tex rz2.tex rz3.tex rz4.tex appA.tex ch12.tex   the chapters, as they were when the job was submitted
# Output: comparison.md and comparison.json.
#
# 1. Archived report against the regenerated one: lines aligned with difflib, then every
#    numeric token of the common lines compared (relative tolerance 1e-9).
# 2. Regenerated report against the run_all output: byte identity and token comparison.
# 3. Book quotations: each row of book_quotes.tsv is located in its chapter and its values
#    compared with the regenerated report at the precision the book prints.
# 4. Reference strings in spot_checks.py ("paper: ...") that no longer occur in the book.
import difflib
import hashlib
import json
import os
import re
import sys

from mpmath import mp, mpf, sqrt, pi, nstr

mp.dps = 50
HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
TOL = mpf('1e-9')
TOKEN = re.compile(r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?')
CHAPTERS = ['rz1.tex', 'rz2.tex', 'rz3.tex', 'rz4.tex', 'appA.tex', 'ch12.tex']


def read(name):
    with open(name, 'rb') as fh:
        return fh.read()


def ulp(s):
    """Unit in the last printed digit of a decimal string such as -0.587340 or 1.99e-13."""
    mant, _, ex = s.strip().lower().partition('e')
    d = len(mant.split('.')[1]) if '.' in mant else 0
    return mpf(10) ** ((int(ex) if ex else 0) - d)


def token_compare(la, lb):
    """Compare the numeric tokens of two lists of lines, position by position."""
    n, worst, bad = 0, mpf(0), []
    for x, y in zip(la, lb):
        ta, tb = TOKEN.findall(x), TOKEN.findall(y)
        if len(ta) != len(tb):
            bad.append((x, y, 'token count differs'))
            continue
        for a, b in zip(ta, tb):
            n += 1
            A, B = mpf(a), mpf(b)
            if A == B:
                continue
            r = abs(A - B) / max(abs(A), abs(B))
            worst = max(worst, r)
            if r > TOL:
                bad.append((x, y, f'{a} vs {b}, rel {nstr(r, 3)}'))
    return n, worst, bad


out, res = [], {}
arch_b, new_b, log_b = read('archived_spot_checks_report.md'), read('spot_checks_report.md'), read('run_all_spot_checks.out')
arch, new, log = arch_b.decode().splitlines(), new_b.decode().splitlines(), log_b.decode().splitlines()

# ---------------------------------------------------------------- 1. archive against rerun
sm = difflib.SequenceMatcher(a=arch, b=new, autojunk=False)
ops = sm.get_opcodes()
common_a, common_b, only_new, only_arch, changed = [], [], [], [], []
for tag, i1, i2, j1, j2 in ops:
    if tag == 'equal':
        common_a += arch[i1:i2]
        common_b += new[j1:j2]
    elif tag == 'insert':
        only_new.append((j1 + 1, j2, new[j1:j2]))
    elif tag == 'delete':
        only_arch.append((i1 + 1, i2, arch[i1:i2]))
    else:
        changed.append((i1 + 1, i2, j1 + 1, j2, arch[i1:i2], new[j1:j2]))
n_tok, worst, bad = token_compare(common_a, common_b)
chg_tok = 0
for c in changed:
    k, w, b = token_compare(c[4], c[5])
    chg_tok += k
    n_tok += k
    worst = max(worst, w)
    bad += b
res['archive_vs_rerun'] = {
    'archived_lines': len(arch), 'rerun_lines': len(new),
    'identical_lines': len(common_a), 'changed_blocks': len(changed),
    'lines_only_in_rerun': sum(i2 - i1 + 1 for i1, i2, _ in only_new),
    'lines_only_in_archive': sum(i2 - i1 + 1 for i1, i2, _ in only_arch),
    'numeric_tokens_compared': n_tok, 'max_rel_diff': nstr(worst, 3), 'beyond_tol': len(bad),
    'archive_is_prefix_plus_summary': arch == new[:len(arch) - 3] + new[-3:],
}
out.append('# spot_checks_report.md: archive, rerun and book\n')
out.append('## 1. Archived report against the regenerated report\n')
a = res['archive_vs_rerun']
out.append(f"Archived report: {a['archived_lines']} lines, {len(arch_b)} bytes. Regenerated report: {a['rerun_lines']} lines, {len(new_b)} bytes.\n")
out.append(f"Identical lines: {a['identical_lines']}. Changed blocks: {a['changed_blocks']}. Lines only in the regenerated report: {a['lines_only_in_rerun']}. Lines only in the archive: {a['lines_only_in_archive']}.\n")
out.append(f"Numeric tokens compared on the common lines (labels included): {n_tok}; beyond relative tolerance 1e-9: {len(bad)}; largest relative difference: {a['max_rel_diff']}.\n")
out.append(f"The archive equals the first {len(arch) - 3} lines of the regenerated report followed by its last 3 lines (blank line, summary header, 'checks failed: 0'): {a['archive_is_prefix_plus_summary']}.\n")
for j1, j2, block in only_new:
    heads = [s for s in block if s.startswith('--- ')]
    out.append(f"Lines {j1} to {j2} of the regenerated report are not in the archive ({j2 - j1 + 1} lines; sections: {', '.join(h.strip('- ') for h in heads)}).\n")
for b in bad[:20]:
    out.append(f"    DIFF: {b}\n")

# ---------------------------------------------------------------- 2. rerun against run_all output
k, w, b2 = token_compare(new, log)
res['rerun_vs_run_all'] = {'byte_identical': new_b == log_b, 'lines': [len(new), len(log)],
                           'numeric_tokens_compared': k, 'max_rel_diff': nstr(w, 3), 'beyond_tol': len(b2)}
out.append('\n## 2. Regenerated report against the full rerun (logs/spot_checks.out)\n')
out.append(f"Byte identical: {new_b == log_b}. Lines {len(new)} and {len(log)}. Numeric tokens compared: {k}; beyond 1e-9: {len(b2)}; largest relative difference: {nstr(w, 3)}.\n")

# ---------------------------------------------------------------- parse the regenerated report
NAMES = {
    "c = zeta'(1/2)/zeta(1/2)": 'c1', "c = (log pi - psi(1/4))/2": 'c2', "c = -theta'(0)": 'c3',
    "b = (euler + log pi)/2": 'b', "tau_0, the minimum of theta": 'tau0', "theta(tau_0)": 'theta_tau0',
    "max M = -theta(tau_0)/pi": 'maxM', "tau_pi with theta(tau_pi) = -pi": 'taupi',
    "M(gamma_1^2 -) = -theta(gamma_1)/pi": 'Mg1', "f(0+) = -zeta(1/2)/4": 'f0',
    "f'(0+) = -zeta(1/2)/16": 'fp0', "tail constant c/(2 sqrt pi)": 'tailc',
    "F(0.003)": 'F003', "F(0.0035)": 'F0035', "F(0.004)": 'F004', "x_*, the zero of W_0": 'xstar',
    "P(0.004), the small-x prime bound": 'P004', "W_0(0.1), closed form vs direct": 'W0_01',
    "w(0.001)": 'w0001', "w(0.01)": 'w001', "w(0.1)": 'w01', "w(1)": 'w1', "w(10)": 'w10',
    "min S over first 300, at gamma_289": 'minS', "max S over first 300, at gamma_213": 'maxS',
    "b_2 = (q^2+q)/2": 'b2', "b_3 = (q^3-q)/3": 'b3', "b_4": 'b4', "b_6": 'b6',
    "negative b_n for 2<=n<=60, p=2": 'neg2', "negative b_n for 2<=n<=60, p=3": 'neg3',
    "negative b_n for 2<=n<=60, p=5": 'neg5', "negative b_n for 2<=n<=60, p=7": 'neg7',
}
V, S = {}, {}          # values (mpf) and the strings they were printed as
refs = []              # (label, printed value, reference string printed by the script)


def put(key, s):
    V[key], S[key] = mpf(s), s


NUMRE = r'([-+]?[\d.]+(?:e[-+]?\d+)?)'
flip_ks, mids, low_ks = [], {}, []
sw = 0
for line in new:
    # lines printed by check(): f"{name:52s} {g:>22s}   paper: {claimed}"
    m = re.match(r'^(.{52}) (.{22})   paper: (.*)$', line)
    if m and m.group(1).strip() in NAMES:
        name, val = m.group(1).strip(), m.group(2).strip()
        refs.append((name, val, m.group(3).strip()))
        put(NAMES[name], val)
        continue
    m = re.match(r'kappa\(rho_(\d)\): xi form\s+\(' + NUMRE + r' ([-+]) ' + NUMRE + r'j\)\s+closed form\s+\(' + NUMRE + r' ([-+]) ' + NUMRE + r'j\)\s+paper: (\S+)', line)
    if m:
        k = m.group(1)
        put(f'k{k}re', m.group(2)); put(f'k{k}im', ('-' if m.group(3) == '-' else '') + m.group(4))
        put(f'c{k}re', m.group(5)); put(f'c{k}im', ('-' if m.group(6) == '-' else '') + m.group(7))
        refs.append((f'kappa(rho_{k})', f'{m.group(2)}{m.group(3)}{m.group(4)}i', m.group(8)))
        continue
    m = re.search(r'sign flips found: (\d+)\s+paper: (.*)$', line)
    if m:
        put('flips', m.group(1)); refs.append(('sign flips', m.group(1), m.group(2).strip()))
        continue
    m = re.match(r'\s+gamma_(\d+) = ' + NUMRE + r'\s+kappa = \(' + NUMRE + r' ([-+]) ' + NUMRE + r'j\)', line)
    if m:
        k = m.group(1); flip_ks.append(int(k))
        put(f'g{k}', m.group(2)); put(f'k{k}re', m.group(3)); put(f'k{k}im', ('-' if m.group(4) == '-' else '') + m.group(5))
        continue
    m = re.search(r'argmax\s+' + NUMRE + r'\s+f\(argmax\)\s+' + NUMRE + r'\s+paper: (.*)$', line)
    if m:
        put('hm_argmax', m.group(1)); put('hm_max', m.group(2)); refs.append(('f_H maximum', f'{m.group(2)} at {m.group(1)}', m.group(3).strip()))
        continue
    m = re.match(r'h=\s*([\d.]+)\s+P\(H>h\) Monte Carlo\s+' + NUMRE + r'\s+closed\s+' + NUMRE, line)
    if m:
        t = {'0.1': '01', '1.0': '1', '2.0': '2'}[m.group(1)]
        put('mc' + t, m.group(2)); put('cl' + t, m.group(3))
        continue
    m = re.match(r'Im w in \(([\d.]+),([\d.]+)\) -> gamma in \(([\d.]+),([\d.]+)\): (\d+) zeros, indices (\d+)\.\.(\d+)\s+table: (\d+)', line)
    if m:
        sw += 1
        put(f'sw{sw}', m.group(5)); put(f'sw{sw}_lo', m.group(6)); put(f'sw{sw}_hi', m.group(7))
        refs.append((f'xi zeros in ({m.group(3)},{m.group(4)})', m.group(5), 'table: ' + m.group(8)))
        continue
    m = re.search(r'gamma_195 = ' + NUMRE + r' < 390, gamma_222 = ' + NUMRE + r' > 430', line)
    if m:
        put('g195', m.group(1)); put('g222', m.group(2))
        continue
    for pat, key in [(r'^log LHS closed form\s+' + NUMRE, 'zr_lhs'), (r'^series truncated at n<3000\s+' + NUMRE, 'zr_trunc'),
                     (r'^truncation error\s+' + NUMRE, 'zr_err'), (r'^weight carried at n=3000\s+' + NUMRE, 'zr_wt'),
                     (r'^least negative minimum on the grid step 0\.025\s+' + NUMRE, 'grid_worst'),
                     (r'largest midpoint in modulus with Re kappa < 0: ' + NUMRE, 'mid_neg_max'),
                     (r'zeros with S\(gamma\^-\) < -1 among the first 460: (\d+)', 'n_low'),
                     (r'checks failed: (\d+)', 'failed'), (r'^int f over \(0,inf\)\s+' + NUMRE, 'int_f')]:
        m = re.search(pat, line)
        if m:
            put(key, m.group(1))
            tail = re.search(r'(paper(?: now)?: .*|row\'s printed gap: .*|round 1 assumed .*)$', line)
            if tail:
                refs.append((key, m.group(1), tail.group(1)))
    m = re.match(r"min Re zeta'/zeta on \(0,60\], sigma=\s*([\d.]+)\s+" + NUMRE + r'\s+paper: (.*)$', line)
    if m:
        key = {'0.55': 'm055', '0.6': 'm06', '0.75': 'm075', '0.9': 'm09', '0.99': 'm099'}[m.group(1)]
        put(key, m.group(2)); refs.append((f'min Re zeta\'/zeta, sigma={m.group(1)}', m.group(2), 'paper: ' + m.group(3).strip()))
    for k, v in re.findall(r'k=(\d+) mid=' + NUMRE, line):
        put(f'mid{k}', v); mids[int(k)] = mpf(v); low_ks.append(int(k))

NS = dict(V)
NS.update(mpf=mpf, sqrt=sqrt, pi=pi, abs=abs, min=min, max=max, flip_ks=flip_ks, mids=mids, low_ks=low_ks)
res['parsed'] = {k: S[k] for k in sorted(S)}


def ulp_of(expr):
    """Printed precision of an expression: the sum of the last digit units of the report values it uses."""
    return sum((ulp(S[n]) for n in set(re.findall(r'[A-Za-z_]\w*', expr)) if n in S), mpf(0))


# ---------------------------------------------------------------- 3. book quotations
tex = {c: read(c).decode().splitlines() for c in CHAPTERS}
res['input_sha256'] = {f: hashlib.sha256(read(f)).hexdigest() for f in
                       CHAPTERS + ['archived_spot_checks_report.md', 'spot_checks_report.md', 'run_all_spot_checks.out', 'book_quotes.tsv']}
rows, n_ok, n_bad, n_vals = [], 0, 0, 0
for raw in open('book_quotes.tsv'):
    if raw.startswith('#') or not raw.strip():
        continue
    rid, f, sub, qty, bvals, exprs, mode = raw.rstrip('\n').split('\t')
    where = [i + 1 for i, s in enumerate(tex[f]) if sub in s]
    loc = f"{f.replace('.tex', '')}:{','.join(map(str, where))}" if where else f'{f}: NOT FOUND'
    if mode == 'expr':
        try:
            ok = bool(eval(exprs, dict(NS)))
        except Exception as e:  # noqa: BLE001
            ok, exprs = False, f'{exprs} raised {e!r}'
        rows.append((rid, loc, qty, bvals, exprs, 'agrees' if ok and where else 'DISAGREES' if where else 'NOT FOUND'))
        n_vals += 1; n_ok += ok and bool(where); n_bad += not (ok and where)
        continue
    assert len(bvals.split(',')) == len(exprs.split(',')), rid
    for bs, ex in zip(bvals.split(','), exprs.split(',')):
        R = eval(ex, dict(NS))
        uR = ulp_of(ex)
        B = mpf(bs)
        if mode == 'int':
            ok = int(R) == int(B)
            shown = str(int(R))
        elif mode == 'round':
            ok = abs(R - B) <= ulp(bs) / 2 + uR / 2 + mpf('1e-40')
            shown = S[ex] if ex in S else nstr(R, 12)
        elif mode == 'trunc':
            ok = (R * B > 0) and abs(B) - uR / 2 <= abs(R) < abs(B) + ulp(bs) + uR / 2
            shown = S[ex] if ex in S else nstr(R, 12)
        else:
            raise SystemExit(f'unknown mode {mode}')
        ok = ok and bool(where)
        n_vals += 1; n_ok += ok; n_bad += not ok
        rows.append((rid, loc, qty, bs, f'{shown}' + ('' if ex in S else f' ({ex})'),
                     'agrees' if ok else ('NOT FOUND' if not where else 'DISAGREES')))
res['book'] = {'values_checked': n_vals, 'agree': n_ok, 'disagree_or_missing': n_bad,
               'rows': [dict(zip(['id', 'where', 'quantity', 'book', 'rerun', 'verdict'], r)) for r in rows]}
out.append('\n## 3. Numbers the book quotes from the report\n')
out.append(f'{n_vals} quoted values or statements checked: {n_ok} agree at the printed precision, {n_bad} do not or were not found.\n')
out.append('A value agrees when the book value is the rerun value correctly rounded at the last digit the book prints (allowing half a unit of the last digit the report prints; the book prints S(gamma_289^-) and gamma_289 truncated in rz3, marked trunc). Line numbers refer to the chapter files in book_snapshot/ (sha256 in comparison.json).\n')
out.append('| id | where | quantity | book | rerun | verdict |')
out.append('|---|---|---|---|---|---|')
for r in rows:
    out.append('| ' + ' | '.join(str(x).replace('|', '\\|') for x in r) + ' |')

# ---------------------------------------------------------------- 4. stale reference strings
alltex = '\n'.join('\n'.join(v) for v in tex.values())
SPECIAL = {'1.99e-13': '1.99\\times10^{-13}', 'exactly three': 'three times', 'now: seven': 'Seven of the first',
           "row's printed gap: -5.00e-4": '5.00\\times10^{-4}', 'below -0.729': '0.729', 'round 1 assumed 1': None}
stale, looked = [], 0
for name, val, ref in refs:
    r = re.sub(r'^(paper|table)(?: now)?: ', '', ref).strip()
    key = next((k for k in SPECIAL if k in ref), None)
    if key is not None:
        needle = SPECIAL[key]
    elif name in ('b_2 = (q^2+q)/2', 'b_3 = (q^3-q)/3', 'b_4', 'b_6'):
        continue  # the script computes this reference itself from the formula in the text
    else:
        needle = r.split(' at ')[0].lstrip('-+')
    looked += 1
    if needle is None:
        stale.append((name, val, ref, 'wording of a reply to a referee; no number to look for'))
        continue
    if key is None and re.fullmatch(r'\d+', needle):
        hit = re.search(r'\$' + needle + r'\$|\b' + needle + r'\b', alltex) is not None
    else:
        hit = needle in alltex
    if not hit:
        stale.append((name, val, ref, f'"{needle}" does not occur in rz1 to rz4, appA or ch12'))
res['stale_references'] = stale
out.append('\n## 4. Reference values printed by spot_checks.py that are no longer in the book\n')
out.append(f'The script prints {len(refs)} reference strings; {len(refs) - looked} of them (b_2, b_3, b_4, b_6) are its own evaluations of formulas in the text, and the other {looked} were looked up in the chapters. {looked - len(stale)} occur in the current text; {len(stale)} do not:\n')
for s in stale:
    out.append(f'- {s[0]}: rerun {s[1]}, script prints "{s[2]}"; {s[3]}')

# ---------------------------------------------------------------- Monte Carlo lines
out.append('\n## 5. Monte Carlo lines\n')
for t, h in [('01', '0.1'), ('1', '1'), ('2', '2')]:
    p = V['cl' + t]
    se = sqrt(p * (1 - p) / 200000)
    out.append(f"- h={h}: Monte Carlo {S['mc' + t]}, closed form {S['cl' + t]}, difference {nstr(V['mc' + t] - p, 3)}, standard error {nstr(se, 3)}, z = {nstr((V['mc' + t] - p) / se, 3)}; the script's tolerance 0.004 is {nstr(mpf('0.004') / se, 3)} standard errors.")
res['failed_checks_reported_by_script'] = S.get('failed')
open('comparison.md', 'w').write('\n'.join(out) + '\n')
json.dump(res, open('comparison.json', 'w'), indent=1)
print('\n'.join(out))
sys.exit(0 if n_bad == 0 else 1)
