# Check for the reconstruction of data/sinh_bench.json (misc/repro/reconstruct/sinh_bench/sinh_bench.py).
# Not an archive script. In mpmath at 40 digits it computes
#  (1) the exact values of the truncated lattice sums that sinh_bench.py evaluates in double precision (atoms and
#      reflections +-1..+-K): mu_K(a,c] = L(a,c) - Z0(a,c) - T_K(a,c), where L is the mass of the full lattice Z
#      (1 for a cell (k-1/2,k+1/2], 1/2 for the strip (0,1/2]), Z0 the kernel of the missing atom 0, and T_K the
#      atoms |j|>K, summed in closed form with log Gamma;
#  (2) a check of that formula against direct arctangent sums at K=300;
#  (3) the closed forms of Section sec:ch3:bench (K infinite), which the text of ch03 and ch15 quotes;
#  (4) the value of K implied by each archived number;
#  (5) several double precision codings of the same sums, to see which reproduce the archived digits
#      ('pair' is the coding of sinh_bench.py);
#  (6) the values printed in ch03 and ch15, against the closed forms and against the archived (truncated) numbers;
# and compares the archived file with the reconstructed file (and the closed form file) number by number.
# Usage: python3 check_sinh_bench.py ARCHIVED.json RECONSTRUCTED.json [CLOSED.json]
import numpy as np, json, sys, hashlib, time, mpmath as mp
mp.mp.dps = 40
book_path, recon_path = sys.argv[1], sys.argv[2]
closed_path = sys.argv[3] if len(sys.argv) > 3 else None
for p in [book_path, recon_path] + ([closed_path] if closed_path else []):
    print('sha256', hashlib.sha256(open(p, 'rb').read()).hexdigest(), p)
book = json.load(open(book_path)); recon = json.load(open(recon_path))
EPS = [0.5, 0.1, 0.05, 0.01, 0.005]; N = 50; K0 = 400000; Q = ['celltv', 'mass']
assert [r['eps'] for r in book] == EPS and [r['eps'] for r in recon] == EPS
B = [mp.mpf(0)] + [mp.mpf(2*j+1)/2 for j in range(N+1)]      # 0, 1/2, 3/2, ..., 101/2
NAT = [0] + [1]*N
LFULL = [mp.mpf(1)/2] + [mp.mpf(1)]*N                         # full lattice Z: (0,1/2] has 1/2, each cell 1
def Z0(a, c, e): return (mp.atan(c/e) - mp.atan(a/e))/mp.pi
def T(a, c, e, K):
    # (1/pi) sum_{j>K} [arctan((c-j)/e)-arctan((a-j)/e)+arctan((c+j)/e)-arctan((a+j)/e)]
    #  = (1/pi) Im[lnG(K+1-a+ie) - lnG(K+1-c+ie) - lnG(K+1+a+ie) + lnG(K+1+c+ie)]
    lg = lambda s: mp.im(mp.loggamma(mp.mpc(K+1+s, e)))
    return (lg(-a) - lg(-c) - lg(a) + lg(c))/mp.pi
def cells(e, K):
    return [LFULL[j] - Z0(B[j], B[j+1], e) - (0 if K is None else T(B[j], B[j+1], e, K)) for j in range(N+1)]
def summ(cm): return dict(celltv=mp.fsum(abs(m-n) for m, n in zip(cm, NAT)), mass=mp.fsum(cm))
def rel(a, b): return abs(float(a) - float(b))/abs(float(b))
def compare(A, Bf, tol=1e-9):
    n = nok = nbit = 0; worst = (0.0, None)
    for i in range(len(EPS)):
        for q in ['eps'] + Q:
            a, b_ = A[i][q], Bf[i][q]; d = abs(a-b_)/abs(a); n += 1; nok += d <= tol; nbit += a == b_
            if d >= worst[0]: worst = (d, '%s at eps=%g' % (q, EPS[i]))
    return dict(numbers=n, agree=nok, tol=tol, bitwise_equal=nbit, max_rel=worst[0], where=worst[1])
out = dict(K0=K0)

# (2) formula against direct sums at K=300
t0 = time.time(); worst = mp.mpf(0)
for e0 in EPS:
    e = mp.mpf(repr(e0)); f = cells(e, 300)
    for j in range(N+1):
        a, c = B[j], B[j+1]
        d = mp.fsum(mp.atan((c-x)/e) - mp.atan((a-x)/e) for kk in range(1, 301) for x in (kk, -kk))/mp.pi
        worst = max(worst, abs(d - f[j]))
out['formula_vs_direct_K300_maxabs'] = float(worst)
print('(2) log Gamma formula against direct sums at K=300: max abs diff %.3e (%.1f s)' % (worst, time.time()-t0))

# (1), (3), (4)
rows = []; CF = []
for i, e0 in enumerate(EPS):
    e = mp.mpf(repr(e0))
    ex = summ(cells(e, K0)); cf = summ(cells(e, None))
    cf2 = dict(celltv=(2*mp.atan(2*e) - mp.atan(e/(N+mp.mpf(1)/2)))/mp.pi, mass=N + mp.atan(e/(N+mp.mpf(1)/2))/mp.pi)
    CF.append(cf2)
    row = dict(eps=e0)
    for q in Q:
        arch = book[i][q]
        fp = (summ(cells(e, K0+1))[q] - summ(cells(e, K0-1))[q])/2   # d value / dK, central difference
        row[q] = dict(archived=arch, reconstructed=recon[i][q], exact_trunc=mp.nstr(ex[q], 20), closed_form=mp.nstr(cf2[q], 20),
                      lattice_limit_minus_closed_form=float(abs(cf[q]-cf2[q])),
                      rel_archived_vs_reconstructed=rel(arch, recon[i][q]), rel_archived_vs_exact_trunc=rel(arch, ex[q]),
                      rel_reconstructed_vs_exact_trunc=rel(recon[i][q], ex[q]), rel_archived_vs_closed=rel(arch, cf2[q]),
                      dvalue_dK=float(fp), K_implied=float(K0 + (mp.mpf(arch) - ex[q])/fp), K_per_ulp=float(np.spacing(arch)/abs(fp)))
    row['mass_defect_over_eps'] = dict(archived=(book[i]['mass'] - N)/e0, exact=float((cf2['mass'] - N)/e))
    rows.append(row)
    print('(1,3,4) eps=%g' % e0)
    for q in Q:
        r = row[q]
        print('   %-6s archived %.17g  reconstructed %.17g  exact(K=4e5) %s  closed(K=inf) %s' % (q, r['archived'], r['reconstructed'], r['exact_trunc'], r['closed_form']))
        print('          rel: arch/recon %.2e  arch/exact %.2e  recon/exact %.2e  arch/closed %.2e  K implied %.3f (1 ulp = %.2g in K)' % (
            r['rel_archived_vs_reconstructed'], r['rel_archived_vs_exact_trunc'], r['rel_reconstructed_vs_exact_trunc'], r['rel_archived_vs_closed'], r['K_implied'], r['K_per_ulp']))
    print('   mass defect/eps: archived %.8f  exact %.8f  first order 1/(50.5 pi) = %.8f' % (row['mass_defect_over_eps']['archived'], row['mass_defect_over_eps']['exact'], 1/(50.5*np.pi)))
out['rows'] = rows
out['archived_max_rel_vs_exact_trunc'] = max(r[q]['rel_archived_vs_exact_trunc'] for r in rows for q in Q)
out['reconstructed_max_rel_vs_exact_trunc'] = max(r[q]['rel_reconstructed_vs_exact_trunc'] for r in rows for q in Q)
out['K_implied_range'] = [min(r[q]['K_implied'] for r in rows for q in Q), max(r[q]['K_implied'] for r in rows for q in Q)]
print('archived vs exact truncated sums: max rel %.2e; reconstructed vs exact: max rel %.2e; K implied in [%.3f, %.3f]' % (
    out['archived_max_rel_vs_exact_trunc'], out['reconstructed_max_rel_vs_exact_trunc'], *out['K_implied_range']))

# comparisons of the files
out['compare_archived_reconstructed'] = compare(book, recon)
print('compare archived vs reconstructed:', out['compare_archived_reconstructed'])
if closed_path:
    closed = json.load(open(closed_path))
    out['closed_file_max_rel_vs_mpmath'] = max(rel(closed[i][q], CF[i][q]) for i in range(len(EPS)) for q in Q)
    out['compare_archived_closed'] = compare(book, closed)
    print('closed form file vs mpmath closed forms: max rel %.2e' % out['closed_file_max_rel_vs_mpmath'])
    print('compare archived vs closed forms:', out['compare_archived_closed'])

# (6) values printed in the book
L50 = N + mp.mpf(1)/2; c50 = (4 - 1/L50)/mp.pi
printed = [('ch03 Sec. sec:ch3:bench, mass defect coefficient to first order', '0.00630', 1/(L50*mp.pi), None),
           ('ch03 Sec. sec:ch3:bench, mass defect at eps=0.5', '3.151e-3', CF[0]['mass'] - N, book[0]['mass'] - N),
           ('ch03 Sec. sec:ch3:bench, cell distance at eps=0.5', '0.4968', CF[0]['celltv'], book[0]['celltv']),
           ('ch03 Sec. sec:ch3:bench, cell distance at eps=0.1', '0.1250', CF[1]['celltv'], book[1]['celltv']),
           ('ch03 Sec. sec:ch3:bench, cell distance at eps=0.05', '0.0631', CF[2]['celltv'], book[2]['celltv']),
           ('ch03 Sec. sec:ch3:bench, cell distance at eps=0.01', '0.01267', CF[3]['celltv'], book[3]['celltv']),
           ('ch03 Sec. sec:ch3:bench, cell distance at eps=0.005', '0.00633', CF[4]['celltv'], book[4]['celltv']),
           ('ch03 Sec. sec:ch3:bench, limit constant (4-1/50.5)/pi', '1.267', c50, book[4]['celltv']/0.005),
           ('ch03 Sec. sec:ch3:bench, per atom 1.267/50', '0.0253', c50/50, book[4]['celltv']/0.005/50),
           ('ch15 Sec. cells, constant for N=50', '1.26694', c50, book[4]['celltv']/0.005),
           ('ch15, first cell (0,3/2]', '0.418', (mp.mpf(4)/3 - 1/L50)/mp.pi, None),
           ('ch15 opening, constant', '1.27', c50, None),
           ('ch03 Exercise 5, mass defect coefficient', '0.00630', 1/(L50*mp.pi), (book[4]['mass'] - N)/0.005)]
out['printed'] = []
print('(6) printed values: printed | closed form | agrees | archived (truncated) data | agrees')
for where, s, v, a in printed:
    m = s.split('e'); dec = len(m[0].split('.')[1]); sc = 10.0**int(m[1]) if len(m) > 1 else 1.0
    ok = abs(float(v)/sc - float(m[0])) <= 0.5*10.0**-dec + 1e-15
    oka = None if a is None else abs(a/sc - float(m[0])) <= 0.5*10.0**-dec + 1e-15
    out['printed'].append(dict(where=where, printed=s, closed_form=float(v), closed_form_rounds_to_printed=bool(ok),
                               archived=a, archived_rounds_to_printed=oka))
    print('   %-62s %-9s %.8g %-5s %s %s' % (where, s, float(v), ok, '-' if a is None else '%.8g' % a, '' if oka is None else oka))

# (5) double precision codings of the same sums
t0 = time.time()
k = np.arange(1, K0+1, dtype=float)
orders = {'k_then_minus_k': np.concatenate([k, -k]), 'sorted': np.concatenate([-k[::-1], k]), 'interleaved': np.column_stack([k, -k]).ravel()}
bf = np.concatenate([[0.0], np.arange(N+1)+0.5]); nat = np.array(NAT)
def v_diff(x, e):
    cm = np.array([np.sum(np.arctan((c-x)/e) - np.arctan((a-x)/e)) for a, c in zip(bf[:-1], bf[1:])])/np.pi
    return np.sum(np.abs(cm-nat)), np.sum(cm)
def v_atan2(x, e):
    cm = np.array([np.sum(np.arctan2(e*(c-a), e*e+(c-x)*(a-x))) for a, c in zip(bf[:-1], bf[1:])])/np.pi
    return np.sum(np.abs(cm-nat)), np.sum(cm)
def v_cum(x, e):
    M = np.array([np.sum(np.arctan((t-x)/e)) for t in bf])/np.pi; cm = np.diff(M)
    return np.sum(np.abs(cm-nat)), M[-1]-M[0]
def v_cum0(x, e):
    M = np.array([np.sum(np.arctan((t-x)/e) + np.arctan(x/e)) for t in bf])/np.pi; cm = np.diff(M)
    return np.sum(np.abs(cm-nat)), M[-1]
def v_pair(x, e):
    kk = x[x > 0]; M = np.array([np.sum(np.arctan((t-kk)/e) + np.arctan((t+kk)/e)) for t in bf])/np.pi; cm = np.diff(M)
    return np.sum(np.abs(cm-nat)), M[-1]
variants = {}
for vn, vf in [('pair', v_pair), ('diff', v_diff), ('atan2', v_atan2), ('cum', v_cum), ('cum0', v_cum0)]:
    for on, x in orders.items():
        if vn == 'pair' and on != 'k_then_minus_k': continue
        vals = []; wa = 0; we = 0; nbit = 0; nok = 0
        for i, e0 in enumerate(EPS):
            tv, ms = vf(x, e0)
            for q, v in [('celltv', float(tv)), ('mass', float(ms))]:
                vals.append(v); a = book[i][q]; ex = float(rows[i][q]['exact_trunc'])
                wa = max(wa, abs(v-a)/abs(a)); we = max(we, abs(v-ex)/abs(ex)); nbit += (v == a); nok += abs(v-a)/abs(a) <= 1e-9
        name = 'pair' if vn == 'pair' else '%s/%s' % (vn, on)
        variants[name] = dict(values=vals, max_rel_vs_archived=wa, max_rel_vs_exact=we, bitwise_equal_archived=nbit, agree_1e9=nok)
        print('(5) %-22s max rel vs archived %.2e, vs exact %.2e, within 1e-9 of archived %2d of 10, bitwise equal %d of 10' % (name, wa, we, nok, nbit), flush=True)
out['variants'] = variants
print('(5) took %.1f s' % (time.time()-t0))
json.dump(out, open('check_sinh_bench.json', 'w'), indent=1)
