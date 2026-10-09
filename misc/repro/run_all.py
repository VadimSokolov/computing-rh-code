# Run the book's scripts in dependency order, in parallel where possible, and log each run.
# Layout under BASE, the environment variable RH_REPRO_BASE or by default /scratch/vsokolov/rh_book_repro: the scripts of
# code/ and misc/repro/zeros.py run in run/; the five programs of verification_reciprocal_zeta/ (enclosures.py,
# spot_checks.py, certify_k9_k2.py, harmonic_decay.py, sawtooth_truncation.py) run in verif/;
# code_bosch/verify_bosch_xi_v7.py runs in bosch/; and the scripts of ferro/ (Chapter 13) run in ferro/, with a copy
# of data/zeros_1700.txt. compare.py then compares the outputs with the archived ones, copied to
# book_data/ (data/), book_verif/ (verification_reciprocal_zeta/) and book_bosch/ (code_bosch/). Not run here: the
# programs of verification_reciprocal_zeta/torus/, which run as Slurm arrays by the steps of torus/README.md, and the
# packages of Part VI and of the Riesz program, which have their own drivers, run6.py and run_riesz.py.
import subprocess, time, json, os, sys
from concurrent.futures import ThreadPoolExecutor
BASE = os.environ.get('RH_REPRO_BASE', '/scratch/vsokolov/rh_book_repro')
RUN = f'{BASE}/run'; LOG = f'{BASE}/logs'
PY = '/scratch/vsokolov/rh_book_repro/venv/bin/python3'
for d in (LOG, f'{RUN}/fig', f'{RUN}/book/fig'):  # most plotting scripts do not create their folders
    os.makedirs(d, exist_ok=True)
TO = 24 * 3600
EPS = ['1.0', '0.5', '0.25', '0.1', '0.05', '0.02', '0.01', '0.005']
T = {}
def task(name, cmd, deps=(), cwd=RUN):
    T[name] = dict(cmd=cmd, deps=list(deps), cwd=cwd)
task('zeros', [PY, '-u', 'zeros.py'])
task('part1', [PY, '-u', 'part1.py'])
for e in EPS: task(f'part2_{e}', [PY, '-u', 'part2.py', e])
P2 = [f'part2_{e}' for e in EPS]
task('analyse', [PY, '-u', 'analyse.py'], ['zeros', 'part1'] + P2)
task('turing_ig', [PY, '-u', 'turing_ig.py'], P2)
task('deconv', [PY, '-u', 'deconv.py'], ['zeros', 'part2_0.1'])
task('clockgrid', [PY, '-u', 'clockgrid.py'], ['part1'])
task('hazard_spec', [PY, '-u', 'hazard_spec.py'], ['part1'])
task('rsboundary', [PY, '-u', 'rsboundary.py'], ['hazard_spec'])
task('hitting', [PY, '-u', 'hitting.py'], ['hazard_spec'])
task('hitting2', [PY, '-u', 'hitting2.py'], ['hazard_spec'])
task('volterra', [PY, '-u', 'volterra.py'], ['hazard_spec', 'clockgrid', 'rsboundary'])
task('mc', [PY, '-u', 'mc.py'], ['rsboundary', 'hazard_spec'])
task('mc2', [PY, '-u', 'mc2.py'], ['volterra', 'rsboundary', 'hazard_spec'])
task('dbn', [PY, '-u', 'dbn.py'], ['part1'])
task('dbnplot', [PY, '-u', 'dbnplot.py'], ['dbn'])
task('zeros_hp', [PY, '-u', 'zeros_hp.py'], ['zeros', 'part1'])
task('arith', [PY, '-u', 'arith.py'], ['zeros_hp'])
for s in ['cm_orders', 'firstfail', 'gaussweil', 'phasefree', 'vdpair']:
    task(s, [PY, '-u', f'{s}.py'], ['zeros_hp'])
task('firstfail_arb', [PY, '-u', 'firstfail_arb.py'], ['firstfail'])
# ball arithmetic certificates that read firstfail.json and compute their zeros with Arb
task('gsmooth_arb', [PY, '-u', 'gsmooth_arb.py'], ['firstfail'])
task('phasefree_arb', [PY, '-u', 'phasefree_arb.py'], ['firstfail'])
for s in ['hankel', 'kenttail', 'kenttail2', 'kmono', 'perturb', 'selberg', 'spacing']:
    task(s, [PY, '-u', f'{s}.py'], ['zeros'])
for s in ['dbn_delta', 'hcm', 'idsd_check', 'idsd_height', 'lchi12', 'lehmer', 'lehmer_dbn', 'lfun', 'li', 'lsin2', 'maxlaw', 'primeside', 'saddle', 'sensitivity', 'sinh_bench', 'turing1000', 'williams_mc', 'tilted']:
    task(s, [PY, '-u', f'{s}.py'])
# ttest.py and ttest2.py take the command lines given in their headers (misc/repro/reconstruct/ttest2/)
task('ttest', [PY, '-u', 'ttest.py', '0.1', '200', '0.004', '4'])
task('ttest2_0.10', [PY, '-u', 'ttest2.py', '0.10', '200', '0.004', '4', '200', '500', '1000'])
task('ttest2_0.02', [PY, '-u', 'ttest2.py', '0.02', '220', '0.001', '4.2', '1000', '2000', '5000'])
task('sinh_bench_closed', [PY, '-u', 'sinh_bench.py', 'closed'])
task('hcm_arb', [PY, '-u', 'hcm_arb.py'], ['hcm'])
task('lplots', [PY, '-u', 'lplots.py'], ['lsin2', 'lfun'])
task('plots', [PY, '-u', 'plots.py'], ['analyse', 'hcm', 'primeside'])
task('plots2', [PY, '-u', 'plots2.py'], ['hazard_spec', 'clockgrid', 'rsboundary', 'mc'])
task('threezones', [PY, '-u', 'threezones.py'], ['volterra', 'rsboundary', 'mc2'])
task('plots3', [PY, '-u', 'plots3.py'], ['zeros', 'arith', 'primeside', 'analyse', 'turing_ig', 'hankel', 'spacing', 'sinh_bench'])
task('plots4', [PY, '-u', 'plots4.py'], ['analyse'])
task('deconvplot', [PY, '-u', 'deconvplot.py'], ['deconv'])
task('lehmerflowplot', [PY, '-u', 'lehmerflowplot.py'])
for s in ['enclosures', 'spot_checks', 'certify_k9_k2', 'harmonic_decay', 'sawtooth_truncation']:
    task(s, [PY, '-u', f'{s}.py'], cwd=f'{BASE}/verif')
task('bosch_v7', [PY, '-u', 'verify_bosch_xi_v7.py'], ['zeros'], cwd=f'{BASE}/bosch')
FERRO = f'{BASE}/ferro'
for s in ['fm_polya', 'fm_interval', 'fm_hcm', 'fm_sizebias']:
    task(s, [PY, '-u', f'{s}.py'], cwd=FERRO)
task('fm_hermite', [PY, '-u', 'fm_hermite.py'] + [str(k) for k in (1, 2, 3, 5, 10, 16, 17, 20, 30, 50, 100, 150, 200, 300, 400, 600, 800)], cwd=FERRO)
task('fm_hermite_all', [PY, '-u', 'fm_hermite.py', 'all', '200'], cwd=FERRO)
task('fm_maxwell', [PY, '-u', 'fm_maxwell.py', '10000000'], cwd=FERRO)
# ball arithmetic certificates of Chapter 13; fm_qcert.py zeros writes fm_zeros_hp.txt (800 zeros and Var X from Arb)
task('fm_cutoff', [PY, '-u', 'fm_cutoff.py'], cwd=FERRO)
task('fm_qcert_zeros', [PY, '-u', 'fm_qcert.py', 'zeros', '800', '9600'], cwd=FERRO)
task('fm_qcert_table', [PY, '-u', 'fm_qcert.py', 'table'], ['fm_qcert_zeros', 'fm_hermite'], cwd=FERRO)
task('fm_qcert_all', [PY, '-u', 'fm_qcert.py', 'all', '200'], ['fm_qcert_zeros'], cwd=FERRO)
task('fm_qcert_minq', [PY, '-u', 'fm_qcert.py', 'minq', '1-200,300,400,600,800'], ['fm_qcert_zeros'], cwd=FERRO)
task('fm_tvcert', [PY, '-u', 'fm_tvcert.py'], ['fm_qcert_zeros'], cwd=FERRO)
task('fm_intcert', [PY, '-u', 'fm_intcert.py'], cwd=FERRO)
status = {}
def run(name):
    t = T[name]
    for d in t['deps']:
        if status.get(d, {}).get('rc') != 0:
            status[name] = dict(rc=None, skipped=f'dependency {d} failed'); return
    t0 = time.time()
    with open(f'{LOG}/{name}.out', 'w') as out:
        try:
            rc = subprocess.run(t['cmd'], cwd=t['cwd'], stdout=out, stderr=subprocess.STDOUT, timeout=TO).returncode
        except subprocess.TimeoutExpired:
            rc = 'timeout'
    status[name] = dict(rc=rc, seconds=round(time.time() - t0, 1))
    print(f'{name}: rc={rc} {status[name]["seconds"]}s', flush=True)
    json.dump(status, open(f'{BASE}/status.json', 'w'), indent=1)
# simple scheduler: repeatedly launch tasks whose dependencies are finished
from concurrent.futures import wait, FIRST_COMPLETED
pending = set(T); running = {}
with ThreadPoolExecutor(max_workers=int(os.environ.get('SLURM_CPUS_PER_TASK', '16'))) as ex:
    while pending or running:
        ready = [n for n in pending if all(d in status for d in T[n]['deps'])]
        for n in ready:
            pending.discard(n); running[ex.submit(run, n)] = n
        if not running: break
        done, _ = wait(list(running), return_when=FIRST_COMPLETED)
        for f in done: running.pop(f)
json.dump(status, open(f'{BASE}/status.json', 'w'), indent=1)
print('ALL DONE', flush=True)
