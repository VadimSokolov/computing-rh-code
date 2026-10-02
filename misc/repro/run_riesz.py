# Reproducibility run for the Riesz program package (riesz_program/, Chapters ch:riesz and ch:almost).
#   python3 run_riesz.py setup | group K | rest      (as run6.py)
# run/code holds the scripts with the zero lists; the archived JSON outputs go to ref/code for compare6.py.
# riesz_explicit.json has no script in the package and is kept as an input of figs_riesz.py.
import subprocess, time, json, os, shutil, glob, sys
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
BASE = '/scratch/vsokolov/rh_book_repro/riesz'
SRC, RUN, REF, LOG, STAT = [f'{BASE}/{d}' for d in ('src', 'run', 'ref', 'logs', 'status')]
PY = '/scratch/vsokolov/rh_book_repro/venv/bin/python3'
TO = 24 * 3600
C = f'{RUN}/code'

def setup():
    for d in (RUN, REF, STAT):
        shutil.rmtree(d, ignore_errors=True)
    os.makedirs(LOG, exist_ok=True); os.makedirs(STAT); os.makedirs(f'{REF}/code')
    shutil.copytree(f'{SRC}/code', C)
    for f in glob.glob(f'{SRC}/zeros/*.npy'):
        shutil.copy(f, C)
    for f in glob.glob(f'{SRC}/data/*.json'):
        shutil.copy(f, f'{REF}/code/')
    shutil.copy(f'{SRC}/data/riesz_explicit.json', C)

T = {}
def py(name, script, *args, deps=()):
    T[name] = dict(cmd=[PY, '-u', script, *args], deps=list(deps))
py('check_large', 'check_large.py'); py('riesz_large_tw', 'riesz_large.py', '1e10', '200', '50')
py('twisted', 'twisted.py'); py('sens', 'sens.py'); py('salem_riesz', 'salem_riesz.py'); py('riesz_thorin', 'riesz_thorin.py')
py('almost', 'almost.py'); py('more2', 'more2.py', deps=['almost']); py('high', 'high.py'); py('ggc_coef', 'ggc_coef.py')
py('killing', 'killing.py'); py('killing2', 'killing2.py', deps=['killing']); py('ks', 'ks.py')
py('figs_riesz', 'figs_riesz.py', deps=['salem_riesz', 'riesz_thorin'])
py('figs_strategy', 'figs_strategy.py', deps=['sens', 'twisted'])
# Added for the book: the certificates of Proposition bp:prop:v1min and Theorem ag:thm:explicit. They are in no group
# and run in the rest step (about 90 s and 10 s).
py('v1min', 'v1min.py'); py('explicit_bounds', 'explicit_bounds.py')
GROUPS = [['check_large'], ['riesz_large_tw'], ['twisted'], ['sens'], ['salem_riesz'], ['riesz_thorin'], ['almost', 'more2'], ['high'], ['ggc_coef'], ['killing', 'killing2'], ['ks']]

def status():
    return {os.path.basename(p)[:-5]: json.load(open(p)) for p in glob.glob(f'{STAT}/*.json')}
def record(name, st):
    json.dump(st, open(f'{STAT}/{name}.json', 'w'))
env = dict(os.environ, OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='2', MKL_NUM_THREADS='2', MPLBACKEND='Agg')
def run(name):
    t = T[name]; st = status()
    for d in t['deps']:
        if st.get(d, {}).get('rc') != 0:
            record(name, dict(rc=None, skipped=f'dependency {d} failed or missing')); return
    t0 = time.time()
    with open(f'{LOG}/{name}.out', 'w') as out:
        try:
            rc = subprocess.run(t['cmd'], cwd=C, stdout=out, stderr=subprocess.STDOUT, timeout=TO, env=env).returncode
        except subprocess.TimeoutExpired:
            rc = 'timeout'
    record(name, dict(rc=rc, seconds=round(time.time() - t0, 1), host=os.uname().nodename))
    print(f'{name}: rc={rc} {round(time.time() - t0, 1)}s', flush=True)
def rest():
    done = {n for n, s in status().items() if s.get('rc') == 0}
    pending = set(T) - done; running = {}
    with ThreadPoolExecutor(max_workers=int(os.environ.get('SLURM_CPUS_PER_TASK', '4'))) as ex:
        while pending or running:
            st = status(); busy = pending | set(running.values())
            ready = [n for n in sorted(pending) if all(d in st and d not in busy for d in T[n]['deps'])]
            for n in ready:
                pending.discard(n); running[ex.submit(run, n)] = n
            if not running:
                break
            fin, _ = wait(list(running), return_when=FIRST_COMPLETED)
            for f in fin:
                running.pop(f)
    for n in sorted(pending):
        record(n, dict(rc=None, skipped='never ready'))
    json.dump(status(), open(f'{BASE}/status.json', 'w'), indent=1)
    print('all done', flush=True)
if __name__ == '__main__':
    m = sys.argv[1]
    if m == 'setup': setup()
    elif m == 'group':
        for n in GROUPS[int(sys.argv[2])]: run(n)
    elif m == 'rest': rest()
