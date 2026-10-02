# Reproducibility run for Part VI of the book: the scripts of basepoint_one/ and sato_bondesson/.
# src/ is a copy of the repository folders (and data/); run/ is src/ without the archived JSON outputs,
# which are kept in ref/ for compare6.py. Scripts run from their own folders.
#   python3 run6.py setup      build run/ and ref/ afresh
#   python3 run6.py group K    run group K of the independent chains (one Slurm array task each)
#   python3 run6.py rest       run every task not yet done, in dependency order, in parallel
# Each task writes status/<task>.json with its return code and time.
import subprocess, time, json, os, shutil, glob, sys
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
BASE = '/scratch/vsokolov/rh_book_repro/part6'
SRC, RUN, REF, LOG, STAT = [f'{BASE}/{d}' for d in ('src', 'run', 'ref', 'logs', 'status')]
PY = '/scratch/vsokolov/rh_book_repro/venv/bin/python3'
TO = 24 * 3600
X, E, S, W, P, B = ('basepoint_one/xi', 'basepoint_one/eta', 'basepoint_one/scale_mixture',
                    'basepoint_one/wiener', 'basepoint_one/process_map', 'sato_bondesson')
FOLDERS = [X, E, S, W, P, B]

def setup():
    for d in (RUN, REF, STAT):
        shutil.rmtree(d, ignore_errors=True)
    os.makedirs(LOG, exist_ok=True); os.makedirs(STAT)
    for f in FOLDERS:
        shutil.copytree(f'{SRC}/{f}', f'{RUN}/{f}')
        os.makedirs(f'{REF}/{f}', exist_ok=True)
        for p in glob.glob(f'{RUN}/{f}/*.json') + glob.glob(f'{RUN}/{f}/*.npy'):
            name = os.path.basename(p)
            if name.startswith('g_'):
                continue  # the zeta ordinates are inputs
            shutil.copy(p, f'{REF}/{f}/{name}')
            if name.endswith('.json'):
                os.remove(p)
    shutil.copytree(f'{SRC}/data', f'{RUN}/data')

T = {}
def py(tag, folder, script, *args, deps=()):
    name = f'{tag}:{script[:-3]}' + ('_' + '_'.join(args) if args else '')
    T[name] = dict(cwd=f'{RUN}/{folder}', cmd=[PY, '-u', script, *args], deps=list(deps))

def stage():
    # inputs that the process map and the Sato note copy from the basepoint one note
    for f in ['alpha1.json', 'rs1.json']:
        shutil.copy(f'{RUN}/{X}/{f}', f'{RUN}/{P}/{f}')
    for f in ['improve.json', 'rs1.json']:
        shutil.copy(f'{RUN}/{X}/{f}', f'{RUN}/{B}/{f}')
    v = '/scratch/vsokolov/rh_book_repro/run/volterra.json'
    shutil.copy(v if os.path.exists(v) else f'{REF}/{P}/volterra.json', f'{RUN}/{P}/volterra.json')
    return 0
T['stage'] = dict(fn=stage, deps=['xi:compute', 'xi:rs1c', 'xi:improve'])

py('xi', X, 'compute.py'); py('xi', X, 'improve.py'); py('xi', X, 'levy.py'); py('xi', X, 'below.py')
py('xi', X, 'errformula.py', deps=['xi:below'])
py('xi', X, 'heat_primes.py')  # needs python-flint; 32 worker processes by default
py('xi', X, 'heat_thresholds.py')  # needs python-flint; 32 worker processes by default
py('xi', X, 'rs1.py'); py('xi', X, 'rs1c.py', deps=['xi:rs1']); py('xi', X, 'rs1d.py', deps=['xi:rs1c'])
py('xi', X, 'mc.py')
py('xi', X, 'figs.py', deps=['xi:compute', 'xi:mc']); py('xi', X, 'figs2.py', deps=['xi:improve'])
py('xi', X, 'figs3.py', deps=['xi:compute', 'xi:levy']); py('xi', X, 'figs4.py', deps=['xi:below'])
py('xi', X, 'figs5.py', deps=['xi:rs1d'])
py('eta', E, 'zeros12.py'); py('eta', E, 'base.py', deps=['eta:zeros12']); py('eta', E, 'arith.py', deps=['eta:base'])
py('eta', E, 'hit.py', deps=['eta:base']); py('eta', E, 'hit2.py', deps=['eta:hit'])
py('eta', E, 'figs.py', deps=['eta:arith', 'eta:hit2', 'xi:compute'])
py('sm', S, 'improve.py'); py('sm', S, 'comp.py'); py('sm', S, 'strip.py'); py('sm', S, 'dens2.py', deps=['sm:strip'])
py('sm', S, 'figs2.py', deps=['sm:comp', 'sm:dens2'])
py('wr', W, 'wiener.py'); py('wr', W, 'salem_riesz.py')
py('pm', P, 'comp.py', deps=['stage']); py('pm', P, 'comp2.py', deps=['stage']); py('pm', P, 'ou.py')
py('pm', P, 'bnd_alpha.py', '0.6'); py('pm', P, 'bnd_alpha.py', '0.75')
py('pm', P, 'figs.py', deps=['pm:ou', 'pm:bnd_alpha_0.6', 'pm:bnd_alpha_0.75', 'stage'])
py('sb', B, 'comp.py', deps=['xi:compute', 'xi:rs1c', 'xi:improve']); py('sb', B, 'figs.py', deps=['sb:comp', 'stage'])
py('sb', B, 'weil.py'); py('sb', B, 'weil.py', 'mp')

# independent chains, one Slurm array task each; everything else runs in the 'rest' job
GROUPS = [['xi:compute'], ['xi:improve'], ['xi:levy'], ['xi:below', 'xi:errformula'],
          ['xi:rs1', 'xi:rs1c', 'xi:rs1d'], ['xi:mc'],
          ['eta:zeros12', 'eta:base', 'eta:arith', 'eta:hit', 'eta:hit2'],
          ['sm:improve'], ['sm:comp'], ['sm:strip', 'sm:dens2'], ['wr:wiener'], ['wr:salem_riesz'],
          ['pm:ou'], ['pm:bnd_alpha_0.6'], ['pm:bnd_alpha_0.75'],
          ['sb:weil'], ['sb:weil_mp'], ['xi:heat_primes'], ['xi:heat_thresholds']]
assert all(n in T for g in GROUPS for n in g)

def status():
    st = {}
    for p in glob.glob(f'{STAT}/*.json'):
        st[os.path.basename(p)[:-5].replace('__', ':')] = json.load(open(p))
    return st

def record(name, st):
    json.dump(st, open(f'{STAT}/{name.replace(":", "__")}.json', 'w'))

env = dict(os.environ, OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='2', MKL_NUM_THREADS='2', MPLBACKEND='Agg')
def run(name):
    t = T[name]; st = status()
    for d in t['deps']:
        if st.get(d, {}).get('rc') != 0:
            record(name, dict(rc=None, skipped=f'dependency {d} failed or missing')); return
    t0 = time.time()
    if 'fn' in t:
        try:
            rc = t['fn']()
        except Exception as e:
            rc = repr(e)
    else:
        with open(f'{LOG}/{name.replace(":", "__")}.out', 'w') as out:
            try:
                rc = subprocess.run(t['cmd'], cwd=t['cwd'], stdout=out, stderr=subprocess.STDOUT, timeout=TO, env=env).returncode
            except subprocess.TimeoutExpired:
                rc = 'timeout'
    record(name, dict(rc=rc, seconds=round(time.time() - t0, 1), host=os.uname().nodename))
    print(f'{name}: rc={rc} {round(time.time() - t0, 1)}s', flush=True)

def rest():
    done = {n for n, s in status().items() if s.get('rc') == 0}
    pending = set(T) - done; running = {}
    with ThreadPoolExecutor(max_workers=int(os.environ.get('SLURM_CPUS_PER_TASK', '8'))) as ex:
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
    mode = sys.argv[1]
    if mode == 'setup':
        setup()
    elif mode == 'group':
        for n in GROUPS[int(sys.argv[2])]:
            run(n)
    elif mode == 'rest':
        rest()
    elif mode == 'ngroups':
        print(len(GROUPS))
