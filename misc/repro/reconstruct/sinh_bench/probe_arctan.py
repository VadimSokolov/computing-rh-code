# Exploration for the reconstruction of data/sinh_bench.json, not an archive script.
# sinh_bench.py (paired coding) reproduces six of the ten archived values bit for bit and the other four to 1.1e-12.
# Question: do the four come from the last bit of numpy's arctan, which on an AVX512 Intel node is the SVML
# routine and elsewhere the C library's atan? This runs the paired coding of sinh_bench.py in child processes
# with numpy's AVX512 (X86_V4) dispatch target on and off (NPY_DISABLE_CPU_FEATURES=X86_V4, the numpy 2 group name)
# and compares each result with the archived file.
# Usage: python3 probe_arctan.py ARCHIVED.json
import os, sys, json, subprocess, numpy as np
from numpy._core._multiarray_umath import __cpu_features__, __cpu_baseline__
child = r'''
import numpy as np, json
from numpy._core._multiarray_umath import __cpu_features__
try:
    from numpy.lib.introspect import opt_func_info
    info = opt_func_info(func_name='^arctan$', signature='float64')
except Exception as ex:
    info = repr(ex)
K = 400000; N = 50; b = np.concatenate([[0.0], np.arange(N+1)+0.5]); nat = np.array([0]+[1]*N)
k = np.arange(1, K+1, dtype=float); res = []
for eps in [0.5, 0.1, 0.05, 0.01, 0.005]:
    M = np.array([np.sum(np.arctan((t-k)/eps)+np.arctan((t+k)/eps)) for t in b])/np.pi; cm = np.diff(M)
    res.append(dict(eps=eps, celltv=float(np.sum(np.abs(cm-nat))), mass=float(M[-1])))
print(json.dumps(dict(numpy=np.__version__, enabled=sorted(f for f, v in __cpu_features__.items() if v), info=str(info), res=res)))
'''
book = json.load(open(sys.argv[1]))
avx512 = ' '.join(f for f, v in __cpu_features__.items() if v and f.startswith('AVX512'))
print('numpy', np.__version__, '; baseline', ' '.join(__cpu_baseline__), '; AVX512 features on this node:', avx512 or 'none', flush=True)
out = {}
for label, env in [('default dispatch', {}), ('X86_V4 (AVX512) target disabled', {'NPY_DISABLE_CPU_FEATURES': 'X86_V4'})]:
    if label != 'default dispatch' and not avx512: continue
    e = dict(os.environ); e.update(env)
    r = subprocess.run([sys.executable, '-c', child], env=e, capture_output=True, text=True)
    if r.returncode != 0: print(label, 'failed:', r.stderr[-2000:]); continue
    d = json.loads(r.stdout.strip().splitlines()[-1])
    cmp = []
    for i, row in enumerate(d['res']):
        for q in ['celltv', 'mass']:
            a = book[i][q]; v = row[q]
            cmp.append(dict(eps=row['eps'], q=q, archived=a, value=v, bitwise=(a == v), rel=abs(v-a)/abs(a)))
    nbit = sum(c['bitwise'] for c in cmp)
    print('%s: %d of 10 bitwise equal to the archived values, max rel %.2e' % (label, nbit, max(c['rel'] for c in cmp)))
    print('   arctan dispatch:', d['info'].replace('\n', ' | ')[:600])
    for c in cmp:
        if not c['bitwise']: print('   differs: %-6s eps=%-6g archived %.17g  here %.17g  rel %.2e' % (c['q'], c['eps'], c['archived'], c['value'], c['rel']))
    out[label] = dict(numpy=d['numpy'], env=env, arctan_dispatch=d['info'], bitwise_equal=nbit, rows=cmp)
json.dump(out, open('probe_arctan.json', 'w'), indent=1)
