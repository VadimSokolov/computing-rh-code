# RECONSTRUCTED SCRIPT, NOT THE AUTHORS' ORIGINAL. No script in the archive wrote data/sinh_bench.json, which
# code/plots3.py (section 12) reads for Figure fig:ch3:bench; this one reproduces it to 1.1e-12. Reconstructed in
# September 2026 (misc/repro/reconstruct/sinh_bench/REPORT.md) and added to code/ for the book.
# It follows the first draft of Section sec:ch3:bench (ch/ch03.tex at the baseline commit 023e266: "The cumulative
# mass and the Voronoi cell masses are sums of arctangents. On the window 0<theta<=50.5, with cells centred on the
# 50 atoms ...") and review U02 (misc/reviews/round1/U02-review.md, B3 and E6), which found that the archived numbers
# are the arctangent sums over the atoms 0<|k|<=4e5. The Thorin density of the sinh clock at resolution eps is
# sum_{k!=0} P_eps(theta-k), P_eps(x)=eps/(pi(eps^2+x^2)); here the lattice is cut at |k|<=K=400000, which is why
# the archived mass is 50.0031113 and not the exact 50+arctan(eps/50.5)/pi=50.0031515 at eps=0.5.
# Output sinh_bench.json, one row per eps:
#   celltv = |mu((0,1/2])-0| + sum_{k=1}^{50} |mu((k-1/2,k+1/2])-1|   (the strip (0,1/2] holds no atom)
#   mass   = mu((0,50.5])
# Usage: python3 sinh_bench.py [K]     K=400000 (default) reproduces the archived file
#        python3 sinh_bench.py closed  the closed forms of Section sec:ch3:bench (K infinite), in sinh_bench_closed.json
# Inputs: none (numpy only); the default run takes about 2 s.
# Reproduces: data/sinh_bench.json to 1.1e-12 relative. 11 of its 15 numbers agree bit for bit; celltv at eps = 0.1
# and 0.01 differs by 5.7e-14 and 1.12e-12 relative, and mass at eps = 0.1 and 0.01 by one and two units in the last
# place. These last bits depend on the rounding of the double precision sums, so other hardware agrees to about 1e-12.
# The closed forms, which the text quotes (Section sec:ch3:bench and Chapter ch:spacing), differ from the truncated
# lattice by at most 7.9e-5 relative; the option closed writes them to sinh_bench_closed.json and leaves
# sinh_bench.json alone.
import numpy as np, json, sys
arg = sys.argv[1] if len(sys.argv) > 1 else '400000'
N = 50                                                    # atoms 1..50 in the window 0<theta<=50.5
b = np.concatenate([[0.0], np.arange(N+1)+0.5])           # boundaries 0, 1/2, 3/2, ..., 50.5
natoms = np.array([0]+[1]*N)                              # atoms in (0,1/2] and in each cell (k-1/2,k+1/2]
res = []
for eps in [0.5, 0.1, 0.05, 0.01, 0.005]:
    if arg == 'closed':
        row = dict(eps=eps, celltv=float((2*np.arctan(2*eps)-np.arctan(eps/(N+0.5)))/np.pi), mass=float(N+np.arctan(eps/(N+0.5))/np.pi))
    else:
        k = np.arange(1, int(arg)+1, dtype=float)         # atoms k and their reflections -k, 1<=k<=K
        # cumulative mass mu((0,t]) = (1/pi) sum_k [arctan((t-k)/eps) + arctan((t+k)/eps)]
        M = np.array([np.sum(np.arctan((t-k)/eps)+np.arctan((t+k)/eps)) for t in b])/np.pi
        cm = np.diff(M)                                   # masses of (0,1/2] and of the cells (k-1/2,k+1/2]
        row = dict(eps=eps, celltv=float(np.sum(np.abs(cm-natoms))), mass=float(M[-1]))
    res.append(row); print(row, flush=True)
json.dump(res, open('sinh_bench_closed.json' if arg == 'closed' else 'sinh_bench.json', 'w'))
