# Independent check of the zero gaps used in Proposition bp:prop:v1min: zeros n in chunk TASK of NTASK,
# n <= NMAX, with mpmath.zetazero (Riemann-Siegel with Gram point bookkeeping), compared with Arb's zeros.
import json, sys
import mpmath as mp
from flint import acb, ctx
TASK, NTASK, NMAX = int(sys.argv[1]), int(sys.argv[2]), 1850
lo = 1 + TASK * NMAX // NTASK; hi = (TASK + 1) * NMAX // NTASK
mp.mp.dps = 20; ctx.prec = 80
zm = [float(mp.zetazero(n).imag) for n in range(lo, hi + 1)]
za = [float(z.imag.mid()) for z in acb.zeta_zeros(lo, hi - lo + 1)]
out = dict(task=TASK, n_lo=lo, n_hi=hi, mpmath=zm, arb=za, max_abs_diff=max(abs(a - b) for a, b in zip(zm, za)))
json.dump(out, open('zeros_%03d.json' % TASK, 'w'))
print(TASK, lo, hi, out['max_abs_diff'])
