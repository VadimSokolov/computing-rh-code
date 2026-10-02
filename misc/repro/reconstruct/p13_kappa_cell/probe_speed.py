# Speed probe of acb.zeta_zeros at larger heights on Hopper (for the tail check of kappa_cell); not part of the book.
import time
from flint import acb, ctx
for prec in (64, 128):
    ctx.prec = prec
    for n in (1701, 10001, 30001, 100001):
        t0 = time.time(); z = acb.zeta_zeros(n, 100); dt = time.time() - t0
        print('prec %d: zeros %d..%d in %.2fs (%.1f ms/zero), last %s' % (prec, n, n + 99, dt, 10 * dt, z[-1].imag.str(20, radius=True)), flush=True)
