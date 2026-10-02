# Probe of the python-flint 0.9 API on Hopper (zeta derivative, zeta zeros, digamma); not part of the book.
import flint, time
from flint import arb, acb, ctx
print('python-flint', flint.__version__)
print([n for n in dir(flint) if 'series' in n or 'poly' in n])
print('acb methods:', [m for m in dir(acb) if 'zeta' in m or 'gamma' in m or 'digamma' in m or 'log' in m])
if hasattr(flint, 'acb_series'):
    print('acb_series methods:', [m for m in dir(flint.acb_series) if not m.startswith('_')])
ctx.prec = 256
s = acb(0.5, 200)
try:
    z = flint.acb_series([s, 1], prec=2).zeta()
    print('acb_series zeta', z)
except Exception as e:
    print('acb_series zeta failed', e)
t0 = time.time(); zz = acb.zeta_zeros(1, 80); print('zeta_zeros(1,80) %.2fs' % (time.time() - t0), zz[0], zz[78], zz[79])
t0 = time.time(); z1 = acb.zeta_zero(1700); print('zeta_zero(1700) %.2fs' % (time.time() - t0), z1)
print('digamma', (s / 2).digamma())
x = arb(-1.5, 1e-10); print('abs', abs(x))
