# code/zeros_hp.py (September 2026, from misc/repro/reconstruct/p13_arith_exact): the ordinates of the first 1700
# zeros of zeta as 70 digit decimal strings, one per line, in zeros_1700.txt, which arith.py, gaussweil.py, firstfail.py, firstfail_arb.py and phasefree.py read. They replace, on the
# zero side of arith.py, the float64 lists g_1_400.npy, g_401_1000.npy, g_1001_1700.npy that the first version of
# arith.py read and that the other scripts still read (float(mpmath.zetazero(n).imag), rounded to double precision).
# Arb (python-flint acb.zeta_zeros) isolates and refines the zeros rigorously at 320 bits. Checks: every ball has
# real part 1/2 and imaginary radius below 1e-80; the ordinates increase; the strings reproduce the balls to a
# relative 1e-69 (the rounding to 70 significant digits). After writing zeros_1700.txt the script checks that the
# float64 lists agree to within half an ulp and that the 79 ordinates of part1.json (45 digits, the book's own
# Newton refined roots below height 200, refined at 50 digits) agree to 1e-41; both files must then be present.
import json, time, numpy as np
from flint import acb, arb, ctx

ctx.prec = 320
t0 = time.time()
Z = acb.zeta_zeros(1, 1700)
print('computed %d zeros in %.1f s' % (len(Z), time.time() - t0), flush=True)
assert len(Z) == 1700
half = arb(1) / 2
maxrad = max(float(z.imag.rad()) for z in Z)
assert all(z.real == half for z in Z), 'a real part is not exactly 1/2'
assert maxrad < 1e-80, maxrad
S = [z.imag.str(70, radius=False) for z in Z]
for s, z in zip(S, Z):
    assert abs(arb(s) - z.imag) < arb('1e-69') * z.imag, s
assert all(Z[k + 1].imag > Z[k].imag for k in range(1699)), 'ordinates not increasing'
open('zeros_1700.txt', 'w').write('\n'.join(S) + '\n')
print('max radius of the ordinates %.2e' % maxrad)
print('gamma_1    =', S[0])
print('gamma_1700 =', S[-1])

# the archived float64 lists
g = np.concatenate([np.load(f) for f in ['g_1_400.npy', 'g_401_1000.npy', 'g_1001_1700.npy']])
assert len(g) == 1700
err = [abs(float((arb(float(x)) - z.imag).mid())) for x, z in zip(g, Z)]
ulp = [float(np.spacing(x)) for x in g]
k = int(np.argmax(err))
print('float64 lists: max |g - gamma| = %.3e at n = %d (gamma = %.3f, ulp %.3e); max error/ulp = %.3f; |g_1 - gamma_1| = %.3e'
      % (err[k], k + 1, g[k], ulp[k], max(e / u for e, u in zip(err, ulp)), err[0]))
assert max(e / u for e, u in zip(err, ulp)) <= 0.5 + 1e-9

# the 79 roots of data/part1.json
R = json.load(open('part1.json'))['roots']
d = max(abs(float((arb(r) - z.imag).mid())) for r, z in zip(R, Z))
print('part1.json roots: %d compared, max |difference| = %.3e' % (len(R), d))
assert d < 1e-41
print('done in %.1f s' % (time.time() - t0))
