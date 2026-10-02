# Added for the book (not part of the authors' package): the twisted Riesz function R_h(x) at h = 50 and x = 1e10, computed
# with riesz_large.py (K = 200) and compared with the explicit formula (ri:eq:twexplicit) over 300 zeros, for the sentence on
# the twisted option of riesz_large.py in Section ri:sec:filters of Chapter ch:riesz.
# riesz_large.py is imported unchanged. The output check_large_twisted.json records:
#   main      R_50(1e10) = riesz(1e10, 200, 50) and its distance |R_h(x) - explicit|/x^(1/4) to the explicit formula;
#   explicit  the explicit formula over the first 300 zeros (both signs of gamma) and three trivial terms, as archived by
#             twisted.py (twisted.json, whose last point is x = 1e10), and recomputed in Arb ball arithmetic from Arb's
#             rigorously isolated zeros, with the sum over the zeros 301 to 400 as a truncation check;
#   rounding  riesz_twisted rounds 1/zeta(2 + 2ih) to double before it enters the tail; the exact effect of that rounding;
#   K_scan    the same computation with K = 100, 200, 400, 800, 1600, which shows the Mobius tail remainder decaying in K;
#   double    the same computation with numpy's longdouble replaced by double, which is what riesz_large.py computes where
#             numpy's longdouble is double (on Apple silicon, for example); the value quoted before this check came from such a run;
#   control   the untwisted R(1e10) against the archived check_large.json (the same code with T = 0).
# The extended precision of riesz_large.py needs numpy's longdouble to be the x87 80 bit type (x86_64 Linux); the script
# records its machine epsilon. Run it from code/ after setup.sh: it reads twisted.json and check_large.json.
import json, math, platform, socket, time
import numpy as np, mpmath as mp
from flint import arb, acb, acb_series, ctx
import riesz_large as rl

X, H, K = 1e10, 50.0, 200
X4 = X ** 0.25
QUOTED = complex(37.441418, -26.042497)            # the value in Section ri:sec:filters before this check
dev = lambda a, b: abs(a - b) / X4                 # the measure of the text: |R_h(x) - explicit| / x^(1/4)
cx = lambda z: [z.real, z.imag]
ld = np.finfo(np.longdouble)
out = dict(x=X, h=H, K=K, N=int(K * np.sqrt(X)), quoted_before=cx(QUOTED))
out['platform'] = dict(host=socket.gethostname(), machine=platform.machine(), python=platform.python_version(),
                       numpy=np.__version__, mpmath=mp.__version__, longdouble_eps=float(ld.eps),
                       longdouble_mantissa_bits=int(ld.nmant), extended=bool(ld.eps < 1e-18))
def save():
    json.dump(out, open('check_large_twisted.json', 'w'), indent=1)
print('platform', out['platform'], flush=True)
if not out['platform']['extended']:
    print('WARNING: numpy longdouble is double here, so the main value is a double precision computation', flush=True)

# explicit formula (ri:eq:twexplicit) as archived by twisted.py, and twisted.py's own Mobius series (N = 1e7, double)
tw = json.load(open('twisted.json'))['series']
assert tw['x'][-1] == X
P_tw = complex(tw['pre'][-1], tw['pim'][-1])
S_tw = complex(tw['re'][-1], tw['im'][-1])

# the explicit formula again in Arb ball arithmetic, from Arb's zeros
t0 = time.time()
ctx.prec = 128
LX = (arb(10) ** 10).log()
ih = acb(0, arb(H))
def zterm(r):
    d = acb_series([r, 1], prec=2).zeta().coeffs()[1]                  # zeta'(rho)
    return (1 - r / 2 + ih).gamma() / (2 * d) * (r / 2 * LX).exp()
def tterm(k):                                                           # zeta'(-2k) = (-1)^k (2k)! zeta(2k+1) / (2 (2 pi)^(2k))
    d = (-1) ** k * arb(math.factorial(2 * k)) * arb(2 * k + 1).zeta() / (2 * (2 * arb.pi()) ** (2 * k))
    return (1 + k + ih).gamma() / (2 * d) * (-k * LX).exp()
zs = acb.zeta_zeros(1, 400)
zt = [zterm(z) + zterm(acb(z.real, -z.imag)) for z in zs]               # rho and its conjugate
tt = [tterm(k) for k in (1, 2, 3)]
P = sum(zt[:300], acb(0)) + sum(tt, acb(0))
P_arb = complex(float(P.real.mid()), float(P.imag.mid()))
out['explicit'] = dict(
    twisted_json=cx(P_tw), arb=cx(P_arb),
    arb_digits=[P.real.str(25, radius=False), P.imag.str(25, radius=False)],
    arb_radius=[float(P.real.rad()), float(P.imag.rad())],
    arb_vs_twisted_json=dev(P_arb, P_tw),
    ordinate_300=float(zs[299].imag.mid()),
    zeros_301_400_abs_upper=float(abs(sum(zt[300:], acb(0))).upper()),
    trivial_abs_upper=[float(abs(t).upper()) for t in tt],
    twisted_py_series=cx(S_tw), twisted_py_series_dev=dev(S_tw, P_arb),
    seconds=time.time() - t0)
print('explicit formula: twisted.json', P_tw, ' Arb', out['explicit']['arb_digits'], ' radius', out['explicit']['arb_radius'],
      ' |Arb - twisted.json|/x^(1/4) %.2e' % out['explicit']['arb_vs_twisted_json'], ' zeros 301-400 %.1e' % out['explicit']['zeros_301_400_abs_upper'], flush=True)
save()

# main computation: riesz_large.py as archived, K = 200
t0 = time.time()
R = rl.riesz(X, K, H)
secs = time.time() - t0
out['main'] = dict(R=cx(R), abs_over_x14=abs(R) / X4, seconds=secs,
                   dev_vs_twisted_json=dev(R, P_tw), dev_vs_arb=dev(R, P_arb),
                   quoted_before_dev_vs_twisted_json=dev(QUOTED, P_tw), quoted_before_minus_R=dev(QUOTED, R))
print('R_50(1e10) = %r  |R|/x^(1/4) = %.6e  (%.1fs)' % (R, abs(R) / X4, secs), flush=True)
print('  |R - explicit|/x^(1/4): %.3e (twisted.json), %.3e (Arb);  quoted value: %.3e from twisted.json, %.3e from R'
      % (out['main']['dev_vs_twisted_json'], out['main']['dev_vs_arb'], out['main']['quoted_before_dev_vs_twisted_json'],
         out['main']['quoted_before_minus_R']), flush=True)
save()

# the rounding of 1/zeta(2 + 2ih) to double in riesz_twisted enters R_h(x) as x e^{ih log x} delta
mp.mp.dps = 40
iz = 1 / mp.zeta(2 + 2j * mp.mpf(H))
delta = mp.mpc(complex(iz)) - iz
eff = complex(mp.mpf(X) * mp.expj(mp.mpf(H) * mp.log(X)) * delta)
out['rounding'] = dict(inv_zeta=cx(complex(iz)), delta=cx(complex(delta)), effect=cx(eff), effect_over_x14=abs(eff) / X4,
                       R_corrected=cx(R - eff), dev_corrected=dev(R - eff, P_arb))
print('rounding of 1/zeta(2+2ih): effect %.3e in |R|/x^(1/4); corrected R deviates by %.3e'
      % (out['rounding']['effect_over_x14'], out['rounding']['dev_corrected']), flush=True)
save()

# control: the untwisted R(1e10) against the archived check_large.json
t0 = time.time()
R0 = rl.riesz(X)
arch = [r for r in json.load(open('check_large.json')) if r['x'] == X][0]
out['control_untwisted'] = dict(R=R0, archived=arch['R'], identical=R0 == arch['R'], seconds=time.time() - t0)
print('control R(1e10) = %r, archived %r, identical %s' % (R0, arch['R'], R0 == arch['R']), flush=True)
save()

# the same twisted computation with numpy's longdouble replaced by double
class _Double:
    longdouble, clongdouble = np.float64, np.complex128
    def __getattr__(self, name):
        return getattr(np, name)
saved = rl.np, rl.LD
rl.np, rl.LD = _Double(), np.float64
try:
    t0 = time.time()
    R64 = rl.riesz(X, K, H)
    s64 = time.time() - t0
finally:
    rl.np, rl.LD = saved
out['double'] = dict(R=cx(R64), dev=dev(R64, P_arb), minus_R=dev(R64, R), quoted_before_minus_double=dev(QUOTED, R64), seconds=s64)
print('double precision variant: %r  |.-explicit|/x^(1/4) %.3e  |.-R|/x^(1/4) %.3e  (%.1fs)'
      % (R64, out['double']['dev'], out['double']['minus_R'], s64), flush=True)
save()

# K scan: the tail remainder is of order K^(-7/2) in R_h(x)/x^(1/4) under square root cancellation
out['K_scan'] = []
for k in (100, 200, 400, 800, 1600):
    if k == K:
        Rk, s = R, secs
    else:
        t0 = time.time()
        Rk = rl.riesz(X, k, H)
        s = time.time() - t0
    row = dict(K=k, N=int(k * np.sqrt(X)), R=cx(Rk), dev=dev(Rk, P_arb), dev_corrected=dev(Rk - eff, P_arb),
               K_power=k ** -3.5, seconds=s)
    out['K_scan'].append(row)
    print('K %5d  N %.1e  R %r  dev %.3e  corrected %.3e  K^-3.5 %.1e  (%.1fs)'
          % (k, row['N'], Rk, row['dev'], row['dev_corrected'], row['K_power'], s), flush=True)
    save()
print(json.dumps(out, indent=1))
