# Reconstruction helper (not an authors' script). Searches for the definitions of the four fields of the archived
# basepoint_one/wiener/riesz_explicit.json (maxdev, signal, amp1, triv) by trying the conventions of the neighbouring
# archived scripts: riesz_program/code/figs_riesz.py and check_large.py (ratio form, mp.dps = 20) and
# riesz_program/code/twisted.py (dev = |R - pred| / x^(1/4) with both signs of gamma, mp.dps = 25).
# Inputs beside it: riesz_explicit.json (archived), salem_riesz.json (archived, reproduced bit for bit on Hopper), g_1_400.npy.
import json, math, numpy as np, mpmath as mp
A = json.load(open('riesz_explicit.json'))
S = json.load(open('salem_riesz.json'))
g400 = np.load('g_1_400.npy'); g = g400[:60]
xs = np.array(S['riesz']['x']); R = np.array(S['riesz']['R'])
m = (xs >= 1e3) & (xs <= 1e10)
print('grid points in [1e3,1e10]:', int(m.sum()), 'first', xs[m][0], 'last', xs[m][-1])
print('g[0] =', repr(g[0]))

def rel(a, b): return abs(a - b) / abs(b)
def show(name, val, ref):
    print('%-70s %-24r rel %.3e %s' % (name, val, rel(val, ref), 'EXACT' if val == ref else ''), flush=True)

print('\n== triv ==')
for dps in [15, 17, 20, 25, 30, 50]:
    mp.mp.dps = dps
    tr = [(k, float(mp.gamma(1 + k) / (2 * mp.zeta(-2 * k, derivative=1)))) for k in range(1, 4)]
    for (k, v), (k0, v0) in zip(tr, A['triv']): show('triv k=%d gamma/(2 zeta\') dps=%d' % (k, dps), v, v0)
    tr2 = [(k, math.factorial(k) / (2 * float(mp.zeta(-2 * k, derivative=1)))) for k in range(1, 4)]
    for (k, v), (k0, v0) in zip(tr2, A['triv']): show('triv k=%d factorial/(2 float zeta\') dps=%d' % (k, dps), v, v0)

print('\n== amp1 ==')
a0 = A['amp1']
for dps in [15, 20, 25, 30, 50]:
    mp.mp.dps = dps
    for src in ['npy', 'zetazero']:
        t = g[0] if src == 'npy' else mp.zetazero(1).imag
        rho = mp.mpc(0.5, t)
        G = mp.gamma(1 - rho / 2); Z = mp.zeta(rho, derivative=1)
        c = complex(G / (2 * Z))
        cands = {
            'float(abs(G/Z))': float(abs(G / Z)),
            '2*abs(complex(G/(2Z)))': 2 * abs(c),
            'abs(2*complex(G/(2Z)))': abs(2 * c),
            'abs(complex(G/Z))': abs(complex(G / Z)),
            'float(2*abs(G/(2Z)))': float(2 * abs(G / (2 * Z))),
            '2*float(abs(G/(2Z)))': 2 * float(abs(G / (2 * Z))),
            'np.abs(np.complex128)*2': float(np.abs(np.complex128(c)) * 2),
            'float(abs(G)/abs(Z))': float(abs(G) / abs(Z)),
        }
        for k, v in cands.items(): show('amp1 %s dps=%d %s' % (k, dps, src), v, a0)

print('\n== maxdev and signal ==')
md0 = A['maxdev']; sg0 = A['signal']
ratio = R[m] / xs[m] ** 0.25
show('signal max|R/x^.25| on [1e3,1e10]', float(np.abs(ratio).max()), sg0)
show('signal max|R|/x^.25 on [1e3,1e10]', float((np.abs(R[m]) / xs[m] ** 0.25).max()), sg0)
for dps in [15, 20, 25, 30]:
    mp.mp.dps = dps
    coef = [complex(mp.gamma(1 - mp.mpc(0.5, t) / 2) / (2 * mp.zeta(mp.mpc(0.5, t), derivative=1))) for t in g]
    coefn = [complex(mp.gamma(1 - mp.mpc(0.5, -t) / 2) / (2 * mp.zeta(mp.mpc(0.5, -t), derivative=1))) for t in g]
    triv = [(k, float(mp.gamma(1 + k) / (2 * mp.zeta(-2 * k, derivative=1)))) for k in range(1, 4)]
    trivJ = A['triv']
    for tname, tv in [('trivcomp', triv), ('trivjson', trivJ)]:
        # A: figs_riesz / check_large ratio form
        predA = lambda x: sum(2 * (c * x ** (0.5j * t)).real for c, t in zip(coef, g)) + sum(r * x ** (-k) for k, r in tv) / x ** 0.25
        pA_np = np.array([predA(x) for x in xs[m]])
        pA_py = np.array([predA(float(x)) for x in xs[m]])
        for nm, p in [('A np.float64 x', pA_np), ('A python float x', pA_py)]:
            show('maxdev %s |R/x^.25-pred| dps=%d %s' % (nm, dps, tname), float(np.abs(ratio - p).max()), md0)
            show('maxdev %s |R-x^.25 pred|/x^.25 dps=%d %s' % (nm, dps, tname), float((np.abs(R[m] - xs[m] ** 0.25 * p) / xs[m] ** 0.25).max()), md0)
            show('signal %s max|pred| dps=%d %s' % (nm, dps, tname), float(np.abs(p).max()), sg0)
        # B: twisted.py style, both signs of gamma from mpmath, complex sum, dev = |R - pred| / x^.25
        co = list(zip(coef, g)) + list(zip(coefn, -g))
        trc = [(complex(r), k) for k, r in tv]
        predB = lambda x: sum(c * x ** (0.25 + 0.5j * t) for c, t in co) + sum(c * x ** (-kk) for c, kk in trc)
        pB = np.array([predB(x) for x in xs[m]])
        devB = np.abs(R[m] - pB) / xs[m] ** 0.25
        show('maxdev B twisted style dps=%d %s' % (dps, tname), float(devB.max()), md0)
        show('signal B max|pred|/x^.25 dps=%d %s' % (dps, tname), float((np.abs(pB) / xs[m] ** 0.25).max()), sg0)
        # B2: twisted style with the conjugate coefficient instead of a second mpmath evaluation
        co2 = list(zip(coef, g)) + [(c.conjugate(), -t) for c, t in zip(coef, g)]
        predB2 = lambda x: sum(c * x ** (0.25 + 0.5j * t) for c, t in co2) + sum(c * x ** (-kk) for c, kk in trc)
        pB2 = np.array([predB2(x) for x in xs[m]])
        show('maxdev B2 conj dps=%d %s' % (dps, tname), float((np.abs(R[m] - pB2) / xs[m] ** 0.25).max()), md0)
        # C: real form of R itself, 2 Re sum c x^(rho/2) + trivial, dev = |R - pred| / x^.25
        predC = lambda x: sum(2 * (c * x ** (0.25 + 0.5j * t)).real for c, t in zip(coef, g)) + sum(r * x ** (-k) for k, r in tv)
        pC = np.array([predC(x) for x in xs[m]])
        show('maxdev C |R-pred|/x^.25 dps=%d %s' % (dps, tname), float((np.abs(R[m] - pC) / xs[m] ** 0.25).max()), md0)
        show('maxdev C |R/x^.25-pred/x^.25| dps=%d %s' % (dps, tname), float(np.abs(ratio - pC / xs[m] ** 0.25).max()), md0)
        show('signal C max|pred|/x^.25 dps=%d %s' % (dps, tname), float((np.abs(pC) / xs[m] ** 0.25).max()), sg0)
    if dps == 20:
        dA = np.abs(ratio - pA_np)
        print('per point (dps 20, form A): x, R/x^.25, pred, dev')
        for x, r, p, d in zip(xs[m], ratio, pA_np, dA): print('  %.6e % .15e % .15e %.6e' % (x, r, p, d))
        print('argmax', xs[m][int(np.argmax(dA))])
# D: high precision reference (60 zeros and 3 trivial terms, all in mpmath at 40 digits)
mp.mp.dps = 40
cm = [mp.gamma(1 - mp.mpc(0.5, t) / 2) / (2 * mp.zeta(mp.mpc(0.5, t), derivative=1)) for t in g]
tm = [mp.gamma(1 + k) / (2 * mp.zeta(-2 * k, derivative=1)) for k in range(1, 4)]
pD = np.array([float((sum(2 * mp.re(c * mp.power(mp.mpf(x), mp.mpc(0.25, 0.5 * t))) for c, t in zip(cm, g)) + sum(r * mp.mpf(x) ** (-k) for k, r in zip(range(1, 4), tm))) / mp.mpf(x) ** 0.25) for x in xs[m]])
show('maxdev D mp40 reference', float(np.abs(ratio - pD).max()), md0)
show('signal D max|pred| mp40', float(np.abs(pD).max()), sg0)
