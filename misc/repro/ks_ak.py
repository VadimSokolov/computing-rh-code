# Recomputes the arithmetic factors a_1, a_2, a_3 that riesz_program/code/ks.py hard-codes
# (its function ak is defined but never called). Same Euler product over the primes below
# 2*10^5 with each local factor written through 2F1(k,k;1;1/p), as the Riesz note says.
import json
import mpmath as mp
import sympy

mp.mp.dps = 30
P = 200000
HARD = {1: 1.0000000000000107, 2: 0.607927332969108, 3: 0.04932184233405926}


def ak(k, P=P):
    s = mp.mpf(0)
    for p in sympy.primerange(2, P):
        p = mp.mpf(p)
        s += k * k * mp.log(1 - 1 / p) + mp.log(mp.hyp2f1(k, k, 1, 1 / p))
    return mp.e ** s


out = {}
for k in (1, 2, 3):
    a = ak(k)
    out[k] = {'a_k': float(a), 'hard_coded_in_ks_py': HARD[k], 'rel_diff': float(abs(a - HARD[k]) / a)}
    print(k, mp.nstr(a, 15), HARD[k], out[k]['rel_diff'], flush=True)
exact2 = 6 / mp.pi ** 2
out['6/pi^2'] = float(exact2)
out['a_2 truncated / (6/pi^2) - 1'] = float(ak(2) / exact2 - 1)
print('6/pi^2', mp.nstr(exact2, 15), 'ratio-1', out['a_2 truncated / (6/pi^2) - 1'])
G = lambda k: mp.barnesg(1 + k) ** 2 / mp.barnesg(1 + 2 * k)
out['leading'] = {k: float(out[k]['a_k'] * G(k)) for k in (1, 2, 3)}
out['42 a_3 / 9!'] = float(42 * out[3]['a_k'] / mp.factorial(9))
print(out['leading'], out['42 a_3 / 9!'])
json.dump(out, open('ks_ak.json', 'w'), indent=1)
