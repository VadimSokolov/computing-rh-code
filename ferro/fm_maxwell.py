# Chapter 13, the Maxwell series of the Polya law, in the book's units.
# M_k are independent symmetric Maxwell variables: density x^2 e^{-x^2/2}/sqrt(2 pi), E e^{i theta M} = (1 - theta^2) e^{-theta^2/2}.
# When every gamma_k is real,
#   E exp(i theta sum_k M_k/gamma_k) = prod_k (1 - theta^2/gamma_k^2) e^{-theta^2/(2 gamma_k^2)} = (Xi(theta)/Xi(0)) e^{-Var X theta^2/4},
# so S = sum_k M_k/gamma_k has the law of X + G, with G normal of variance Var X/2 independent of X.
# The first 1700 terms are sampled and the rest replaced by a normal of variance 3 sum_{k>1700} gamma_k^-2. The empirical law
# of S is compared with p * N(0, Var X/2) (Kolmogorov distance), and kappa_4(S) with kappa_4(X) = -12 sum gamma_k^-4.
# Usage: python3 fm_maxwell.py [n]   Reads zeros_1700.txt; writes fm_maxwell.json.
import json, sys, time
import numpy as np, mpmath as mp
from scipy.special import ndtr
from scipy.stats import kstwobign

t0 = time.time()
n = int(sys.argv[1]) if len(sys.argv) > 1 else 1000000
mp.mp.dps = 30
xi = lambda s: s*(s - 1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
var = float(mp.diff(xi, mp.mpf(1)/2, 2)/xi(mp.mpf(1)/2)); X0 = float(xi(mp.mpf(1)/2))
gam = np.array([float(l) for l in open('zeros_1700.txt') if l.strip()])
G = gam[-1]
tail2 = var/2 - np.sum(gam**-2.)                      # sum_{k>1700} gamma_k^-2, from Var X = 2 sum gamma_k^-2
p2 = np.sum(gam**-4.) + (3*np.log(G/(2*np.pi)) + 1)/(18*np.pi*G**3)    # sum gamma^-4 with the smooth tail
rng = np.random.default_rng(20261001)
S = np.empty(n); chunk = 20000; w = 1/gam
for i in range(0, n, chunk):
    m = min(chunk, n - i)
    M = np.sqrt(rng.standard_gamma(1.5, size=(m, len(gam)))*2)*(2*rng.integers(0, 2, size=(m, len(gam))) - 1)
    S[i:i + m] = M @ w + np.sqrt(3*tail2)*rng.standard_normal(m)

pi = np.pi
def Phi(t):
    t = np.abs(t); s = 0.0
    for k in range(1, 8):
        s = s + 2*(2*pi**2*k**4*np.exp(4.5*t) - 3*pi*k**2*np.exp(2.5*t))*np.exp(-pi*k**2*np.exp(2*t))
    return s
t = np.linspace(-2.4, 2.4, 4801); dt = t[1] - t[0]; pt = Phi(t)/X0; pt /= np.sum(pt)*dt
sg = np.sqrt(var/2)
s_grid = np.linspace(-2.0, 2.0, 4001)
F = np.array([np.sum(pt*ndtr((s - t)/sg))*dt for s in s_grid])
Ss = np.sort(S); Fs = np.interp(Ss, s_grid, F)
i = np.arange(1, n + 1)
D = float(max(np.max(i/n - Fs), np.max(Fs - (i - 1)/n)))
c = S - S.mean(); m2 = np.mean(c**2); m4 = np.mean(c**4)
k4 = m4 - 3*m2**2
se_k4 = float(np.sqrt(np.var(c**4 - 6*m2*c**2)/n))     # delta method, leading term
R = dict(n=n, var_X=var, var_S=float(m2), var_S_exact=1.5*var, tail_var=float(3*tail2),
         KS=D, sqrt_n_KS=D*np.sqrt(n), KS_pvalue=float(kstwobign.sf(D*np.sqrt(n))),
         kappa4_S=float(k4), kappa4_se=se_k4, kappa4_exact=float(-12*p2), seconds=time.time() - t0)
print(json.dumps(R, indent=1))
json.dump(R, open('fm_maxwell.json', 'w'), indent=1)
