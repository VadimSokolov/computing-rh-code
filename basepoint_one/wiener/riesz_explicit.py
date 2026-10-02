# RECONSTRUCTED script, not the authors' own. It reproduces the archived file basepoint_one/wiener/riesz_explicit.json byte for
# byte (that file is byte identical to riesz_program/data/riesz_explicit.json). Reconstructed in September 2026
# (misc/repro/reconstruct/riesz_explicit/REPORT.md) and added here for the book. The archive had no script for the file:
# basepoint_one/wiener/README.md says the computation was "inline in the paper's description" and used g_1_400.npy.
# The conventions are those of the neighbouring archived scripts riesz_program/code/figs_riesz.py and check_large.py:
# mpmath at 20 digits, the first 60 ordinates of g_1_400.npy, the explicit formula for R(x)/x^(1/4) written as in those
# scripts, and the grid of salem_riesz.json (output of salem_riesz.py, Mobius sieve to N = 1e7) restricted to [1e3, 1e10],
# the range quoted in Section wr:sec:riesz of ch/strip.tex.
# Fields written:
#   maxdev = max |R(x)/x^(1/4) - E(x)| over the 36 grid points x in [1e3, 1e10], where E(x) is the explicit formula divided by x^(1/4):
#            E(x) = sum over the first 60 zeros of 2 Re[Gamma(1-rho/2)/(2 zeta'(rho)) x^(i gamma/2)] + sum_{k=1..3} triv_k x^(-k) / x^(1/4)
#   signal = max |R(x)/x^(1/4)| over the same points (attained at x = 1e3)
#   amp1   = |Gamma(1-rho_1/2)/zeta'(rho_1)|, the amplitude of the first zero's term in R(x)/x^(1/4)
#   triv   = [k, k!/(2 zeta'(-2k))] for k = 1, 2, 3
# Inputs: salem_riesz.json (written by salem_riesz.py in this folder) and g_1_400.npy, both in basepoint_one/wiener/.
# Run in basepoint_one/wiener/:  python3 riesz_explicit.py   -> riesz_explicit.json   (a few seconds; to recompute the
# input as well: python3 salem_riesz.py; python3 riesz_explicit.py)
# Reproduces: riesz_explicit.json byte for byte (182 bytes), 9 of 9 numbers identical.
import json, numpy as np, mpmath as mp
mp.mp.dps=20
S=json.load(open('salem_riesz.json'))
g=np.load('g_1_400.npy')[:60]
coef=[complex(mp.gamma(1-mp.mpc(0.5,t)/2)/(2*mp.zeta(mp.mpc(0.5,t),derivative=1))) for t in g]
triv=[(k,float(mp.gamma(1+k)/(2*mp.zeta(-2*k,derivative=1)))) for k in range(1,4)]
pred=lambda x: sum(2*(c*x**(0.5j*t)).real for c,t in zip(coef,g))+sum(r*x**(-k) for k,r in triv)/x**0.25
xs=np.array(S['riesz']['x']); R=np.array(S['riesz']['R']); m=(xs>=1e3)&(xs<=1e10)
ratio=R[m]/xs[m]**0.25
pr=np.array([pred(x) for x in xs[m]])
dev=np.abs(ratio-pr)
for x,r,p,d in zip(xs[m],ratio,pr,dev): print('x %.3e  R/x^1/4 % .9e  explicit % .9e  |diff| %.2e'%(x,r,p,d))
amp1=2*abs(coef[0])
out=dict(maxdev=float(dev.max()),signal=float(np.abs(ratio).max()),amp1=float(amp1),triv=triv)
print('grid points in [1e3,1e10]:',int(m.sum()),'  max deviation %.6e at x = %.3e'%(out['maxdev'],xs[m][int(np.argmax(dev))]),
      '  signal %.6e at x = %.3e'%(out['signal'],xs[m][int(np.argmax(np.abs(ratio)))]))
print('amp1 = |Gamma(1-rho_1/2)/zeta\'(rho_1)| =',out['amp1'],'  trivial coefficients k!/(2 zeta\'(-2k)):',triv)
json.dump(out,open('riesz_explicit.json','w'))
