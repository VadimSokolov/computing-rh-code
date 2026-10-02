# Supplementary check, not the authors' code (misc/repro/reconstruct/bosch_v7). Block 6 of
# code_bosch/verify_bosch_xi_v7.py (archive lines 183 to 199) sums the residue series of the clock density with
# 200 zeros at 90 digits; at t=0.003 about 92 digits cancel, so its value there (output line 60) is rounding error,
# as Appendix app:bosch says. Chapter ch:boschB and the appendix quote, from reviewers' recomputations with no script
# in the archive, a density of about 3.3e-72 and terms up to 2.7e20 at t=0.003. This array task computes, with the
# formula of block 6 at 150 digits, the zeros n = TASK, TASK+NTASK, ... up to NMAX and their coefficients
# c_n = Re(-2 gamma_n xi(1/2) / (i xi'(rho_n))). residue_hp_collect.py sums the series.
# Usage: python residue_hp_task.py TASK NTASK NMAX
import sys, json, time
import mpmath as mp
TASK,NTASK,NMAX=int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])
mp.mp.dps=150
def xi(z): return z*(z-1)/2*mp.pi**(-z/2)*mp.gamma(z/2)*mp.zeta(z)
half=mp.mpf(1)/2; x0=xi(half)
def c_coef(r):
    g=mp.im(r); dxi=r*(r-1)/2*mp.pi**(-r/2)*mp.gamma(r/2)*mp.zeta(r,derivative=1)
    return g, mp.re(-2*g*x0/(1j*dxi))
res=[]
for n in range(1+TASK,NMAX+1,NTASK):
    t0=time.time(); g,c=c_coef(mp.zetazero(n))
    res.append([n,mp.nstr(g,150),mp.nstr(c,150)])
    print(n,mp.nstr(g,20),mp.nstr(c,20),"%.0fs"%(time.time()-t0),flush=True)
json.dump(res,open("residue_%03d.json"%TASK,"w"))
