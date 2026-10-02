# Command line, recovered in September 2026 (misc/repro/reconstruct/ttest2/REPORT.md). Run beside tilted.py; the arguments
# are delta, the working precision in bits, the step h and the cut |x| <= xmax (the heights are fixed below):
#   python3 ttest.py 0.1 200 0.004 4    (2001 nodes, Horner evaluator xi_pair; at the heights 200 to 2000 it prints nan,
#                                        because the Horner ball of xi contains zero there)
import time, mpmath as mp, sys
from flint import acb, arb, ctx
import tilted
mp.mp.dps=40
def dl(s): return 1/s+1/(s-1)-mp.log(mp.pi)/2+mp.digamma(s/2)/2+mp.zeta(s,derivative=1)/mp.zeta(s)
def tomp(z): return mp.mpc(mp.mpf(z.real.mid().str(40,radius=False)), mp.mpf(z.imag.mid().str(40,radius=False)))
def relrad(z): return float((z.real.rad()+z.imag.rad()))/max(1e-300,float(abs(z).mid().log()) and 1)
delta,prec,h,xmax = [float(a) for a in sys.argv[1:5]]
t=time.time(); S=tilted.setup(delta,int(prec),h,xmax); print('setup',round(time.time()-t,2),'nodes',2*S['J']+1)
for th in [14,200,500,1000,2000]:
    t=time.time(); A,B=tilted.xi_pair(S,acb(0.05,th)); el=time.time()-t
    r=B/A
    err=abs(tomp(r)-dl(mp.mpf(0.55)+1j*mp.mpf(th)))
    rr = r.real.rad()
    print(th, mp.nstr(err,3), 'ball rad of ratio', mp.nstr(mp.mpf(rr.str(5,radius=False)) if False else float(rr.mid().str(5,radius=False)) if False else rr,3), round(el,4), flush=True)
