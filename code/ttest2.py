# Command lines of Table tab:tilt (Chapter ch:polya), recovered in September 2026 (misc/repro/reconstruct/ttest2/REPORT.md).
# Run beside tilted.py; the arguments are delta, the working precision in bits, the step h, the cut |x| <= xmax, and the heights:
#   python3 ttest2.py 0.10 200 0.004 4 200 500 1000        (2001 nodes: the rows with delta = 0.10)
#   python3 ttest2.py 0.02 220 0.001 4.2 1000 2000 5000    (8401 nodes: the rows with delta = 0.02)
import time, mpmath as mp, sys
from flint import acb, arb
import tilted
mp.mp.dps=50
def dl(s): return 1/s+1/(s-1)-mp.log(mp.pi)/2+mp.digamma(s/2)/2+mp.zeta(s,derivative=1)/mp.zeta(s)
def tomp(z): return mp.mpc(mp.mpf(z.real.mid().str(45,radius=False)), mp.mpf(z.imag.mid().str(45,radius=False)))
delta,prec,h,xmax = float(sys.argv[1]),int(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4])
S=tilted.setup(delta,prec,h,xmax); print('nodes',2*S['J']+1)
for th in [int(t) for t in sys.argv[5:]]:
    t=time.time(); A,B=tilted.xi_pair_direct(S,acb(arb('0.05'),th)); el=time.time()-t
    r=B/A; err=abs(tomp(r)-dl(mp.mpf('0.55')+1j*mp.mpf(th)))
    print(th, 'err', mp.nstr(err,3), 'relrad', mp.nstr(mp.mpf(r.real.rad().str(5,radius=False)),3), 'sec', round(el,3), flush=True)
