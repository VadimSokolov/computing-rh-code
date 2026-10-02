# Reconstructed diagnostic, not the authors' script. It is meant to accompany the archived code/ttest2.py,
# which produces Table tab:tilt of Chapter ch:polya (ch/ch04.tex), and uses the archived code/tilted.py
# (copied beside it) with the recovered command lines:
#   python ttest2.py 0.10 200 0.004 4 200 500 1000
#   python ttest2.py 0.02 220 0.001 4.2 1000 2000 5000
# For every row it recomputes, as ttest2.py does, err = |tomp(B/A) - dl(0.55 + i theta)| at mp.dps = 50,
# where tomp rounds the tilted value to 45 significant digits, and splits it into
#   err_true   the error of the tilted midpoint, against a 110 digit reference;
#   ref50_err  the error of the 50 digit mpmath reference;
#   conv_err   the change made by rounding the tilted value to 45 significant digits;
# and gives the error of the real part alone, which is what the column header "error of Re(xi'/xi)" names,
# the ball radii, and the time per point (median of 20 evaluations). Writes floor_check.json.
import time, json, statistics, platform, subprocess
import mpmath as mp
from flint import acb, arb, ctx
import tilted, flint

RUNS = [dict(delta=0.10, prec=200, h=0.004, xmax=4.0, heights=[200, 500, 1000]),
        dict(delta=0.02, prec=220, h=0.001, xmax=4.2, heights=[1000, 2000, 5000])]

def dl(s): return 1/s+1/(s-1)-mp.log(mp.pi)/2+mp.digamma(s/2)/2+mp.zeta(s,derivative=1)/mp.zeta(s)
def tomp(z, n): return mp.mpc(mp.mpf(z.real.mid().str(n,radius=False)), mp.mpf(z.imag.mid().str(n,radius=False)))
def e(x): return mp.nstr(x, 6)

cpu = subprocess.run("lscpu | grep 'Model name'", shell=True, capture_output=True, text=True).stdout.strip()
out = dict(host=platform.node(), cpu=cpu, flint=flint.__version__, mpmath=mp.__version__, rows=[])
for R in RUNS:
    mp.mp.dps = 50
    t0 = time.time(); S = tilted.setup(R['delta'], R['prec'], R['h'], R['xmax']); tset = time.time() - t0
    for th in R['heights']:
        A, B = tilted.xi_pair_direct(S, acb(arb('0.05'), th))
        r = B / A
        ts = []
        for k in range(20):
            t = time.perf_counter(); tilted.xi_pair_direct(S, acb(arb('0.05'), th)); ts.append(time.perf_counter() - t)
        mp.mp.dps = 50                                   # exactly as in ttest2.py
        ref50 = dl(mp.mpf('0.55') + 1j*mp.mpf(th))
        r45 = tomp(r, 45)
        err = abs(r45 - ref50)
        err_re = abs(r45.real - ref50.real)
        mp.mp.dps = 110                                  # the split
        ref = dl(mp.mpf('0.55') + 1j*mp.mpf(th))
        rmid = tomp(r, 100)
        row = dict(delta=R['delta'], prec=R['prec'], h=R['h'], xmax=R['xmax'], nodes=2*S['J']+1, theta=th,
                   err=e(err), err_re=e(err_re),
                   err_true=e(abs(rmid - ref)), err_true_re=e(abs(rmid.real - ref.real)),
                   rel_err_true_re=e(abs(rmid.real - ref.real) / abs(ref.real)),
                   ref50_err=e(abs(ref50 - ref)), conv_err=e(abs(r45 - rmid)),
                   conv_err_re=e(abs(r45.real - rmid.real)), conv_err_im=e(abs(r45.imag - rmid.imag)),
                   re_ref=mp.nstr(ref.real, 15), im_ref=mp.nstr(ref.imag, 15),
                   rad_re=e(mp.mpf(r.real.rad().str(10, radius=False))), rad_im=e(mp.mpf(r.imag.rad().str(10, radius=False))),
                   rad_A_rel=e(mp.mpf(max(A.real.rad(), A.imag.rad()).str(10, radius=False)) / mp.mpf(abs(A).mid().str(30, radius=False))),
                   ms_median=round(1000*statistics.median(ts), 2), ms_min=round(1000*min(ts), 2), setup_s=round(tset, 2))
        mp.mp.dps = 50
        out['rows'].append(row)
        print(json.dumps(row), flush=True)
json.dump(out, open('floor_check.json', 'w'), indent=1)
print('host', out['host'], out['cpu'], 'python-flint', out['flint'], 'mpmath', out['mpmath'])
