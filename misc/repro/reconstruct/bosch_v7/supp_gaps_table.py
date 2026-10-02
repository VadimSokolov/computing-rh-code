# Reconstruction, not the authors' code (misc/repro/reconstruct/bosch_v7). Block 13 of code_bosch/verify_bosch_xi_v7.py
# (archive lines 393 to 413: fails, threshold, freqs) with the P and make_g of block 12 (archive lines 347 to 367),
# all copied verbatim at 30 digits; only the case list (archive line 414) and the output file name differ.
# Cases: the ten of the stale output lines 93 to 102 (an earlier version of block 13 that is not in the archive), in
# their printed order, then the empty set, so that together with block 13 every computed row of Table bs:tab:gaps
# in Chapter ch:boschB is made by the printed code.
import mpmath as mp, time
mp.mp.dps=30
out=[]
def P(*a):
    s=" ".join(str(x) for x in a); print(s,flush=True); out.append(s)
# L(y) = (1/4) sech^2(y/2): density of log W for freqs 1,2,3,... (char fn pi t/sinh pi t)
# removing frequency k multiplies chi by (1+t^2/k^2): g = prod_k (1 - D^2/k^2) L
# Taylor expansion approach: represent L^(j) via mp.diff
def make_g(removed):
    # polynomial in D: prod (1 - D^2/k^2) -> coefficients c_j of D^j
    poly=[mp.mpf(1)]
    for k in removed:
        new=[mp.mpf(0)]*(len(poly)+2)
        for j,c in enumerate(poly):
            new[j]+=c; new[j+2]-=c/mp.mpf(k)**2
        poly=new
    L=lambda y: mp.sech(y/2)**2/4
    def g_and_gp(y):
        ders=mp.diffs(L,y,len(poly))
        ders=list(ders)
        g=sum(c*ders[j] for j,c in enumerate(poly))
        gp=sum(c*ders[j+1] for j,c in enumerate(poly))
        return g,gp
    return g_and_gp
def fails(f,v,xs):
    for x in xs:
        g,gp=f(mp.mpc(x,v))
        if abs(g)<mp.mpf('1e-25') or mp.im(-gp/g)< -mp.mpf('1e-20'): return True
    return False
def threshold(removed,vmax=mp.pi):
    f=make_g(removed); xs=[mp.mpf(k)/5 for k in range(0,81)]
    lo=mp.mpf('0.01'); hi=vmax*mp.mpf('0.9995')
    # coarse scan then bisection
    grid=[hi*k/40 for k in range(1,41)]
    found=False
    for v in grid:
        if fails(f,v,xs): hi=v; found=True; break
        lo=v
    if not found: return None
    for _ in range(14):
        m=(lo+hi)/2
        if fails(f,m,xs): hi=m
        else: lo=m
    return (lo+hi)/2
def freqs(rem,n=12): return [k for k in range(1,n+1) if k not in rem][:7]
cases=[[2],[2,3],[2,4],[2,3,6],[2,3,5,6],[2,3,4],[3,4],[4,5],[2,5],[3],[]]
for rem in cases:
    fr=freqs(rem); gaps=[fr[i+1]-fr[i] for i in range(len(fr)-1)]
    t0=time.time(); th=threshold(rem)
    thv = th if th is not None else mp.pi
    P("removed %s: freqs %s..., gaps %s...: first gap %d, max gap %d; threshold strip %s; pi/threshold = %s; pi/first gap = %s (%.0fs)"%(rem,fr,gaps,gaps[0],max(gaps),mp.nstr(thv,6),mp.nstr(mp.pi/thv,6),mp.nstr(mp.pi/gaps[0],5),time.time()-t0))
open("gaps_table_output.txt","w").write("\n".join(out)+"\n")
