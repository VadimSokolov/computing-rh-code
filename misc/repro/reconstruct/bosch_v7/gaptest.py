# Reconstructed file, not the authors' (misc/repro/reconstruct/bosch_v7): the gaptest.py that archive lines 390
# and 392 of code_bosch/verify_bosch_xi_v7.py read. It is block 12 of that archived script (lines 344 to 388),
# copied verbatim, on the evidence that block 12 writes gaptest_output.txt. Used only to run the archived
# script unchanged, as a check that the fixed copy verify_bosch_xi_v7_fixed.py changes nothing else.
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
def h_star(removed, vmax=mp.pi, nv=60, xs=None):
    f=make_g(removed)
    xs=xs or [mp.mpf(k)/4 for k in range(0,49)]
    last_ok=0
    for iv in range(1,nv+1):
        v=vmax*iv/nv*mp.mpf('0.999')
        bad=False
        for x in xs:
            g,gp=f(mp.mpc(x,v))
            if abs(g)<mp.mpf('1e-25') or mp.im(-gp/g)< -mp.mpf('1e-18'):
                bad=True; break
        if bad: return last_ok, v
        last_ok=v
    return last_ok, None
cases=[("freqs 1,2,3,... (logistic)",[],1),("freqs 1,3,4,5,...",[2],2),("freqs 1,4,5,6,...",[2,3],3),
       ("freqs 2,3,4,...",[1],1),("freqs 1,2,5,6,... (interior gap 3)",[3,4],1),("freqs 1,2,3,6,7,... (interior gap 3)",[4,5],1)]
for name,rem,gap in cases:
    t0=time.time()
    ok,bad=h_star(rem)
    P("%s: first gap %d, predicted strip pi/gap = %s; Pick holds up to v=%s, first failure at v=%s (%.0fs)"%(name,gap,mp.nstr(mp.pi/gap,5),mp.nstr(ok,5),mp.nstr(bad,5) if bad else "none below pi",time.time()-t0))
open("gaptest_output.txt","w").write("\n".join(out)+"\n")
