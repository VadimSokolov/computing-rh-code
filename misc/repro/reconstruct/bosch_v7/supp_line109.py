# Reconstruction, not the authors' code (misc/repro/reconstruct/bosch_v7). Output line 109 of
# code_bosch/verify_bosch_xi_v7_output.txt ("first gap ... largest later gap among n<1000 ...") has no generating
# code in the archive (Appendix app:bosch says so); Chapter ch:boschB quotes it before Question bs:q:gaprefined.
# Ordinates by mpmath.zetazero at 20 digits, as in blocks 3 and 5 of the archived script. The printed line uses the
# zeros with n < 1000; the extra lines below (not in the archive) repeat it for the first 1000 zeros, which is what
# the chapter says, and give the height of zero 1000 that the chapter quotes.
import mpmath as mp
from multiprocessing import Pool
def z(n):
    mp.mp.dps=20
    return mp.im(mp.zetazero(n))
if __name__=='__main__':
    with Pool(4) as pool:
        gam=pool.map(z,range(1,1001),chunksize=5)
    mp.mp.dps=20
    out=[]
    def P(*a):
        s=" ".join(str(x) for x in a); print(s,flush=True); out.append(s)
    g=gam[:999]  # zeros with n<1000
    gaps=[g[i+1]-g[i] for i in range(len(g)-1)]
    i=max(range(1,len(gaps)),key=lambda k:gaps[k])
    P("first gap",mp.nstr(gaps[0],8)," largest later gap among n<1000:",mp.nstr(gaps[i],6),"at n=",i+1)
    open("line109_output.txt","w").write("\n".join(out)+"\n")
    gaps=[gam[i+1]-gam[i] for i in range(len(gam)-1)]
    i=max(range(1,len(gaps)),key=lambda k:gaps[k])
    later=sorted(range(1,len(gaps)),key=lambda k:-gaps[k])[:3]
    print("(extra) first 1000 zeros: first gap %s, largest later gap %s at n=%d, next ones %s; gamma_1000 = %s"%(mp.nstr(gaps[0],10),mp.nstr(gaps[i],8),i+1,[(k+1,mp.nstr(gaps[k],6)) for k in later[1:]],mp.nstr(gam[-1],12)),flush=True)
