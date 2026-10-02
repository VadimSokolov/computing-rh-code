# Supplementary check, not the authors' code (misc/repro/reconstruct/bosch_v7): collects the array tasks of
# residue_hp_task.py and sums the residue series of block 6 of code_bosch/verify_bosch_xi_v7.py at 150 digits,
# with N = 200 (as in block 6), 220 and 250 zeros, at the three times of block 6 (t = 0.02, 0.005, 0.003).
# It reports the sum, the largest term, the number of digits that cancel (log10 of largest term over the sum), and
# the size of the last term (the truncation). The ordinates are checked against Arb (python-flint) at 520 bits.
import json, glob
import mpmath as mp
mp.mp.dps=150
rows=sorted(r for f in glob.glob("residue_*.json") for r in json.load(open(f)))
ns=[r[0] for r in rows]
assert ns==list(range(1,len(ns)+1)), "missing zeros"
gam=[mp.mpf(r[1]) for r in rows]; c=[mp.mpf(r[2]) for r in rows]
res={"nzeros":len(ns),"dps":150}
try:
    from flint import acb, ctx
    ctx.prec=520
    zs=acb.zeta_zeros(1,len(ns))
    d=max(abs(mp.mpf(z.imag.str(150,radius=False))-g) for z,g in zip(zs,gam))
    res["max |gamma_mpmath - gamma_arb|"]=mp.nstr(d,5)
except Exception as e:
    res["arb check"]="failed: %r"%e
res["c_1..c_4"]=[mp.nstr(x,12) for x in c[:4]]
res["signs alternate"]=all(c[i]*c[i+1]<0 for i in range(len(c)-1))
res["series"]={}
for ts in ["0.02","0.005","0.003"]:
    t=mp.mpf(ts)
    for N in [200,220,250]:
        if N>len(c): continue
        terms=[c[i]*mp.e**(-gam[i]**2*t) for i in range(N)]
        S=mp.fsum(terms)
        imax=max(range(N),key=lambda i:abs(terms[i]))
        big=abs(terms[imax])
        res["series"]["t=%s N=%d"%(ts,N)]=dict(sum=mp.nstr(S,15),max_term=mp.nstr(big,6),max_term_n=imax+1,
            digits_cancelled=mp.nstr(mp.log10(big/abs(S)),4),last_term=mp.nstr(terms[-1],4))
json.dump(res,open("residue_hp.json","w"),indent=1)
print(json.dumps(res,indent=1))
