# Comparison for the bosch_v7 gap (misc/repro/reconstruct/bosch_v7); not the authors' code. Run on Hopper in the
# agent folder: python src/compare_bosch_v7.py. It reads the archived output (src/verify_bosch_xi_v7_output.txt,
# a copy of code_bosch/verify_bosch_xi_v7_output.txt), the rerun of the fixed copy (fixed/), the rerun of the
# unchanged archived script beside the reconstructed gaptest.py (equiv/), the supplementary reconstructions
# (supp37/, supp109/, suppgaps/), and the 150 digit checks (residue/, talbot/). It writes compare_bosch_v7.json.
# Part 1 compares the archived output with the rerun line by line and number by number (relative tolerance 1e-9;
# timings, printed as "(Ns)", are listed apart). Part 2 compares the fixed copy with the unchanged archived script.
# Part 3 checks every number that Chapters ch:boschA, ch:boschB and Appendix app:bosch quote from this output, and
# a few that follow from it by arithmetic or by the script's own functions; each check names its source line.
import re, json, glob, os
import mpmath as mp
import numpy as np
mp.mp.dps=40
def lines(path): return open(path).read().rstrip("\n").split("\n")
arch=lines("src/verify_bosch_xi_v7_output.txt")
fixed=lines("fixed/joined_output.txt"); equiv=lines("equiv/joined_output.txt")
s37=lines("supp37/lines37_39_output.txt"); s109=lines("supp109/line109_output.txt"); sgap=lines("suppgaps/gaps_table_output.txt")
assert (len(arch),len(fixed),len(s37),len(s109),len(sgap))==(109,95,3,1,11),(len(arch),len(fixed),len(s37),len(s109),len(sgap))
# the rerun counterpart of the 109 archived lines
full=fixed[:36]+s37+fixed[36:89]+sgap[:10]+fixed[89:95]+s109
prov=(["fixed script, blocks 1 to 3"]*36+["reconstruction supp_lines37_39.py"]*3+["fixed script, blocks 4 to 12"]*53
      +["printed block 13 code on the ten stale cases (supp_gaps_table.py)"]*10+["fixed script, block 13"]*6+["reconstruction supp_line109.py"])
NUM=re.compile(r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?')
TIM=re.compile(r'\((\d+)s\)')
def split(line):
    t=TIM.findall(line); core=TIM.sub('(#s)',line)
    return NUM.findall(core),NUM.sub('#',core),t
def mpf(s): return mp.mpf(s.rstrip('.'))
def rel(a,b):
    A,B=mpf(a),mpf(b)
    if A==B: return mp.mpf(0)
    return abs(A-B)/max(abs(A),abs(B))
# ---------------- Part 1 ----------------
per_line=[]; groups={}
for i,(a,r) in enumerate(zip(arch,full)):
    na,ska,ta=split(a); nr,skr,tr=split(r)
    d=dict(line=i+1,source=prov[i],identical=(a==r),timings_archive=ta,timings_rerun=tr,text_same=(ska==skr),n_numbers=len(na))
    if len(na)!=len(nr): d["number_count"]=(len(na),len(nr)); diffs=[(k,x,'') for k,x in enumerate(na)]; mr=None
    else:
        rs=[rel(x,y) for x,y in zip(na,nr)]; mr=max(rs) if rs else mp.mpf(0)
        diffs=[(k,x,y,mp.nstr(q,3)) for k,(x,y,q) in enumerate(zip(na,nr,rs)) if q>mp.mpf('1e-9')]
    d["max_rel"]=None if mr is None else mp.nstr(mr,3); d["differences"]=diffs
    if not d["identical"]: d["archive"]=a; d["rerun"]=r
    per_line.append(d)
    g=groups.setdefault(prov[i],dict(lines=0,identical_lines=0,identical_but_timing=0,numbers=0,numbers_agree=0,numbers_differ=0,max_rel_diff=mp.mpf(0),timing_pairs=[]))
    g["lines"]+=1; g["numbers"]+=len(na)
    g["identical_lines"]+=int(a==r); g["identical_but_timing"]+=int(a!=r and ska==skr and not diffs)
    g["numbers_differ"]+=len(diffs); g["numbers_agree"]+=len(na)-len(diffs)
    if mr is not None: g["max_rel_diff"]=max(g["max_rel_diff"],mr)
    if ta or tr: g["timing_pairs"].append((i+1,ta,tr))
for g in groups.values(): g["max_rel_diff"]=mp.nstr(g["max_rel_diff"],3)
# ---------------- Part 2 ----------------
p2=dict(lines_fixed=len(fixed),lines_equiv=len(equiv),identical=sum(x==y for x,y in zip(fixed,equiv)),
        differ_only_in_timing=sum(x!=y and split(x)[:2]==split(y)[:2] for x,y in zip(fixed,equiv)),
        other_differences=[(i+1,x,y) for i,(x,y) in enumerate(zip(fixed,equiv)) if split(x)[:2]!=split(y)[:2]])
# per block output file line counts of the fixed run, against Table tab:app:bosch
files=["verify_v2_output.txt","verify_v3_output.txt","thorin_p_output.txt","open_problems_output.txt","op5_output.txt","gW_output.txt","gpos_output.txt","gaptest_output.txt","gaptest3_output.txt"]
counts={f:len(lines("fixed/"+f)) for f in files}
src=lines("src/verify_bosch_xi_v7.py")
starts=[i+1 for i,l in enumerate(src) if l.startswith("import mpmath") or l.startswith("import numpy")]
# ---------------- Part 3 ----------------
checks=[]
def L(n): return full[n-1]
def grab(n,pat,k=1):
    m=re.search(pat,L(n)); assert m,(n,pat,L(n)); return m.group(k)
def quantum(q):
    m=re.fullmatch(r'([-+]?)(\d*)(?:\.(\d*))?(?:e([-+]?\d+))?',q); assert m,q
    return mp.mpf(10)**(int(m.group(4) or 0)-len(m.group(3) or ''))
def rounds(where,what,quoted,value,scale=1,src=""):
    v=mp.mpf(value)/scale; q=mp.mpf(quoted)
    ok=abs(v-q)<=quantum(quoted)/2*(1+mp.mpf('1e-9'))
    checks.append(dict(where=where,what=what,quoted=quoted+("" if scale==1 else " x %s"%mp.nstr(scale,3)),value=mp.nstr(mp.mpf(value),12),source=src,ok=bool(ok)))
def truth(where,what,quoted,cond,value,src=""):
    checks.append(dict(where=where,what=what,quoted=quoted,value=value,source=src,ok=bool(cond)))
MIN=r"min Im\(-x h'/h\)=([-0-9.e+]+)"
# boschB, Section bs:sec:bosch and after Corollary bs:cor:uncond
for n,aq in [(11,"(0.4,0.6)"),(12,"(0.4,0.55)"),(14,"(0.5,0.5)"),(17,"(0.3,0.7)")]:
    v=mp.mpf(grab(n,MIN)); truth("boschB:72","min Im eps_h positive for (alpha,q)="+aq,"> 0",v>0,mp.nstr(v,4),"line %d"%n)
rounds("boschB:72","min Im eps_h for (0.6,0.3)","-0.0114",grab(16,MIN),src="line 16")
rounds("boschB:72","min Im eps_h for (0.5,0.6)","-191",grab(15,MIN),src="line 15")
rounds("boschB:72","min Im eps_h for (0.4,0.75)","-1.1e16",grab(13,MIN),src="line 13")
truth("boschB:72, app:bosch","grid of 49 x 199 points (radii 10^(k/4), k=-24..24; angles j pi/200)","49 x 199",src.count("rs=[mp.mpf(10)**(k/4) for k in range(-24,25)]")==1 and src.count("ths=[mp.pi*j/200 for j in range(1,200)]")==1,"49 x 199","script lines 56, 57")
rounds("boschB:149, app:bosch","ray scan: number of points","885",grab(10,r"over (\d+) points"),src="line 10")
rounds("boschB:149","ray scan: minimum of Im eps_phi(r e^{i theta})/r","3.8e-4",grab(10,r"points: ([-0-9.e+]+) at"),src="line 10")
rounds("boschB:149, app:bosch","toy p_* = 1-0.2/pi","0.9363",grab(18,r"= ([0-9.]+)$"),src="line 18")
for n,pp in [(19,"0.80"),(20,"0.90"),(21,"0.93")]:
    v=mp.mpf(grab(n,MIN)); truth("boschB:149","toy p=%s: no poles in Omega, scan positive"%pp,"False, > 0","poles in slit plane=False" in L(n) and v>0,mp.nstr(v,4),"line %d"%n)
for n,pp in [(22,"0.94"),(23,"0.97")]:
    v=mp.mpf(grab(n,MIN)); truth("boschB:149, app:bosch","toy p=%s: poles in Omega, scan stays positive"%pp,"True, > 0","poles in slit plane=True" in L(n) and v>0,mp.nstr(v,4),"line %d"%n)
rounds("boschB:149","toy p=1: minimum","-2.08",grab(24,MIN),src="line 24")
rounds("boschB:149","toy p=1: radius of the minimum","178",grab(24,r"at \('([0-9.e+]+)'"),src="line 24")
rounds("boschB:149","toy p=1: angle/pi of the minimum","0.965",grab(24,r"', '([0-9.]+)'\)"),src="line 24")
truth("app:bosch","sign test: no violation at p=0.95","[None, None, None, None]",L(8).endswith("[None, None, None, None]"),L(8).split(":")[-1].strip(),"line 8")
truth("app:bosch","sign test: violation at p=1","first entry 29",L(9).endswith("[29, None, None, None]"),L(9).split(":")[-1].strip(),"line 9")
# toy pole positions against the grid of block 2 (arithmetic, not in the output)
for pp in ["0.94","0.97"]:
    p=mp.mpf(pp); ang=(mp.pi-mp.mpf('0.2'))/p/mp.pi; rad=mp.mpf(196)**(1/p)
    lastray=mp.mpf(199)/200; kr=mp.floor(4*mp.log10(rad)); r1,r2=mp.mpf(10)**(kr/4),mp.mpf(10)**((kr+1)/4)
    between=r1<rad<r2 and rad not in (r1,r2)
    truth("app:bosch","toy p=%s: pole angle/pi %s, radius %s; beyond the last ray 0.995 or between grid radii %s and %s"%(pp,mp.nstr(ang,5),mp.nstr(rad,5),mp.nstr(r1,4),mp.nstr(r2,4)),
          "between grid points or beyond the last ray",ang>lastray or between,"angle beyond last ray" if ang>lastray else "radius between grid radii","arithmetic from script lines 57, 78")
# boschB after Theorem bs:thm:ggc: block 3 and the reconstructed lines 37 to 39
for n in [33,34,35,36]:
    a,b=re.findall(r"= ?(-?[0-9.]+)",L(n))[-2:]
    truth("boschB:172","log phi(s^p) reproduced to all ten printed digits (line %d)"%n,"equal",a==b,a+" | "+b,"line %d"%n)
rounds("boschB:172","log phi(s^p) at p=0.9, s=10","-0.1823802136",grab(36,r"log phi\(s\^p\)=(-?[0-9.]+)"),src="line 36")
rounds("boschB:172","p=0.5, tau=100: zero sum, 30 zeros","2.74",grab(26,r"30 zeros = ([0-9.e-]+)"),mp.mpf('1e-4'),"line 26")
rounds("boschB:172","p=0.5, tau=100: zero sum, 100 zeros","3.18",grab(37,r"100 zeros ([0-9.e-]+)"),mp.mpf('1e-4'),"line 37 (reconstructed)")
rounds("boschB:172","p=0.5, tau=100: zero sum, 200 zeros","3.35",grab(37,r"200 zeros ([0-9.e-]+)"),mp.mpf('1e-4'),"line 37 (reconstructed)")
rounds("boschB:172","p=0.5, tau=100: exact u_p","3.67",grab(26,r"xi'/xi = ([0-9.e-]+)"),mp.mpf('1e-4'),"line 26")
rounds("boschB:172","p=0.9, tau=gamma_1^(2/p): zero sum, 30 zeros","3.504",grab(31,r"30 zeros = ([0-9.e-]+)"),mp.mpf('1e-3'),"line 31")
rounds("boschB:172","p=0.9, tau=gamma_1^(2/p): zero sum, 100 zeros","3.644",grab(38,r"100 zeros ([0-9.e-]+)"),mp.mpf('1e-3'),"line 38 (reconstructed)")
rounds("boschB:172","p=0.9, tau=gamma_1^(2/p): zero sum, 200 zeros","3.696",grab(38,r"200 zeros ([0-9.e-]+)"),mp.mpf('1e-3'),"line 38 (reconstructed)")
rounds("boschB:172","p=0.9, tau=gamma_1^(2/p): exact u_p","3.797",grab(31,r"xi'/xi = ([0-9.e-]+)"),mp.mpf('1e-3'),"line 31")
v=mp.mpf(grab(39,r"min u_p = ([-0-9.e+]+)")); truth("boschB:172, app:bosch","u_p positive on 148 points","148 points, min > 0",grab(39,r"grid: (\d+) points")=="148" and v>0,"148 points, min %s"%mp.nstr(v,4),"line 39 (reconstructed)")
# boschB Theorem bs:thm:eta
for n in [40,41,42,43]:
    a,b=re.findall(r"([0-9.]+(?:e-?\d+)?)",L(n).split("density")[1])[0],L(n).split()[-1]
    truth("boschB:306","series and eta form agree to fifteen digits (line %d)"%n,"equal",a==b,a+" | "+b,"line %d"%n)
rounds("boschB:306","int e^{-2t} f_C","0.214383752348",grab(44,r"f = ([0-9.]+)"),src="line 44")
rounds("boschB:306","sech(pi/sqrt 2)","0.214383752348",grab(44,r"= ([0-9.]+)$"),src="line 44")
# boschB Table bs:tab:prime and the remark after it
TE={"0.5":"0.707","0.9":"3.20","0.99":"31.8","0.999":"318","0.9999":"3183"}
T1={"0.5":"4.07","0.9":"11.1","0.99":"69.6","0.999":"553","0.9999":"4818"}
RA={"0.5":"5.75","0.9":"3.47","0.99":"2.19","0.999":"1.74","0.9999":"1.51"}
TA={"0.5":"274","0.9":"210","0.99":"5.29e3","0.999":"3.10e5","0.9999":"2.33e7"}
for n,pp in zip(range(61,66),["0.5","0.9","0.99","0.999","0.9999"]):
    t1=grab(n,r"\|z\|>([0-9.e+]+)"); ta=grab(n,r"tau_1=([0-9.e+]+)"); te=grab(n,r"strip exit ([0-9.e+]+)")
    rounds("boschB:273-277","Table bs:tab:prime, p=%s: strip exit t_e"%pp,TE[pp],te,src="line %d"%n)
    rounds("boschB:273-277","Table bs:tab:prime, p=%s: t_1"%pp,T1[pp],t1,src="line %d"%n)
    rounds("boschB:273-277","Table bs:tab:prime, p=%s: t_1/t_e"%pp,RA[pp],mpf(t1)/mpf(te),src="line %d"%n)
    rounds("boschB:273-277","Table bs:tab:prime, p=%s: tau_1"%pp,TA[pp],ta,src="line %d"%n)
truth("boschB:265","-zeta'/zeta(sigma) < 1/(sigma-1) at 80 points, 10^-5.9 <= sigma-1 <= 10^2","True",L(66).endswith("True") and any("range(-20,60)" in l and "10)**(-k/10)" in l for l in src),"True","line 66, script line 217")
# remark after Table bs:tab:prime: predicted ratio and the sign change of m_p at p=0.9999 (block 7's L at 25 digits)
with mp.workdps(25):
    def mzl(sig): return -mp.zeta(sig,derivative=1)/mp.zeta(sig)
    def Lp(t,th):
        z=t*mp.expj(th); w=mp.mpf(1)/2+z; sw=mp.re(w)
        main=mp.im(z*(1/(2*w)+1/(w-1)+mp.log(w/(2*mp.pi))/2))
        return main - t/(6*sw**2) - t*mzl(sw)
    p=mp.mpf('0.9999'); th=p*mp.pi/2; c=mp.cos(th); te=1/(2*c)
    tmin=(mp.mpf(1)/2)/c*mp.mpf('1.0001')
    ts=[tmin*mp.mpf(10)**(k/50) for k in range(0,50*14)]
    vals=[Lp(t,th) for t in ts]; kb=max(k for k,vv in enumerate(vals) if vv<=0)
    lo,hi=ts[kb],ts[kb+1]
    for _ in range(60):
        m=(lo+hi)/2
        if Lp(m,th)<=0: lo=m
        else: hi=m
    rounds("boschB:285","p=0.9999: leading term ratio 1+4/log(t_e/2 pi)","1.64",1+4/mp.log(te/(2*mp.pi)),src="arithmetic")
    rounds("boschB:285","p=0.9999: sign change of m_p, as a ratio to t_e","1.53",lo/te,src="recomputed with block 7's L, bisection between grid points")
    rounds("boschB:285","p=0.9999: last grid point with m_p <= 0, as a ratio to t_e","1.51",ts[kb]/te,src="recomputed with block 7's L")
# boschB after Theorem bs:thm:nb
rounds("boschB:341","c_1","10163.5",grab(50,r"\['([0-9.]+)'"),src="line 50")
rounds("boschB:341","c_1/gamma_1^2","50.87",grab(50,r"gamma_1\^2 = ([0-9.]+)"),src="line 50")
rounds("boschA:128, boschB:341","A_1","50.9",grab(50,r"gamma_1\^2 = ([0-9.]+)"),src="line 50")
truth("boschB:341","c_n alternate in sign for n <= 200","True",L(51).endswith("True"),"True","line 51")
rounds("boschB:341","smallest gap gamma_{n+1}^2-gamma_n^2","159",grab(52,r"n<200: ([0-9.]+)"),src="line 52")
truth("boschB:341","smallest gap at n=4","n=4",grab(52,r"at n=(\d+)")=="4","n="+grab(52,r"at n=(\d+)"),"line 52")
LR=re.findall(r"'([0-9.e-]+)'",L(53))
for q,v,nn in zip(["0.0462","0.0144","0.00524","0.00322","0.00195"],LR,[1,10,50,100,200]):
    rounds("boschB:341","log|c_n|/gamma_n^2 at n=%d"%nn,q,v,src="line 53")
PI4=re.findall(r"'([0-9.e-]+)'",L(54))
truth("boschB:341","log|c_n|/gamma_n^2 tracks pi/(4 gamma_n)","ratio near 1",all(abs(mpf(a)/mpf(b)-1)<mp.mpf('0.2') for a,b in zip(LR,PI4)),", ".join(mp.nstr(mpf(a)/mpf(b),3) for a,b in zip(LR,PI4)),"lines 53, 54")
a,b=grab(58,r"90 digits\) ([0-9.e+-]+)"),grab(58,r"Talbot ([0-9.e+-]+)")
rounds("boschB:341","t=0.02: residue series","80.10419773",a,src="line 58"); rounds("boschB:341","t=0.02: Talbot","80.10419773",b,src="line 58")
a,b=grab(59,r"90 digits\) ([0-9.e+-]+)"),grab(59,r"Talbot ([0-9.e+-]+)")
rounds("boschB:341","t=0.005: residue series","6.537917064e-24",a,src="line 59")
truth("boschB:341","t=0.005: series and Talbot agree to nine digits","rel diff < 1e-9 but not to ten digits",rel(a,b)<mp.mpf('1e-9') and a!=b,"rel diff "+mp.nstr(rel(a,b),3),"line 59")
for n in [56,60]:
    v=mp.mpf(grab(n,r"Talbot[A-Za-z ]*? (-?[0-9.e+-]+)")); truth("app:bosch","t=0.003: the Talbot value is negative (line %d)"%n,"< 0",v<0,mp.nstr(v,5),"line %d"%n)
# 150 digit residue series and Talbot (supplementary)
R=json.load(open("residue/residue_hp.json")) if os.path.exists("residue/residue_hp.json") else None
talbot={json.load(open(f))["t"]:json.load(open(f)) for f in glob.glob("talbot/talbot_*.json")}
if R:
    s=R["series"]; key=max((k for k in s if k.startswith("t=0.003")),key=lambda k:int(k.split("N=")[1]))
    rounds("boschB:341, app:bosch","t=0.003: density (150 digits, %s)"%key.split()[1],"3.3e-72",s[key]["sum"],src="residue_hp.json")
    rounds("boschB:341, app:bosch","t=0.003: largest term","2.7e20",s[key]["max_term"],src="residue_hp.json")
    rounds("boschB:341, app:bosch","t=0.003: digits that cancel","92",s[key]["digits_cancelled"],src="residue_hp.json")
    truth("boschB:341, app:bosch","t=0.003: more digits cancel than the 90 of line 60 and the 25 of line 56","> 90",mp.mpf(s[key]["digits_cancelled"])>90,s[key]["digits_cancelled"],"residue_hp.json")
    for t,q,n in [("0.02","80.10419773",58),("0.005","6.537917064e-24",59)]:
        k2="t=%s N=200"%t; rounds("boschB:341","t=%s: 150 digit series, 200 zeros, against line %d"%(t,n),q,s[k2]["sum"],src="residue_hp.json")
    for t,v in talbot.items():
        k2="t=%s N=%d"%(t,max(int(k.split("N=")[1]) for k in s if k.startswith("t=%s "%t)))
        truth("boschB:341","t=%s: Talbot at 150 digits against the 150 digit series (%s)"%(t,k2.split()[1]),"agree",rel(v["talbot"],s[k2]["sum"])<mp.mpf('1e-8'),"%s vs %s, rel diff %s"%(v["talbot"][:14],s[k2]["sum"],mp.nstr(rel(v["talbot"],s[k2]["sum"]),3)),"talbot_%s.json"%t)
# boschB Theorem bs:thm:gpos
truth("boschB:375","Sigma_off < 2e-11","< 2e-11",mp.mpf(grab(85,r"eps_off <= ([0-9.e+-]+)"))<mp.mpf('2e-11'),grab(85,r"eps_off <= ([0-9.e+-]+)"),"line 85")
truth("boschB:375","kappa < 8e-9","< 8e-9",mp.mpf(grab(85,r"kappa=([0-9.e+-]+)"))<mp.mpf('8e-9'),grab(85,r"kappa=([0-9.e+-]+)"),"line 85")
truth("boschB:375","A > 0.99989","> 0.99989",mp.mpf(grab(85,r"A >= ([0-9.]+)"))>mp.mpf('0.99989'),grab(85,r"A >= ([0-9.]+)"),"line 85")
truth("boschB:377","each pair term below exp(-1.8e13)","log bound < -1.8e13",mp.mpf(grab(86,r"b=1/2: ([-0-9.e+]+)"))< -mp.mpf('1.8e13'),grab(86,r"b=1/2: ([-0-9.e+]+)"),"line 86")
for n in [81,82,83]:
    a,b=re.findall(r"=([0-9.]+)",L(n))[-2:]; truth("boschB:412","pair kernel: FT equals q(t^2) (line %d)"%n,"equal",a==b,a+" | "+b,"line %d"%n)
# boschB Section on W
truth("boschB:390","(log g)'' <= -2.2 on [0,1]","max <= -2.2",mp.mpf(grab(77,r": ([-0-9.e+]+)$"))<=mp.mpf('-2.2'),grab(77,r": ([-0-9.e+]+)$"),"line 77")
rounds("boschA:192, boschB:397","gamma_2-gamma_1","6.8873",grab(68,r"gamma_2-gamma_1 = ([0-9.]+)"),src="line 68")
B=re.findall(r"'(-?[0-9.]+)'",L(68))
for q,v,k in zip(["359.52","-28007","392918"],B,[1,2,3]): rounds("boschB:412","b_%d"%k,q,v,src="line 68")
rounds("boschB:412","g(0)","1.8994210510",grab(67,r"quadrature: ([0-9.]+)"),src="line 67")
y1=[mp.mpf(x) for x in re.findall(r"g=\(([-0-9.e]+)\+0j\)",L(69))]
ca=[complex(x.replace(' ','')) for x in re.findall(r"g=(\([-0-9.e+j]+\))",L(70))]
r1=abs(y1[0]-y1[1])/abs(y1[0]); r2=abs(ca[0]-ca[1])/abs(ca[0])
truth("boschB:412","quadrature and series agree to eleven digits at u=1 and u=1.2+0.3i","rel diff <= 5e-11",r1<=5e-11 and r2<=5e-11,"%s, %s"%(mp.nstr(r1,3),mp.nstr(r2,3)),"lines 69, 70")
rounds("boschB:412","v_* = pi/(gamma_2-gamma_1)","0.45614",grab(68,r"pi/\(gamma_2-gamma_1\) = ([0-9.]+)"),src="line 68")
truth("boschB:412","scan of 4840 points: no value below -1.3e-15","4840 points, min >= -1.3e-15",grab(76,r"\((\d+) points\)")=="4840" and mp.mpf(grab(76,r"= ([-0-9.e+]+) at"))>=mp.mpf('-1.3e-15'),grab(76,r"= ([-0-9.e+]+) at"),"line 76")
rounds("boschB:412","min at v=0.999 v_* (40 digits)","1.0e-27",grab(78,r"= ([-0-9.e+]+) at"),src="line 78")
rounds("boschB:412","min at v=1.001 v_*","-6.5e-8",grab(79,r"= ([-0-9.e+]+) at"),src="line 79")
rounds("boschB:412","min at v=1.01 v_*","-3.6e-5",grab(74,r"= ([-0-9.e+]+) at"),src="line 74")
rounds("boschB:412","min at v=1.1 v_*","-2.1e-2",grab(75,r"= ([-0-9.e+]+) at"),src="line 75")
a,b=re.findall(r"([0-9.]+e-\d+)",L(80)); truth("boschB:412","asymptotic formula at x=8, v=v_*/2 to eight digits","equal to 8 digits",a==b and len(a.split('e')[0].replace('.',''))>=8,a+" | "+b,"line 80")
# boschA: A_n = 2 b_n/gamma_n and the three term accuracy of the series (block 9's quadrature, chi_nodes.npz)
with mp.workdps(30):
    gam=[mp.im(mp.zetazero(n)) for n in (1,2,3)]
for q,v,k in zip(["50.9","-2664.6","31420"],B,[0,1,2]): rounds("boschA:128","A_%d = 2 b_%d/gamma_%d"%(k+1,k+1,k+1),q,2*mp.mpf(v)/gam[k],src="line 68")
d=np.load("fixed/chi_nodes.npz"); T=d["T"]; Wt=d["W"]; lc=d["logchi"]
def g_int(y):  # block 9, verbatim
    y=complex(y)
    e1=np.exp(lc+1j*T*y); e2=np.exp(lc-1j*T*y)
    g=np.sum(Wt*(e1+e2)/2)/np.pi
    gp=np.sum(Wt*T*(e1-e2)/(2j))*(-1)/np.pi
    return g,gp
for y,q in [(1.0,"3e-3"),(1.2,"1.4e-4")]:
    gi=mp.mpf(g_int(y)[0].real); s3=sum(mp.mpf(bb)*mp.e**(-gg*y) for bb,gg in zip(B,gam))
    rounds("boschA:128","three terms against the Fourier integral at y=%s, relative error"%y,q,abs(s3-gi)/abs(gi),src="block 9 quadrature (fixed/chi_nodes.npz), b_n of line 68")
# boschB Table bs:tab:gaps: beta_* from block 13 (lines 103 to 108) and from the printed code on the stale cases
BETA={"[]":"1.000","[2]":"2.000","[2, 4]":"2.000","[2, 5]":"2.000","[2, 3]":"3.000","[2, 3, 6]":"3.000","[2, 3, 5, 6]":"3.000","[2, 3, 4]":"4.000","[3]":"1.269","[4, 5]":"1.317","[3, 4]":"1.483"}
def gapline(l):
    return re.search(r"removed (\[[0-9, ]*\])",l).group(1),re.search(r"first gap (\d+), max gap (\d+)",l).groups(),re.search(r"pi/threshold = ([0-9.]+)",l).group(1)
seen={}
for l,where in [(x,"block 13 (fixed run)") for x in fixed[89:95]]+[(x,"printed block 13 code, supp_gaps_table.py") for x in sgap]:
    rem,(fg,mg),beta=gapline(l)
    if rem in seen:
        truth("boschB:437","Table bs:tab:gaps %s: block 13 and supp_gaps_table.py agree"%rem,"same line",seen[rem]==TIM.sub('',l),beta,where); continue
    seen[rem]=TIM.sub('',l)
    rounds("boschB:426-437","Table bs:tab:gaps, removed %s (first gap %s, largest gap %s): beta_*"%(rem,fg,mg),BETA[rem],beta,src=where)
truth("boschB:418","bisection grid 0 <= x <= 16","xs = k/5, k=0..80",any("xs=[mp.mpf(k)/5 for k in range(0,81)]" in l for l in src),"xs=[mp.mpf(k)/5 for k in range(0,81)]","script line 399")
# boschB before Question bs:q:gaprefined: line 109 (reconstructed) and arithmetic
rounds("boschB:444","first gap among the first 1000 zeros","6.8873",grab(109,r"first gap ([0-9.]+)"),src="line 109 (reconstructed)")
rounds("boschB:444","next largest gap gamma_4-gamma_3","5.414",grab(109,r"n<1000: ([0-9.]+)"),src="line 109 (reconstructed)")
truth("boschB:444","the next largest gap is at n=3 (gamma_4-gamma_3)","n=3",L(109).endswith("at n= 3"),L(109).split(":")[-1].strip(),"line 109 (reconstructed)")
ex=[l for l in lines("supp109/stdout.txt") if l.startswith("(extra)")]
if ex:
    g1000=re.search(r"gamma_1000 = ([0-9.]+)",ex[0]).group(1); rounds("boschB:444","height of zero 1000","1419",g1000,src="supp_line109.py extra line")
    truth("boschB:444","first 1000 zeros (not n<1000): same largest later gap","5.414 at n=3",re.search(r"largest later gap ([0-9.]+) at n=(\d+)",ex[0]).group(2)=="3",ex[0][9:120],"supp_line109.py extra line")
with mp.workdps(30):
    h=mp.mpf(grab(68,r"gamma_2-gamma_1 = ([0-9.]+)"))
    F=lambda T: T/(2*mp.pi)*mp.log(T/(2*mp.pi*mp.e)); E=lambda T: mp.mpf('0.112')*mp.log(T)+mp.mpf('0.278')*mp.log(mp.log(T))+mp.mpf('2.510')+mp.mpf('0.2')/T
    Dt=lambda T: F(T+h)-F(T)-E(T)-E(T+h)
    lo,hi=mp.mpf(2000),mp.mpf(10)**6
    for _ in range(100):
        m=(lo+hi)/2
        if Dt(m)>0: hi=m
        else: lo=m
    rounds("boschB:444","Trudgian: every (T, T+6.8873] contains an ordinate for T >= this","1.34e4",hi,src="arithmetic")
    truth("boschB:444","mean spacing 2 pi/log(T/2 pi) below 1.2 on [1419, 1.34e4]","< 1.2",2*mp.pi/mp.log(mp.mpf(1419)/(2*mp.pi))<mp.mpf('1.2'),mp.nstr(2*mp.pi/mp.log(mp.mpf(1419)/(2*mp.pi)),4),"arithmetic")
# app:bosch: Table tab:app:bosch (block starts in the script, output lines per block)
truth("app:bosch","block starts in the script: 1,51,87,117,149,183,200,220,241,295,319,344,389","as in Table tab:app:bosch",starts==[1,51,87,117,149,183,200,220,241,295,319,344,389],str(starts),"src/verify_bosch_xi_v7.py")
truth("app:bosch","output lines per file of the fixed run: 10,14,12,10,17,14,6,6,6","1-10, 11-24, 25-36, 40-49, 50-66, 67-80, 81-86, 87-92, 103-108",[counts[f] for f in files]==[10,14,12,10,17,14,6,6,6],str([counts[f] for f in files]),"fixed/*_output.txt")
res=dict(part1=dict(groups=groups,lines=per_line),part2=p2,part3=checks,residue=R,talbot=talbot)
json.dump(res,open("compare_bosch_v7.json","w"),indent=1,default=str)
print("== Part 1: archived output against the rerun")
for k,g in groups.items(): print(" ",k,{x:y for x,y in g.items() if x!="timing_pairs"})
for d in per_line:
    if not d["identical"]: print("   line %d (%s): numbers %d, differing %d, text same %s, timings %s -> %s"%(d["line"],d["source"],d["n_numbers"],len(d["differences"]),d["text_same"],d["timings_archive"],d["timings_rerun"]))
print("== Part 2: fixed copy against the unchanged archived script:",{k:v for k,v in p2.items()})
print("== Part 3: %d checks, %d fail"%(len(checks),sum(not c["ok"] for c in checks)))
for c in checks: print("  %s %-16s %s: quoted %s, got %s (%s)"%("ok  " if c["ok"] else "FAIL",c["where"],c["what"],c["quoted"],c["value"],c["source"]))
