import mpmath as mp
out=[]
def P(*a):
    s=" ".join(str(x) for x in a); print(s); out.append(s)

def first_violation(pair, w0, order, dps=60):
    with mp.workdps(dps):
        c = mp.taylor(pair, mp.mpf(w0), order)
        for n,cn in enumerate(c):
            if (-1)**n*cn < 0: return n
    return None

def pair_of(h, u):
    def F(w):
        v = (w+mp.sqrt(w*w-4))/2
        return h(u*v)*h(u/v)
    return F

# 1. Bosch family: density of T_alpha^beta, q=1/beta, h(x)=x^(q-1)/(x^(2q)+2 cos(pi a) x^q+1)
def bosch_h(a,q):
    c=mp.cos(mp.pi*a)
    return lambda x: x**(q-1)/(x**(2*q)+2*c*x**q+1)
for (a,q) in [(0.4,0.60),(0.4,0.75),(0.4,0.9),(0.6,0.3)]:
    a=mp.mpf(a); q=mp.mpf(q)
    res=[first_violation(pair_of(bosch_h(a,q),mp.mpf(u)),2.5,60) for u in [0.3,1,3,10]]
    P("Bosch alpha=",a," q=1/beta=",q," threshold 1-alpha=",1-a," first violations (u=0.3,1,3,10):",res)

# 2. toy reciprocal Xi with one off line pair gamma=14 e^{i delta}
delta=mp.mpf('0.1'); g=14*mp.expj(delta)
phi_toy=lambda s: 1/((1+s/g**2)*(1+s/mp.conj(g)**2)*(1+s/21.0**2))
pstar=1-2*delta/mp.pi
P("toy: delta=0.1, predicted p* = 1-2delta/pi =",mp.nstr(pstar,6))
for p in [0.85,0.92,0.95,1.0]:
    p=mp.mpf(p); h=lambda x,p=p: mp.re(phi_toy(x**p))
    res=[first_violation(pair_of(h,mp.mpf(u)),2.5,60) for u in [50,150,200,400]]
    P("  p=",p," first violations (u=50,150,200,400):",res)

# 3. ray monotonicity of the true phi: Im(e^{i th} psi(r e^{i th})) >= 0
mp.mp.dps=20
def dlogxi(z):
    return 1/z+1/(z-1)-mp.log(mp.pi)/2+mp.digamma(z/2)/2+mp.zeta(z,derivative=1)/mp.zeta(z)
def psi(s):
    r=mp.sqrt(s); return dlogxi(mp.mpf(1)/2+r)/(2*r)
mn=mp.inf; cnt=0
for R in [0.5,1,5,20,50,100,150,190,200,210,250,400,800,1500,3000]:
    for j in range(1,60):
        th=mp.pi*j/60; val=mp.im(mp.expj(th)*psi(R*mp.expj(th))); cnt+=1
        if val<mn: mn=val; where=(R,j)
P("ray monotonicity: min of Im(e^{i theta} psi(r e^{i theta})) over",cnt,"points:",mp.nstr(mn,5),"at",where)
open("verify_v2_output.txt","w").write("\n".join(out)+"\n")
import mpmath as mp
mp.mp.dps=25
out=[]
def P(*a):
    s=" ".join(str(x) for x in a); print(s); out.append(s)
rs=[mp.mpf(10)**(k/4) for k in range(-24,25)]
ths=[mp.pi*j/200 for j in range(1,200)]
def min_elasticity(E):
    m=mp.inf
    for r in rs:
        for t in ths:
            v=mp.im(E(r*mp.expj(t)))
            if v<m: m=v; at=(mp.nstr(r,3),mp.nstr(t/mp.pi,3))
    return m,at
# Bosch family: E(x) = -x h'/h, h = x^(q-1)/D, D = x^(2q)+2c x^q+1
def E_bosch(a,q):
    c=mp.cos(mp.pi*a)
    def E(x):
        X=x**q
        return -(q-1) + (2*q*X*X+2*c*q*X)/(X*X+2*c*X+1)
    return E
for a,q in [(0.4,0.6),(0.4,0.55),(0.4,0.75),(0.5,0.5),(0.5,0.6),(0.6,0.3),(0.3,0.7)]:
    a=mp.mpf(a);q=mp.mpf(q)
    poles = (1-a)/q < 1
    m,at=min_elasticity(E_bosch(a,q))
    P("Bosch alpha=%s q=%s: zeros of D in slit plane=%s, min Im(-x h'/h)=%s at (r,theta/pi)=%s, HCM predicted=%s"%(a,q,poles,mp.nstr(m,4),at, (a<=0.5 and q<=1-a)))
# power family of a toy reciprocal Xi with one off line pair
delta=mp.mpf('0.1'); g=14*mp.expj(delta); cs=[g**2,mp.conj(g)**2,mp.mpf(21)**2]
P("toy with gamma=14 e^{0.1 i} and 21: p* = 1 - 2 delta/pi =",mp.nstr(1-2*delta/mp.pi,6))
for p in [0.8,0.9,0.93,0.94,0.97,1.0]:
    p=mp.mpf(p)
    E=lambda x,p=p: p*sum(x**p/(x**p+c) for c in cs)
    poles = any(abs(mp.arg(-c))<p*mp.pi for c in cs if mp.im(c)!=0)
    m,at=min_elasticity(E)
    P("  p=%s: poles in slit plane=%s, min Im(-x h'/h)=%s at %s"%(p,poles,mp.nstr(m,4),at))
open("verify_v3_output.txt","w").write("\n".join(out)+"\n")
import mpmath as mp, time
mp.mp.dps=20
out=[]
def P(*a):
    s=" ".join(str(x) for x in a); print(s,flush=True); out.append(s)
def dlogxi(z):
    return 1/z+1/(z-1)-mp.log(mp.pi)/2+mp.digamma(z/2)/2+mp.zeta(z,derivative=1)/mp.zeta(z)
def xi(z): return z*(z-1)/2*mp.pi**(-z/2)*mp.gamma(z/2)*mp.zeta(z)
x0=xi(mp.mpf(1)/2)
def eps_phi(s):
    z=mp.sqrt(s); return z/2*dlogxi(mp.mpf(1)/2+z)
def u(p,tau):  # Thorin density of phi(s^p), computed from xi'/xi only
    return p/(mp.pi*tau)*mp.im(eps_phi(tau**p*mp.expj(p*mp.pi)))
gam=[mp.im(mp.zetazero(n)) for n in range(1,31)]
def u_zeros(p,tau,N=30):  # Lamperti superposition over first N zeros
    return p/(mp.pi*tau)*sum(g**2*tau**p*mp.sin(p*mp.pi)/(tau**(2*p)+2*g**2*tau**p*mp.cos(p*mp.pi)+g**4) for g in gam[:N])
for p in [mp.mpf('0.5'),mp.mpf('0.9')]:
    for tau in [1,100,(gam[0]**(2/p)),1e4]:
        tau=mp.mpf(tau)
        P("p=%s tau=%s  u_p from xi'/xi = %s   Lamperti sum over 30 zeros = %s"%(mp.nstr(p,3),mp.nstr(tau,6),mp.nstr(u(p,tau),8),mp.nstr(u_zeros(p,tau),8)))
# Laplace transform reconstruction: log phi(s^p) = -int log(1+s/tau) u_p(tau) dtau
for p in [mp.mpf('0.5'),mp.mpf('0.9')]:
    for s in [mp.mpf(1),mp.mpf(10)]:
        t0=time.time()
        bps=[0,1e-3,1,10,100]+[g**(2/p) for g in gam[:6]]+[1e5,1e7,mp.inf]
        bps=sorted(set([mp.mpf(b) for b in bps]))
        I=mp.quad(lambda t: mp.log(1+s/t)*u(p,t), bps, maxdegree=8)
        lhs=mp.log(x0/xi(mp.mpf(1)/2+mp.sqrt(s**p)))
        P("p=%s s=%s: log phi(s^p)=%s, -int log(1+s/t)u_p = %s (%.0fs)"%(mp.nstr(p,3),s,mp.nstr(lhs,10),mp.nstr(-I,10),time.time()-t0))
open("thorin_p_output.txt","w").write("\n".join(out)+"\n")
import mpmath as mp
mp.mp.dps=25
out=[]
def P(*a):
    s=" ".join(str(x) for x in a); print(s,flush=True); out.append(s)
# (a) Cauchy analogue X_C = sum E_k/(2k-1)^2: density = (4/pi) sum (-1)^(k-1)(2k-1)e^{-(2k-1)^2 t} = (4/pi) eta(4it/pi)^3
def f_series(t): return 4/mp.pi*mp.nsum(lambda k:(-1)**(k-1)*(2*k-1)*mp.e**(-(2*k-1)**2*t),[1,mp.inf])
def eta(tau):
    q=mp.expj(2*mp.pi*tau)
    return mp.expj(2*mp.pi*tau/24)*mp.qp(q)
for t in [mp.mpf('0.05'),mp.mpf('0.3'),mp.mpf(1),mp.mpf(3)]:
    P("t=%s  series density %s   (4/pi) eta(4it/pi)^3 %s"%(t,mp.nstr(f_series(t),15),mp.nstr(mp.re(4/mp.pi*eta(4j*t/mp.pi)**3),15)))
P("check LT at s=2: int e^{-st} f = %s, sech(pi sqrt2/2) = %s"%(mp.nstr(mp.quad(lambda t: mp.e**(-2*t)*f_series(t),[0,0.05,1,mp.inf]),12),mp.nstr(mp.sech(mp.pi*mp.sqrt(2)/2),12)))
# (b) prime region sufficient condition
def G(w): return 1/w+1/(w-1)-mp.log(mp.pi)/2+mp.digamma(w/2)/2
def mzl(sig): return -mp.zeta(sig,derivative=1)/mp.zeta(sig)
def dlogxi(z): return G(z)+mp.zeta(z,derivative=1)/mp.zeta(z)
for p in ['0.5','0.9','0.99','0.999']:
    p=mp.mpf(p); th=p*mp.pi/2; c=mp.cos(th)
    t_min=(mp.mpf(1)/2)/c*mp.mpf('1.0001')
    ts=[t_min*mp.mpf(10)**(k/40) for k in range(0,40*14)]
    last_bad=None; exact_min=mp.inf
    for t in ts:
        z=t*mp.expj(th); w=mp.mpf(1)/2+z
        D=mp.im(z*G(w))-t*mzl(mp.mpf(1)/2+t*c)
        if D<=0: last_bad=t
    t1=last_bad if last_bad else t_min
    P("p=%s: prime region starts at |z|=%s; sufficient inequality holds for all |z|>%s  (tau_1 = %s), checked to |z|=%s"%(mp.nstr(p,4),mp.nstr(t_min,5),mp.nstr(t1,5),mp.nstr(t1**(2/p),5),mp.nstr(ts[-1],3)))
# asymptotic constant check: D/t vs 0.5[sin th log(t/2pi)+th cos th] at large t
p=mp.mpf('0.9'); th=p*mp.pi/2; t=mp.mpf(10)**8; z=t*mp.expj(th)
P("asymptotic check p=0.9 |z|=1e8: Im[zG]/|z| = %s, 0.5[sin th log(|z|/2pi)+th cos th] = %s"%(mp.nstr(mp.im(z*G(mp.mpf(1)/2+z))/t,10),mp.nstr((mp.sin(th)*mp.log(t/(2*mp.pi))+th*mp.cos(th))/2,10)))
open("open_problems_output.txt","w").write("\n".join(out)+"\n")
import mpmath as mp, time
mp.mp.dps=30
out=[]
def P(*a):
    s=" ".join(str(x) for x in a); print(s,flush=True); out.append(s)
def xi(z): return z*(z-1)/2*mp.pi**(-z/2)*mp.gamma(z/2)*mp.zeta(z)
half=mp.mpf(1)/2; x0=xi(half)
N=200
rho=[mp.zetazero(n) for n in range(1,N+1)]
gam=[mp.im(r) for r in rho]
def c_coef(r):
    g=mp.im(r)
    dxi=r*(r-1)/2*mp.pi**(-r/2)*mp.gamma(r/2)*mp.zeta(r,derivative=1)   # xi'(rho)
    Xip=1j*dxi                                                          # Xi'(gamma)
    return mp.re(-2*g*x0/Xip)
cs=[c_coef(r) for r in rho]
P("c_1..c_4 =",[mp.nstr(c,8) for c in cs[:4]],"  (c_1/gamma_1^2 =",mp.nstr(cs[0]/gam[0]**2,8),", compare A_1 ~ 50.92)")
P("signs alternate for n<=200:", all(cs[i]*cs[i+1]<0 for i in range(N-1)))
gaps=[gam[i+1]**2-gam[i]**2 for i in range(N-1)]
i=min(range(N-1),key=lambda k:gaps[k])
P("min gap of lambda_n=gamma_n^2 over n<200: %s at n=%d (gamma gap %s)"%(mp.nstr(gaps[i],6),i+1,mp.nstr(gam[i+1]-gam[i],6)))
rat=[mp.log(abs(c))/g**2 for c,g in zip(cs,gam)]
P("log|c_n|/lambda_n at n=1,10,50,100,200:",[mp.nstr(rat[k],5) for k in [0,9,49,99,199]])
P("pi/(4 gamma_n) at same n:",[mp.nstr(mp.pi/(4*gam[k]),5) for k in [0,9,49,99,199]])
# density series vs numerical inverse Laplace transform
mp.mp.dps=25
phi=lambda s: x0/xi(half+mp.sqrt(s))
for t in [mp.mpf('0.01'),mp.mpf('0.003')]:
    ser=sum(c*mp.e**(-g**2*t) for c,g in zip(cs,gam))
    t0=time.time(); inv=mp.invertlaplace(phi,t,method='talbot')
    P("t=%s: residue series %s, Talbot inverse Laplace %s (%.0fs)"%(t,mp.nstr(ser,12),mp.nstr(inv,12),time.time()-t0))
# Cauchy side: lambda=(2k-1)^2, gaps 8k; sigma_c = limsup log(2k-1)/(2k-1)^2 = 0
P("Cauchy side: gaps (2k+1)^2-(2k-1)^2 = 8k >= 8; log|c_k|/lambda_k -> 0")
open("op5_output.txt","w").write("\n".join(out)+"\n")
import mpmath as mp
mp.mp.dps=90
def xi(z): return z*(z-1)/2*mp.pi**(-z/2)*mp.gamma(z/2)*mp.zeta(z)
half=mp.mpf(1)/2; x0=xi(half)
rho=[mp.zetazero(n) for n in range(1,201)]
def c_coef(r):
    g=mp.im(r); dxi=r*(r-1)/2*mp.pi**(-r/2)*mp.gamma(r/2)*mp.zeta(r,derivative=1)
    return g, mp.re(-2*g*x0/(1j*dxi))
gc=[c_coef(r) for r in rho]
phi=lambda s: x0/xi(half+mp.sqrt(s))
out=[]
for t in ['0.02','0.005','0.003']:
    t=mp.mpf(t)
    ser=sum(c*mp.e**(-g**2*t) for g,c in gc)
    with mp.workdps(40): inv=mp.invertlaplace(phi,t,method='talbot')
    s="t=%s: residue series (200 zeros, 90 digits) %s, Talbot %s"%(t,mp.nstr(ser,10),mp.nstr(inv,10)); print(s); out.append(s)
open("op5_output.txt","a").write("\n".join(out)+"\n")
import mpmath as mp
mp.mp.dps=25
out=[]
def mzl(sig): return -mp.zeta(sig,derivative=1)/mp.zeta(sig)
def L(t,th):
    z=t*mp.expj(th); w=mp.mpf(1)/2+z; sw=mp.re(w)
    main=mp.im(z*(1/(2*w)+1/(w-1)+mp.log(w/(2*mp.pi))/2))
    return main - t/(6*sw**2) - t*mzl(sw)
for p in ['0.5','0.9','0.99','0.999','0.9999']:
    p=mp.mpf(p); th=p*mp.pi/2; c=mp.cos(th)
    tmin=(mp.mpf(1)/2)/c*mp.mpf('1.0001')
    ts=[tmin*mp.mpf(10)**(k/50) for k in range(0,50*14)]
    bad=[t for t in ts if L(t,th)<=0]
    t1=bad[-1] if bad else tmin
    # check also closed form bound -zeta'/zeta(s) < 1/(s-1) at these sigma
    s="p=%s: rigorous bound positive for |z|>%s, tau_1=%s (strip exit %s)"%(mp.nstr(p,5),mp.nstr(t1,5),mp.nstr(t1**(2/p),5),mp.nstr(tmin,5))
    print(s); out.append(s)
ok=all(mzl(mp.mpf(1)+mp.mpf(10)**(-k/10))<1/(mp.mpf(10)**(-k/10)) for k in range(-20,60))
s="check -zeta'/zeta(sigma) < 1/(sigma-1) on sigma-1 in [1e-6,100]: %s"%ok; print(s); out.append(s)
open("op5_output.txt","a").write("\n".join(out)+"\n")
import mpmath as mp, numpy as np, time
mp.mp.dps=30
t0=time.time()
half=mp.mpf(1)/2
def logxi(z): return mp.log(abs(z*(z-1)/2))-z/2*mp.log(mp.pi)+mp.loggamma(z/2)+mp.log(abs(mp.zeta(z)))
lx0=logxi(half)
# composite Gauss-Legendre on [0,Tmax]
Tmax=700.0; panels=1400; n=16
xg,wg=np.polynomial.legendre.leggauss(n)
edges=np.linspace(0,Tmax,panels+1)
T=[];W=[]
for a,b in zip(edges[:-1],edges[1:]):
    T+=list((b-a)/2*xg+(a+b)/2); W+=list((b-a)/2*wg)
T=np.array(T);W=np.array(W)
def lchi(t):
    z=half+mp.mpf(t)
    if abs(z-1)<mp.mpf('1e-6'): z=z+mp.mpf('2e-6')
    return float(lx0-logxi(z))
logchi=np.array([lchi(t) for t in T])
print("chi nodes computed", round(time.time()-t0,1),"s; logchi(700)=",logchi[-1])
np.savez("chi_nodes.npz",T=T,W=W,logchi=logchi)
import numpy as np, mpmath as mp, time
d=np.load("chi_nodes.npz"); T=d["T"]; Wt=d["W"]; lc=d["logchi"]
def g_int(y):   # g(y)=(1/pi) int chi(t) cos(t y) dt, and g'(y)
    y=complex(y)
    e1=np.exp(lc+1j*T*y); e2=np.exp(lc-1j*T*y)
    g=np.sum(Wt*(e1+e2)/2)/np.pi
    gp=np.sum(Wt*T*(e1-e2)/(2j))*(-1)/np.pi
    return g,gp
out=[]
def P(*a):
    s=" ".join(str(x) for x in a); print(s,flush=True); out.append(s)
g0,_=g_int(0); P("g(0) from quadrature:",g0.real,"(first version: 1.89942105102)")
# series coefficients b_n = -Xi(0)/Xi'(gamma_n)
mp.mp.dps=30
half=mp.mpf(1)/2
def xi(z): return z*(z-1)/2*mp.pi**(-z/2)*mp.gamma(z/2)*mp.zeta(z)
x0=xi(half)
rho=[mp.zetazero(n) for n in range(1,201)]
gam=[mp.im(r) for r in rho]
b=[]
for r in rho:
    dxi=r*(r-1)/2*mp.pi**(-r/2)*mp.gamma(r/2)*mp.zeta(r,derivative=1)
    b.append(mp.re(-x0/(1j*dxi)))
P("b_1,b_2,b_3 =",[mp.nstr(v,8) for v in b[:3]],"; gamma_2-gamma_1 =",mp.nstr(gam[1]-gam[0],10),"; pi/(gamma_2-gamma_1) =",mp.nstr(mp.pi/(gam[1]-gam[0]),10))
def g_ser(y):
    y=mp.mpc(y); g=sum(bb*mp.e**(-gg*y) for bb,gg in zip(b,gam)); gp=sum(-gg*bb*mp.e**(-gg*y) for bb,gg in zip(b,gam))
    return complex(g),complex(gp)
for y in [1.0,1.2+0.3j]:
    a=g_int(y); c=g_ser(y); P("overlap y=%s: integral g=%s, series g=%s"%(y,a[0],c[0]))
# Pick scan
vstar=float(mp.pi/(gam[1]-gam[0]))
def pick(y):
    g,gp=(g_int(y) if y.real<1.1 else g_ser(y))
    return (-gp/g).imag, abs(g)
for frac in [0.5,0.9,0.99,1.01,1.1]:
    v=frac*vstar; mn=1e99; mng=1e99
    for x in np.linspace(0,6,241):
        val,ag=pick(complex(x,v))
        if val<mn: mn=val; at=x
        mng=min(mng,ag)
    P("v=%.4f (=%.2f v*): min over x in [0,6] of Im(-g'/g) = %.4e at x=%.3f; min |g| = %.3e"%(v,frac,mn,at,mng))
mn=1e99
for v in np.linspace(0.005,0.99*vstar,40):
    for x in np.linspace(0,6,121):
        val,_=pick(complex(x,v))
        if val<mn: mn=val; at=(x,v)
P("full strip scan 0<v<0.99 v*, 0<=x<=6 (4840 points): min Im(-g'/g) = %.4e at (x,v)=(%.3f,%.4f)"%(mn,at[0],at[1]))
# log concavity on real line: (log g)'' <= 0
lcv=[]
for x in np.linspace(0,1.0,101):
    h=1e-3; f=lambda u: np.log(g_int(u)[0].real)
    lcv.append((f(x+h)-2*f(x)+f(x-h))/h**2)
P("max (log g)'' on [0,1]:",max(lcv))
open("gW_output.txt","w").write("\n".join(out)+"\n")
import mpmath as mp
mp.mp.dps=40
half=mp.mpf(1)/2
def xi(z): return z*(z-1)/2*mp.pi**(-z/2)*mp.gamma(z/2)*mp.zeta(z)
x0=xi(half)
rho=[mp.zetazero(n) for n in range(1,121)]
gam=[mp.im(r) for r in rho]
b=[mp.re(-x0/(1j*r*(r-1)/2*mp.pi**(-r/2)*mp.gamma(r/2)*mp.zeta(r,derivative=1))) for r in rho]
vs=mp.pi/(gam[1]-gam[0])
def im_el(y):
    g=sum(bb*mp.e**(-gg*y) for bb,gg in zip(b,gam)); gp=sum(-gg*bb*mp.e**(-gg*y) for bb,gg in zip(b,gam))
    return mp.im(-gp/g)
out=[]
for frac in ['0.999','1.001']:
    v=mp.mpf(frac)*vs; mn=mp.inf
    for k in range(0,81):
        x=mp.mpf(1.1)+k*mp.mpf('0.1')
        val=im_el(mp.mpc(x,v))
        if val<mn: mn=val; at=x
    s="v=%s v*: min over x in [1.1,9.1] of Im(-g'/g) (40 digits, 120 zeros) = %s at x=%s"%(frac,mp.nstr(mn,6),mp.nstr(at,3)); print(s); out.append(s)
# asymptotic prediction of Im at large x: (g2-g1)|b2/b1| e^{-D x} sin(D v)
x=mp.mpf(8); v=mp.mpf('0.5')*vs; D=gam[1]-gam[0]
s="x=8,v=v*/2: Im(-g'/g)=%s, leading prediction %s"%(mp.nstr(im_el(mp.mpc(x,v)),8),mp.nstr(D*abs(b[1]/b[0])*mp.e**(-D*x)*mp.sin(D*v),8)); print(s); out.append(s)
open("gW_output.txt","a").write("\n".join(out)+"\n")
import mpmath as mp
mp.mp.dps=25
out=[]
def P(*a):
    s=" ".join(str(x) for x in a); print(s); out.append(s)
# pair kernel check
rho_,d=mp.mpf(14),mp.mpf('0.05'); g=rho_*mp.expj(d); a,bb=mp.re(g),mp.im(g)
k=lambda x: rho_*mp.e**(-a*abs(x))*mp.sin(bb*abs(x)+d)/(2*mp.sin(2*d))
for t in [0,3,10]:
    ft=2*mp.quad(lambda x: k(x)*mp.cos(t*x),[0,1,5,mp.inf])
    q=1/((1+t**2/g**2)*(1+t**2/mp.conj(g)**2))
    P("pair kernel: t=%s  FT=%s  q(t^2)=%s"%(t,mp.nstr(ft,12),mp.nstr(mp.re(q),12)))
P("first negative lobe starts at x=(pi-delta)/Im(gamma) =",mp.nstr((mp.pi-d)/bb,6), " k there:",mp.nstr(k((mp.pi-d)/bb+0.01),4))
# unconditional constants with T0=3e12
T0=mp.mpf('3e12'); g1=mp.mpf('14.134725141734693790')
eps=2*(mp.log(T0)+1)/T0
kap=2*g1**2*eps
A=mp.e**kap*mp.erfc(mp.sqrt(kap))
P("eps_off <= %s, kappa=%s, A >= %s"%(mp.nstr(eps,4),mp.nstr(kap,4),mp.nstr(A,10)))
# B bound for worst case: a=T0, b=1/2
a=T0; b=mp.mpf(1)/2; dl=mp.atan(b/a); rho=mp.sqrt(a*a+b*b)
xg=(mp.pi-dl)/b
logB=mp.log(rho/mp.sin(2*dl)/(a-g1))-(a-g1)*xg
P("log of single pair C(k^-) bound at a=T0,b=1/2: %s  (so B < exp(-9e12))"%mp.nstr(logB,6))
open("gpos_output.txt","w").write("\n".join(out)+"\n")
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
import mpmath as mp, time
mp.mp.dps=30
out=[]
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
cases=[[2,3],[2,3,5,6],[3,4],[4,5],[3],[]]
for rem in cases:
    fr=freqs(rem); gaps=[fr[i+1]-fr[i] for i in range(len(fr)-1)]
    t0=time.time(); th=threshold(rem)
    thv = th if th is not None else mp.pi
    P("removed %s: freqs %s..., gaps %s...: first gap %d, max gap %d; threshold strip %s; pi/threshold = %s; pi/first gap = %s (%.0fs)"%(rem,fr,gaps,gaps[0],max(gaps),mp.nstr(thv,6),mp.nstr(mp.pi/thv,6),mp.nstr(mp.pi/gaps[0],5),time.time()-t0))
open("gaptest3_output.txt","w").write("\n".join(out)+"\n")
