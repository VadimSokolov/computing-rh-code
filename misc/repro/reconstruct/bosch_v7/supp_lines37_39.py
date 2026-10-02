# Reconstruction, not the authors' code (misc/repro/reconstruct/bosch_v7). Output lines 37 to 39 of
# code_bosch/verify_bosch_xi_v7_output.txt have no generating code in the archive (Appendix app:bosch says so);
# Chapter ch:boschB quotes them after Theorem bs:thm:ggc. This script recomputes them with the functions of
# block 3 of code_bosch/verify_bosch_xi_v7.py (archive lines 87 to 116: mpmath at 20 digits, u from xi'/xi,
# u_zeros as the Lamperti sum), copied verbatim except that the zero list has 200 zeros instead of 30.
# Line 39 needs a grid of 148 points in p and tau; the grid below (4 values of p times 37 values
# tau = 10^(k/3), k = -12..24, which spans [1e-4, 1e8]) is a guess that matches the printed description.
# The printed minimum sits at the end point tau = 1e8, so it does not depend on the interior of the grid.
import mpmath as mp
mp.mp.dps=20
out=[]
def P(*a):
    s=" ".join(str(x) for x in a); print(s,flush=True); out.append(s)
def dlogxi(z):
    return 1/z+1/(z-1)-mp.log(mp.pi)/2+mp.digamma(z/2)/2+mp.zeta(z,derivative=1)/mp.zeta(z)
def eps_phi(s):
    z=mp.sqrt(s); return z/2*dlogxi(mp.mpf(1)/2+z)
def u(p,tau):  # Thorin density of phi(s^p), computed from xi'/xi only
    return p/(mp.pi*tau)*mp.im(eps_phi(tau**p*mp.expj(p*mp.pi)))
gam=[mp.im(mp.zetazero(n)) for n in range(1,201)]
def u_zeros(p,tau,N=30):  # Lamperti superposition over first N zeros
    return p/(mp.pi*tau)*sum(g**2*tau**p*mp.sin(p*mp.pi)/(tau**(2*p)+2*g**2*tau**p*mp.cos(p*mp.pi)+g**4) for g in gam[:N])
# lines 37 and 38: the two points of block 3 that Chapter ch:boschB quotes, with 30, 100 and 200 zeros
for p,tau in [(mp.mpf('0.5'),100),(mp.mpf('0.9'),gam[0]**(2/mp.mpf('0.9')))]:
    tau=mp.mpf(tau)
    P("p=%s tau=%s: xi'/xi %s | 30 zeros %s | 100 zeros %s | 200 zeros %s"%(mp.nstr(p,3),mp.nstr(tau,6),mp.nstr(u(p,tau),6),mp.nstr(u_zeros(p,tau,30),6),mp.nstr(u_zeros(p,tau,100),6),mp.nstr(u_zeros(p,tau,200),6)))
# line 39: positivity of u_p on a grid
mn=mp.inf; cnt=0; neg=0
for p in ['0.5','0.9','0.99','0.999']:
    p=mp.mpf(p)
    for k in range(-12,25):
        tau=mp.mpf(10)**(mp.mpf(k)/3); v=u(p,tau); cnt+=1
        if v<=0: neg+=1
        if v<mn: mn=v; at=(p,tau)
P("positivity grid: %d points, p in {0.5,0.9,0.99,0.999}, tau in [1e-4,1e8]: min u_p = %s at p=%s tau=%s"%(cnt,mp.nstr(mn,4),mp.nstr(at[0],3),mp.nstr(at[1],3)))
print("(not in the archived line) points with u_p <= 0:",neg,flush=True)
open("lines37_39_output.txt","w").write("\n".join(out)+"\n")
