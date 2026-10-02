# Keiper Li coefficients from the Polya density: Taylor coefficients of log xi(1+w) by a Cauchy integral
import mpmath as mp, json
from flint import acb, arb, ctx
from polya_flint import setup, xi_pair
S=setup(600,200,3); mp.mp.dps=150
N=512; r=mp.mpf(6)
vals=[]
for j in range(N):
    w=r*mp.expjpi(mp.mpf(2*j)/N)
    s=mp.mpf(0.5)+w                     # xi(1/2+s) with s=1/2+w gives xi(1+w)
    A,_=xi_pair(S,acb(arb(mp.nstr(mp.re(s),140)),arb(mp.nstr(mp.im(s),140))))
    z=mp.mpc(mp.mpf(A.real.mid().str(140,radius=False)),mp.mpf(A.imag.mid().str(140,radius=False)))
    vals.append(mp.log(z))
# unwrap is unnecessary: log xi analytic in |w|<14, principal branch continuous since xi(1+w) stays near positive axis? check
a=[mp.fsum(vals[j]*mp.expjpi(-mp.mpf(2*j*k)/N) for j in range(N))/N/r**k for k in range(0,121)]
lam=[]
for n in range(1,101):
    lam.append(n*mp.fsum(mp.binomial(n-1,j)*a[n-j] for j in range(n)))
exact1=1+mp.euler/2-mp.log(4*mp.pi)/2
print('lambda1',mp.nstr(mp.re(lam[0]),25),'exact',mp.nstr(exact1,25))
for n in [1,2,3,5,10,20,50,100]:
    asym=n/2*(mp.log(n)-1-mp.log(2*mp.pi)+mp.euler)
    print(n, mp.nstr(mp.re(lam[n-1]),20), mp.nstr(mp.im(lam[n-1]),3), 'asym',mp.nstr(asym,10))
json.dump(dict(lam=[mp.nstr(mp.re(x),30) for x in lam]),open('li.json','w'))
