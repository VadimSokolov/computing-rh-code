import mpmath as mp
from flint import acb, arb
from polya_flint import setup, xi_pair
S = setup(420, 200, 3)
mp.mp.dps = 40
XI0 = xi_pair(S, acb(0))[0].real
import tilted
ST = tilted.setup(0.1, 320, 0.003, 4.5)
def to_acb(z): return acb(arb(mp.nstr(mp.re(z),45)), arb(mp.nstr(mp.im(z),45)))
def to_mp(z): return mp.mpc(mp.mpf(z.real.mid().str(45,radius=False)), mp.mpf(z.imag.mid().str(45,radius=False)))
def Lap(u):            # E exp(-u T) = xi(1/2)/xi(1/2+sqrt u), from the Polya density only
    s = mp.sqrt(mp.mpc(u))
    if abs(mp.im(s)) < 120: A,_ = xi_pair(S, to_acb(s))
    else: A,_ = tilted.xi_pair_direct(ST, to_acb(s))
    return to_mp(acb(XI0)/A)
def dens(t):  return mp.invertlaplace(Lap, t, method='talbot')
def surv(t):  return mp.invertlaplace(lambda u: (1-Lap(u))/u, t, method='talbot')
