"""Polya route to xi(1/2+s) and xi'(1/2+s): trapezoid rule on R for
   xi(1/2+s) = int Phi(u) e^{su} du  (Phi even),  xi'(1/2+s) = int u Phi(u) e^{su} du,
   Phi(u) = 2 sum_n (2 pi^2 n^4 e^{9u/2} - 3 pi n^2 e^{5u/2}) exp(-pi n^2 e^{2u})  (Titchmarsh 10.1).
   Nodes u_j = j h, |j| <= J. Ball arithmetic (python-flint) for the sums."""
from flint import arb, acb, acb_poly, ctx
def setup(prec_bits=240, h_inv=200, umax=3):
    ctx.prec = prec_bits
    h = arb(1)/h_inv; J = int(umax*h_inv)
    pi = arb.pi()
    Phi=[]
    for j in range(J+1):
        u = j*h; e2=(2*u).exp(); s=arb(0); n=1
        while True:
            t = (2*pi**2*n**4*(9*u/2).exp() - 3*pi*n**2*(5*u/2).exp())*(-pi*n*n*e2).exp()
            s += t
            if n>2 and abs(float(t.mid())) < 2.0**(-prec_bits-20): break
            n += 1
        Phi.append(2*s)
    # P(z) = h*[Phi0/2 + sum_{j>=1} Phi_j z^j]; xi = P(E)+P(1/E)
    cA = [h*Phi[0]/2] + [h*Phi[j] for j in range(1,J+1)]
    cB = [arb(0)] + [h*j*h*Phi[j] for j in range(1,J+1)]
    return dict(h=h, PA=acb_poly(cA), PB=acb_poly(cB))
def xi_pair(S, s):
    s = acb(s) if not isinstance(s, acb) else s
    E = (s*S['h']).exp(); Ei = 1/E
    A = S['PA'](E) + S['PA'](Ei)
    B = S['PB'](E) - S['PB'](Ei)   # u Phi odd weight: int u Phi e^{su} over R
    return A, B
