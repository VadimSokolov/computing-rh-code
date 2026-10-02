"""Tilted Polya contour: xi(1/2+s) = int_{Im u = eta} Phi(u) e^{su} du, eta = pi/4 - delta.
   Phi is analytic in |Im u| < pi/4; on the tilted line |e^{su}| carries e^{-theta*eta},
   so the cancellation drops from e^{-pi theta/4} to e^{-delta theta}."""
from flint import arb, acb, acb_poly, ctx
import math
def Phi_c(u, tol_bits):
    # Phi(u) for Re u >= 0 (use evenness otherwise)
    if u.real < 0: u = -u
    pi = arb.pi(); e2 = (2*u).exp(); s = acb(0); n = 1
    e92 = (9*u/2).exp(); e52 = (5*u/2).exp()
    while True:
        t = (2*pi**2*n**4*e92 - 3*pi*n**2*e52)*(-pi*n*n*e2).exp()
        s += t
        if n > 2 and abs(complex(t)) < 2.0**(-tol_bits): break
        n += 1
    return 2*s
def setup(delta, prec_bits, h, xmax):
    ctx.prec = prec_bits
    eta = arb.pi()/4 - arb(delta)
    J = int(math.ceil(xmax/h)); hh = arb(h)
    c0 = []; c1 = []
    for j in range(-J, J+1):
        u = acb(j*hh, eta)
        p = Phi_c(u, prec_bits+20)
        c0.append(hh*p); c1.append(hh*u*p)
    return dict(J=J, h=hh, eta=eta, P=acb_poly(c0), Q=acb_poly(c1))
def xi_pair(S, s):
    s = acb(s)
    E = (s*S['h']).exp(); pref = (s*acb(0, S['eta'])).exp() * E**(-S['J'])
    return pref*S['P'](E), pref*S['Q'](E)
def xi_pair_direct(S, s):
    s = acb(s); J=S['J']; h=S['h']
    cP = S['P'].coeffs(); cQ = S['Q'].coeffs()
    A = acb(0); B = acb(0)
    for k in range(len(cP)):
        e = (s*((k-J)*h)).exp()
        A += cP[k]*e; B += cQ[k]*e
    pref = (s*acb(0, S['eta'])).exp()
    return pref*A, pref*B
