"""Second condition of Theorem rz:thm:torus without logarithms: the function G, in numpy floating point.

With T = log(y/pi) - 1/y, a1 = |P(3/2-sigma,w)| + r(3/2-sigma)/y, a2 = |P(1+sigma,w)| - r(1+sigma)/y and
E(sigma) = exp(T (2 sigma - 1/2)), put
    G(sigma, w) = a2 E(sigma) - a1.
For 0.3 <= sigma <= 1/2 the denominator 2 sigma - 1/2 is positive and a1 > 0, so G > 0 at a point is equivalent to
a2 > 0 together with q = (log a1 - log a2)/(2 sigma - 1/2) < T. Hence G > 0 on the compact set [0.3, 1/2] x torus is
equivalent to the second condition log(y/pi) - 1/y > B_n(y) with |P(1+sigma,w)| > r(1+sigma)/y. G has no singularity
where P(3/2-sigma, w) vanishes, which is what made the search for q slow.

Lower bound of G over a box (the larger of the two):
  zero order:  inf(a2) E(sigma_lo) - sup(a1)   (with E(sigma_hi) if inf(a2) < 0), from the balls of P;
  mean value:  G(c) - sum_j D_j h_j, with D_j >= sup |d_j G| along coordinate paths in the box:
     w_p:   D = E(sigma_hi) sup|d_p |P(x2)|| + sup|d_p |P(x1)||,
     sigma: D = (sup|d_x |P(x2)|| + r'(x2)_max/y) E(sigma_hi) + sup|a2| 2T E(sigma_hi) + sup|d_x |P(x1)|| + r'(x1)_max/y,
  where sup|d |P|| comes from the balls of P and its derivatives (|P| is Lipschitz with constant sup|dP| even
  where P vanishes, and the sharper centre plus radius form is used where |P| stays away from 0).
"""
import numpy as np
from bnb_q import EncQ


def dabs(Pc, Pr, Dc, Dr):
    """Upper bound of sup |d|P|| over the box from the balls (Pc, Pr) of P and (Dc, Dr) of a partial derivative."""
    aP = np.abs(Pc)
    Alo = aP - Pr
    ok = Alo > 0
    with np.errstate(divide="ignore", invalid="ignore"):
        cen = np.real(np.conj(Pc)[:, None] * Dc) / aP[:, None]
        rad_num = aP[:, None] * Dr + np.abs(Dc) * Pr[:, None] + Pr[:, None] * Dr
        var = rad_num / Alo[:, None] + np.abs(Dc) * (aP / Alo - 1.0)[:, None]
        tight = np.abs(cen) + var
    crude = np.abs(Dc) + Dr
    return np.where(ok[:, None], np.minimum(tight, crude), crude)


class EncG(EncQ):
    def __init__(self, T, y):
        super().__init__(T, y)
        self.Tq = float(np.log(y / np.pi) - 1.0 / y)

    def bound(self, C, H):
        y, Tq = self.y, self.Tq
        sc, hs = C[:, 0], H[:, 0]
        wc, hw = C[:, 1:], H[:, 1:]
        x1c, x2c = 1.5 - sc, 1.0 + sc
        P1, P1r, P1x, P1xr, P1w, P1wr = self.P_ball(x1c, hs, wc, hw)
        P2, P2r, P2x, P2xr, P2w, P2wr = self.P_ball(x2c, hs, wc, hw)
        r1c, r1lo, r1hi, rp1lo, rp1hi = self.r_ball(x1c, hs)
        r2c, r2lo, r2hi, rp2lo, rp2hi = self.r_ball(x2c, hs)
        Ec = np.exp(Tq * (2 * sc - 0.5))
        Elo = np.exp(Tq * (2 * (sc - hs) - 0.5))
        Ehi = np.exp(Tq * (2 * (sc + hs) - 0.5))
        A1c, A2c = np.abs(P1), np.abs(P2)
        Gc = (A2c - r2c / y) * Ec - (A1c + r1c / y)
        a2lo = A2c - P2r - r2hi / y
        a1hi = A1c + P1r + r1hi / y
        G0 = np.where(a2lo >= 0, a2lo * Elo, a2lo * Ehi) - a1hi
        s1w = dabs(P1, P1r, P1w, P1wr)
        s2w = dabs(P2, P2r, P2w, P2wr)
        s1x = dabs(P1, P1r, P1x[:, None], P1xr[:, None])[:, 0]
        s2x = dabs(P2, P2r, P2x[:, None], P2xr[:, None])[:, 0]
        a2abs = A2c + P2r + r2hi / y
        Ds = (s2x + rp2hi / y) * Ehi + a2abs * 2 * Tq * Ehi + s1x + rp1hi / y
        Dw = Ehi[:, None] * s2w + s1w
        contrib = np.concatenate([(Ds * hs)[:, None], Dw * hw], axis=1)
        G1 = Gc - contrib.sum(1)
        return Gc, np.maximum(G0, G1), contrib
