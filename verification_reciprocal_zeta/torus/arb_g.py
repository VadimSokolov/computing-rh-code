"""Arb version of bnb_g.EncG: a rigorous lower bound of G(sigma, w) = a2 E(sigma) - a1 over a box (see bnb_g.py).

Same formulas as the floating point code, in ball arithmetic at 128 bits (arb_disc.ArbDisc supplies the balls of P
and of its partial derivatives, r(x) and r'(x) over the box, and the bound dabs of sup|d|P||); T = log(y/pi) - 1/y is
the exact value. The box may be given by floats or by arb balls (centres and half widths)."""
from flint import arb
from arb_disc import ArbDisc, exact


def _lo(a):
    return a.lower()


def _min_lo(a, b):
    a, b = a.lower(), b.lower()
    return a if a < b else b


class ArbG(ArbDisc):
    def bound_g(self, c, h, y):
        yy = arb(y)
        Tq = (yy / arb.pi()).log() - 1 / yy
        half = arb(1) / 2
        sc, hs = exact(c[0]), exact(h[0])
        wc = [exact(t) for t in c[1:]]
        hw = [exact(t) for t in h[1:]]
        x1c = arb(3) / 2 - sc
        x2c = 1 + sc
        P1, P1r, _, _, dP1, dP1r, _, _ = self.balls(x1c, hs, wc, hw)
        P2, P2r, _, _, dP2, dP2r, _, _ = self.balls(x2c, hs, wc, hw)
        r1c, r1lo, r1hi, rp1hi = self.r_ball(x1c, hs)
        r2c, r2lo, r2hi, rp2hi = self.r_ball(x2c, hs)
        Ec = (Tq * (2 * sc - half)).exp()
        Elo = (Tq * (2 * (sc - hs) - half)).exp()
        Ehi = (Tq * (2 * (sc + hs) - half)).exp()
        A1c, A2c = abs(P1), abs(P2)
        Gc = (A2c - r2c / yy) * Ec - (A1c + r1c / yy)
        a2lo = A2c - P2r - r2hi / yy
        a1hi = A1c + P1r + r1hi / yy
        G0 = (_min_lo(a2lo * Elo, a2lo * Ehi) - a1hi).lower()
        s1x = self.dabs(P1, P1r, dP1[0], dP1r[0])
        s2x = self.dabs(P2, P2r, dP2[0], dP2r[0])
        a2abs = A2c + P2r + r2hi / yy
        tot = ((s2x + rp2hi / yy) * Ehi + a2abs * 2 * Tq * Ehi + s1x + rp1hi / yy) * hs
        for j in range(self.K):
            s1w = self.dabs(P1, P1r, dP1[j + 1], dP1r[j + 1])
            s2w = self.dabs(P2, P2r, dP2[j + 1], dP2r[j + 1])
            tot += (Ehi * s2w + s1w) * hw[j]
        G1 = (Gc - tot).lower()
        return G0 if G0 > G1 else G1

    def certify_g(self, c, h, y):
        b = self.bound_g(c, h, y)
        return bool(b > 0), float(b.mid())
