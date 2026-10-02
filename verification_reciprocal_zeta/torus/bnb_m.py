"""Enclosures for the ratio kinds M and l that use a certified lower bound m_low of |P_n| on the band.

Same as bnb.Enc for kinds M and l, except that the lower bound of |P| over a box is A_lo = max(m_low, |P_c| - P_r)
instead of |P_c| - P_r, which must then be positive. This is valid because |P_n(x, w)| >= m_low everywhere on the band
and the torus (certified by the kind m runs of the full branch and bound, with the same band). With it,
  |N/P - N_c/P_c| <= (N_r |P_c| + |N_c| P_r) / (A_lo |P_c|),
  |num/P^2 - num_c/P_c^2| <= (num_r |P_c|^2 + |num_c| r2) / (A_lo^2 |P_c|^2),  r2 = 2|P_c| P_r + P_r^2,
and no box is rejected for P_r being large compared with |P_c|.
"""
import numpy as np
from bnb import Enc


class EncM(Enc):
    def __init__(self, T, m_low):
        super().__init__(T)
        self.m_low = float(m_low)

    def bound(self, kind, C, H):
        assert kind in ("M", "l")
        P_c, P_r, N_c, N_r, dP_c, dP_r, dN_c, dN_r = self.balls(C, H)
        aP = np.abs(P_c)
        A_lo = np.maximum(self.m_low, aP - P_r)
        R_c = N_c / P_c
        c2 = P_c ** 2
        a2 = aP ** 2
        r2 = 2 * aP * P_r + P_r ** 2
        num_c = dN_c * P_c[:, None] - N_c[:, None] * dP_c
        num_r = (np.abs(dN_c) * P_r[:, None] + aP[:, None] * dN_r + dN_r * P_r[:, None]
                 + np.abs(N_c)[:, None] * dP_r + np.abs(dP_c) * N_r[:, None] + N_r[:, None] * dP_r)
        z_c = num_c / c2[:, None]
        z_r = (num_r * a2[:, None] + np.abs(num_c) * r2[:, None]) / (A_lo ** 2 * a2)[:, None]
        R_r = (N_r * aP + np.abs(N_c) * P_r) / (A_lo * aP)
        if kind == "M":
            F = R_c.real
            g = np.abs(z_c.real) + z_r
            contrib = g * H
            ub = np.minimum(F + contrib.sum(1), R_c.real + R_r)
        else:
            F = np.abs(R_c) ** 2
            cz = np.conj(R_c)[:, None] * z_c
            rad = np.abs(R_c)[:, None] * z_r + np.abs(z_c) * R_r[:, None] + R_r[:, None] * z_r
            g = 2 * (np.abs(cz.real) + rad)
            contrib = g * H
            ub = np.minimum(F + contrib.sum(1), (np.abs(R_c) + R_r) ** 2)
        return F, ub, contrib
