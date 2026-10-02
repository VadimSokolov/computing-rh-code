"""Branch and bound for the second condition of Theorem rz:thm:torus, in numpy floating point.

Proves  q(sigma, w) = [log(|P(3/2-sigma,w)| + r(3/2-sigma)/y) - log(|P(1+sigma,w)| - r(1+sigma)/y)] / (2 sigma - 1/2)
        < target   for all 0.3 <= sigma <= 1/2 and w in the torus,
so that B_n(y) <= target. Mean value form in (sigma, w), with ball enclosures of P and its derivatives over each box
at x1 = 3/2 - sigma and x2 = 1 + sigma, as in bnb.py. A box whose lower bound of |P(x2,w)| - r(x2)/y is not positive
is split (or, if its centre value is not positive, reported: then B_n(y) is infinite). SLACK = 1e-9.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, time, json, math
import numpy as np
from torus_common import Torus

SLACK = 1e-9


class EncQ:
    def __init__(self, T, y):
        self.T = T
        self.L = T.logk
        self.beta = T.beta
        self.absA = np.abs(T.alpha)
        self.V = T.V
        self.K = T.K
        self.y = y

    def P_ball(self, xc, hx, wc, hw):
        """P and grad (x then w) at centre, with radii over the box; x-derivative is d/dx."""
        L, V = self.L, self.V
        Ec = self.beta * np.exp(-2.0 * np.outer(xc, L))
        phi = wc @ V.T
        dphi = hw @ V.T
        gx = np.exp(2.0 * np.outer(hx, L))
        nu = (gx - 1.0) + gx * np.minimum(2.0, dphi)
        t = Ec * np.exp(-1j * phi)
        Ab = Ec * nu
        P_c, P_r = t.sum(1), Ab.sum(1)
        Px_c, Px_r = (-2 * L * t).sum(1), (2 * L * Ab).sum(1)
        Pw_c, Pw_r = -1j * (t @ V), Ab @ V
        return P_c, P_r, Px_c, Px_r, Pw_c, Pw_r

    def r_ball(self, xc, hx):
        """r(x), r'(x) = -dr/dx: values at the centre and ranges over [xc-hx, xc+hx] (both decreasing in x)."""
        L = self.L
        e_c = np.exp(-2.0 * np.outer(xc, L))
        e_lo = np.exp(-2.0 * np.outer(xc + hx, L))
        e_hi = np.exp(-2.0 * np.outer(xc - hx, L))
        r_c = (self.absA * e_c).sum(1)
        r_lo, r_hi = (self.absA * e_lo).sum(1), (self.absA * e_hi).sum(1)
        rp_lo, rp_hi = (2 * L * self.absA * e_lo).sum(1), (2 * L * self.absA * e_hi).sum(1)
        return r_c, r_lo, r_hi, rp_lo, rp_hi

    def bound(self, C, H):
        """C[:,0] = sigma centre, C[:,1:] = w centre; H half widths. Returns centre value, upper bound, contributions,
        and the lower bound of |P(x2)| - r(x2)/y over the box (and its centre value)."""
        y = self.y
        sc, hs = C[:, 0], H[:, 0]
        wc, hw = C[:, 1:], H[:, 1:]
        x1c, x2c = 1.5 - sc, 1.0 + sc
        P1, P1r, P1x, P1xr, P1w, P1wr = self.P_ball(x1c, hs, wc, hw)
        P2, P2r, P2x, P2xr, P2w, P2wr = self.P_ball(x2c, hs, wc, hw)
        r1c, r1lo, r1hi, rp1lo, rp1hi = self.r_ball(x1c, hs)
        r2c, r2lo, r2hi, rp2lo, rp2hi = self.r_ball(x2c, hs)
        A1c, A2c = np.abs(P1), np.abs(P2)
        a1c = A1c + r1c / y
        a2c = A2c - r2c / y
        den_c = 2 * sc - 0.5
        # lower bounds over the box
        A1lo = np.maximum(A1c - P1r, 0.0)
        A2lo = A2c - P2r
        a1lo = A1lo + r1lo / y
        a2lo = A2lo - r2hi / y
        den_lo = 2 * (sc - hs) - 0.5
        valid_c = a2c > 0
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            numc = np.log(a1c) - np.log(a2c)
            q_c = numc / den_c
            # |d|P|/dw| <= |P_w| (Cauchy Schwarz), and over the box |Re(conj(P) P_w)/|P|| <= |P_w|; use the sharper
            # centre-plus-radius form when |P| stays away from 0
            def dabs(Pc, Pr, Dc, Dr):
                Alo = np.abs(Pc) - Pr
                ok = Alo > 0
                cen = np.real(np.conj(Pc)[:, None] * Dc) / np.abs(Pc)[:, None]
                # variation of conj(P)D/|P| over the box: bound by the ball of conj(P)D divided by |P| range
                rad_num = np.abs(Pc)[:, None] * Dr + np.abs(Dc) * Pr[:, None] + Pr[:, None] * Dr
                # |conj(P)D/|P| - conj(Pc)Dc/|Pc|| <= rad_num/Alo + |Dc| |Pc| * (1/Alo - 1/|Pc|)
                var = rad_num / Alo[:, None] + np.abs(Dc) * (np.abs(Pc) / Alo - 1.0)[:, None]
                sup = np.where(ok[:, None], np.minimum(np.abs(cen) + var, np.abs(Dc) + Dr), np.abs(Dc) + Dr)
                return sup
            s1w = dabs(P1, P1r, P1w, P1wr)            # sup |d|P(x1)|/dw_p|
            s2w = dabs(P2, P2r, P2w, P2wr)
            s1x = dabs(P1, P1r, P1x[:, None], P1xr[:, None])[:, 0]  # sup |d|P|/dx| at x1
            s2x = dabs(P2, P2r, P2x[:, None], P2xr[:, None])[:, 0]
            ok = (a2lo > 0) & (den_lo > 0)
            # d q / d w_p = (d|P1|/a1 - d|P2|/a2)/den
            gw = (s1w / a1lo[:, None] + s2w / a2lo[:, None]) / den_lo[:, None]
            # d q / d sigma = [ (-d_x|P1| + r'(x1)/y)/a1 - (d_x|P2| + r'(x2)/y)/a2 ]/den - 2 num/den^2
            num_hi = np.log(A1c + P1r + r1hi / y) - np.log(a2lo)
            num_lo = np.log(a1lo) - np.log(A2c + P2r - r2lo / y)
            num_abs = np.maximum(np.abs(num_hi), np.abs(num_lo))
            gs = ((s1x + rp1hi / y) / a1lo + (s2x + rp2hi / y) / a2lo) / den_lo + 2 * num_abs / den_lo ** 2
            contrib = np.concatenate([(gs * hs)[:, None], gw * hw], axis=1)
            ub = q_c + contrib.sum(1)
        ub = np.where(ok & valid_c, ub, np.inf)
        contrib = np.where((ok & valid_c)[:, None], contrib, np.inf)
        q_c = np.where(valid_c, q_c, np.inf)
        return q_c, ub, contrib, a2lo, a2c


def run_q(n, y, target, batch=200000, max_boxes=10 ** 11, verbose=True):
    T = Torus(n)
    E = EncQ(T, y)
    C0 = np.concatenate([[0.4], np.full(T.K, np.pi)])[None, :]
    H0 = np.concatenate([[0.1], np.full(T.K, np.pi)])[None, :]
    stackC, stackH = [C0], [H0]
    processed = leaves = 0
    best, best_pt = -np.inf, None
    t0 = time.time()
    while stackC:
        C = np.concatenate(stackC)
        H = np.concatenate(stackH)
        newC, newH = [], []
        for i in range(0, C.shape[0], batch):
            Cb, Hb = C[i:i + batch], H[i:i + batch]
            qc, ub, contrib, a2lo, a2c = E.bound(Cb, Hb)
            processed += Cb.shape[0]
            j = np.argmax(qc)
            if qc[j] > best:
                best, best_pt = qc[j], Cb[j].copy()
            if not np.all(a2c > 0):
                return {"status": "B infinite", "n": n, "y": y, "at": Cb[np.argmin(a2c)].tolist(), "processed": processed}
            if np.any(qc >= target):
                return {"status": "violated", "n": n, "y": y, "target": target, "value": float(qc.max()),
                        "at": Cb[np.argmax(qc)].tolist(), "processed": processed}
            done = ub <= target - SLACK
            leaves += int(done.sum())
            rest = ~done
            if rest.any():
                Cr, Hr, cr = Cb[rest], Hb[rest], contrib[rest]
                inf_rows = ~np.isfinite(cr).all(1)
                cr = np.where(np.isfinite(cr), cr, 0.0)
                jj = np.argmax(cr, axis=1)
                if inf_rows.any():
                    sc = Hr.copy()
                    sc[:, 0] *= 30.0
                    jj = np.where(inf_rows, np.argmax(sc, axis=1), jj)
                rows = np.arange(Cr.shape[0])
                Hn = Hr.copy()
                Hn[rows, jj] *= 0.5
                C1, C2 = Cr.copy(), Cr.copy()
                C1[rows, jj] -= Hn[rows, jj]
                C2[rows, jj] += Hn[rows, jj]
                newC += [C1, C2]
                newH += [Hn, Hn.copy()]
        stackC, stackH = newC, newH
        if newC and verbose:
            npend = sum(c.shape[0] for c in newC)
            print(f"  q n={n} y={y} target={target:.6f}: processed={processed} leaves={leaves} pending={npend} best={best:.6f} t={time.time()-t0:.1f}s", flush=True)
        if processed > max_boxes:
            return {"status": "gave up", "processed": processed, "best": float(best)}
    return {"status": "certified", "n": n, "y": y, "target": target, "processed": processed, "leaves": leaves,
            "max_centre_value": float(best), "at": best_pt.tolist(), "seconds": time.time() - t0}


if __name__ == "__main__":
    n = int(sys.argv[1]); y = float(sys.argv[2])
    target = math.log(y / math.pi) - 1 / y
    print(json.dumps(run_q(n, y, target)))
