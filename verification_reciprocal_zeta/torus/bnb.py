"""Branch and bound over the torus for the extremes of Theorem rz:thm:torus, in numpy floating point.

Quantities, over w in T^{pi(n)} and x in a band [a, b]:
  kind "m": prove |P_n(x,w)|^2 >= target^2       (certified lower bound target for m_n)
  kind "M": prove -Re(P_n'/P_n)(x,w) <= target  (certified upper bound target for M_n)
  kind "l": prove |P_n'/P_n|^2 <= target^2       (certified upper bound target for l_n)
Each box (centre c, half widths h) is bounded by the mean value form F(c) + sum_j sup_box|d_j F| h_j, with the
derivative enclosures computed in ball arithmetic (centre value plus a radius that bounds the variation of every
term k over the box: |beta_k k^{-2x} e^{-i phi_k} - centre| <= beta_k k^{-2x_c} [(k^{2h_x}-1) + k^{2h_x} min(2, dphi_k)]),
and by the first order bound where it is better. A box is discarded when its bound clears the target by SLACK = 1e-9,
far above the floating point error of these sums; otherwise it is split in the coordinate j with the largest
sup|d_j F| h_j. If a box centre violates the target, the target is false and the run stops.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import sys, time, json
import numpy as np
from torus_common import Torus

SLACK = 1e-9


class Enc:
    def __init__(self, T):
        self.T = T
        self.L = T.logk
        self.beta = T.beta
        self.V = T.V
        self.K = T.K

    def balls(self, C, H):
        L, V = self.L, self.V
        xc, hx = C[:, 0], H[:, 0]
        wc, hw = C[:, 1:], H[:, 1:]
        Ec = self.beta * np.exp(-2.0 * np.outer(xc, L))
        phi = wc @ V.T
        dphi = hw @ V.T
        gx = np.exp(2.0 * np.outer(hx, L))
        nu = (gx - 1.0) + gx * np.minimum(2.0, dphi)
        t = Ec * np.exp(-1j * phi)
        Ab = Ec * nu
        B = C.shape[0]
        K1 = 1 + self.K
        P_c, P_r = t.sum(1), Ab.sum(1)
        t2, A2 = 2 * L * t, 2 * L * Ab
        N_c, N_r = t2.sum(1), A2.sum(1)
        dP_c = np.empty((B, K1), complex)
        dP_r = np.empty((B, K1))
        dN_c = np.empty((B, K1), complex)
        dN_r = np.empty((B, K1))
        dP_c[:, 0], dP_r[:, 0] = -N_c, N_r
        dN_c[:, 0] = -(4 * L * L * t).sum(1)
        dN_r[:, 0] = (4 * L * L * Ab).sum(1)
        dP_c[:, 1:] = -1j * (t @ V)
        dP_r[:, 1:] = Ab @ V
        dN_c[:, 1:] = -1j * (t2 @ V)
        dN_r[:, 1:] = A2 @ V
        return P_c, P_r, N_c, N_r, dP_c, dP_r, dN_c, dN_r

    def bound(self, kind, C, H):
        """Return (value at centre, bound over box, per coordinate slack contributions).
        For "m" the bound is a lower bound of |P|^2, for "M" an upper bound of Re(N/P), for "l" of |N/P|^2."""
        P_c, P_r, N_c, N_r, dP_c, dP_r, dN_c, dN_r = self.balls(C, H)
        aP = np.abs(P_c)
        if kind == "m":
            F = aP ** 2
            cc = np.conj(P_c)[:, None] * dP_c
            rad = aP[:, None] * dP_r + np.abs(dP_c) * P_r[:, None] + P_r[:, None] * dP_r
            g = 2 * (np.abs(cc.real) + rad)
            contrib = g * H
            lb = F - contrib.sum(1)
            fo = np.where(aP > P_r, (aP - P_r) ** 2, 0.0)
            return F, np.maximum(lb, fo), contrib
        # ratio R = N/P and its derivatives
        ok = aP > P_r
        R_c = N_c / P_c
        c2 = P_c ** 2
        a2 = np.abs(c2)
        r2 = 2 * aP * P_r + P_r ** 2
        num_c = dN_c * P_c[:, None] - N_c[:, None] * dP_c
        num_r = (np.abs(dN_c) * P_r[:, None] + aP[:, None] * dN_r + dN_r * P_r[:, None]
                 + np.abs(N_c)[:, None] * dP_r + np.abs(dP_c) * N_r[:, None] + N_r[:, None] * dP_r)
        good = ok & (a2 > r2)
        with np.errstate(divide="ignore", invalid="ignore"):
            z_c = num_c / c2[:, None]
            z_r = (num_r * a2[:, None] + np.abs(num_c) * r2[:, None]) / (a2 * (a2 - r2))[:, None]
            R_r = (N_r * aP + np.abs(N_c) * P_r) / (aP * (aP - P_r))
        if kind == "M":
            F = R_c.real
            g = np.abs(z_c.real) + z_r
            contrib = g * H
            ub = F + contrib.sum(1)
            fo = R_c.real + R_r
            ub = np.minimum(ub, fo)
        else:
            F = np.abs(R_c) ** 2
            cz = np.conj(R_c)[:, None] * z_c
            rad = np.abs(R_c)[:, None] * z_r + np.abs(z_c) * R_r[:, None] + R_r[:, None] * z_r
            g = 2 * (np.abs(cz.real) + rad)
            contrib = g * H
            ub = F + contrib.sum(1)
            fo = (np.abs(R_c) + R_r) ** 2
            ub = np.minimum(ub, fo)
        ub = np.where(good, ub, np.inf)
        contrib = np.where(good[:, None], contrib, np.inf)
        return F, ub, contrib


def run(n, a, b, kind, target, batch=100000, max_boxes=10 ** 12, verbose=True, time_limit=None, report_every=30.0,
        leaf_sink=None):
    """Prove the target over the band [a,b] x torus, depth first (bounded memory).
    Returns a dict with counts and the extreme centre value seen. leaf_sink(C, H) receives certified leaves."""
    T = Torus(n)
    E = Enc(T)
    C0 = np.concatenate([[0.5 * (a + b)], np.full(T.K, np.pi)])[None, :]
    H0 = np.concatenate([[0.5 * (b - a)], np.full(T.K, np.pi)])[None, :]
    stack = [(C0, H0)]
    tgt = target ** 2 if kind in ("m", "l") else target
    processed = 0
    leaves = 0
    best = None
    best_pt = None
    t0 = time.time()
    last = t0
    maxdepth = 0
    while stack:
        C, H = stack.pop()
        if C.shape[0] > batch:
            stack.append((C[batch:], H[batch:]))
            C, H = C[:batch], H[:batch]
        F, bd, contrib = E.bound(kind, C, H)
        processed += C.shape[0]
        if kind == "m":
            j = int(np.argmin(F))
            if best is None or F[j] < best:
                best, best_pt = F[j], C[j].copy()
            if F[j] < tgt:
                return {"status": "violated", "n": n, "band": [a, b], "kind": kind, "target": target,
                        "at": C[j].tolist(), "value": float(np.sqrt(F[j])), "processed": processed}
            done = bd >= tgt + SLACK
        else:
            j = int(np.argmax(F))
            if best is None or F[j] > best:
                best, best_pt = F[j], C[j].copy()
            if F[j] > tgt:
                v = F[j]
                return {"status": "violated", "n": n, "band": [a, b], "kind": kind, "target": target,
                        "at": C[j].tolist(), "value": float(np.sqrt(v) if kind == "l" else v), "processed": processed}
            done = bd <= tgt - SLACK
        leaves += int(done.sum())
        if leaf_sink is not None and done.any():
            leaf_sink(C[done], H[done])
        rest = ~done
        if rest.any():
            Cr, Hr, cr = C[rest], H[rest], contrib[rest]
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
            C1 = Cr.copy()
            C2 = Cr.copy()
            C1[rows, jj] -= Hn[rows, jj]
            C2[rows, jj] += Hn[rows, jj]
            stack.append((np.concatenate([C1, C2]), np.concatenate([Hn, Hn])))
        now = time.time()
        if verbose and now - last > report_every:
            last = now
            npend = sum(c.shape[0] for c, _ in stack)
            bv = best if kind == "M" else np.sqrt(best)
            print(f"  n={n} [{a},{b}] {kind} target={target}: processed={processed:.3e} leaves={leaves:.3e} "
                  f"pending={npend} best={bv:.6f} t={now-t0:.0f}s", flush=True)
        if processed > max_boxes or (time_limit is not None and now - t0 > time_limit):
            return {"status": "gave up", "n": n, "band": [a, b], "kind": kind, "target": target,
                    "processed": processed, "leaves": leaves,
                    "best": float(best if kind == "M" else np.sqrt(best)), "seconds": now - t0}
    return {"status": "certified", "n": n, "band": [a, b], "kind": kind, "target": target, "processed": processed,
            "leaves": leaves, "extreme_centre_value": float(best if kind == "M" else np.sqrt(best)),
            "extreme_centre_point": best_pt.tolist(), "seconds": time.time() - t0}


if __name__ == "__main__":
    n = int(sys.argv[1]); a = float(sys.argv[2]); b = float(sys.argv[3]); kind = sys.argv[4]; target = float(sys.argv[5])
    r = run(n, a, b, kind, target)
    print(json.dumps(r))
