"""Arb (python-flint) versions of the branch and bound enclosures of bnb.py (kinds m, M, l) and bnb_q.py (kind q).

The same disc enclosures as the floating point code, evaluated in ball arithmetic at 128 bits, so that every
rounding error is accounted for. Every quantity that the enclosure formulas treat as exact (centre values, radii,
lower bounds of denominators) is an Arb ball containing it or an exact endpoint on the safe side, and the final bound
is the safe endpoint of the resulting ball. The formulas themselves are the ones derived in the report (N3.md):
  term k over the box:  |beta_k k^{-2x} e^{-i phi_k} - centre| <= beta_k k^{-2x_c} [(k^{2h_x} - 1) + k^{2h_x} min(2, dphi_k)],
  mean value form F(box) <= F(c) + sum_j sup|d_j F| h_j (upper kinds), >= F(c) - sum_j ... (kind m),
  first order form from the ball of P, N/P.
Also: a naive Arb check for kind q (Arb interval evaluation of the partial derivatives over the box, no hand made
radius formulas), bisecting along the coordinate with the largest contribution, as check2 in arb_leafcheck.py.
"""
from fractions import Fraction as Fr
from flint import arb, acb, fmpq, ctx
from torus_common import exact_coeffs, primes_upto, vp

ctx.prec = 128


def exact(v):
    """Exact arb of a Python float (binary64 values are dyadic rationals); an arb is passed through unchanged, so the
    enclosures also accept boxes whose centres and half widths are balls (used by the replay, bnb7.py)."""
    if isinstance(v, arb):
        return v
    f = Fr(v)
    return arb(fmpq(f.numerator, f.denominator))


def dec(v):
    """Exact arb of the shortest decimal representation of a float target (e.g. 1.062 -> 1062/1000)."""
    f = Fr(repr(float(v)))
    return arb(fmpq(f.numerator, f.denominator))


def amax(a, b):
    return a if a > b else b


def amin_up(u1, u2):
    """Smaller of two exact upper endpoints."""
    return u1 if u1 < u2 else u2


class ArbDisc:
    def __init__(self, n):
        self.n = n
        ab = exact_coeffs(n)
        self.beta = [arb(fmpq(b.numerator, b.denominator)) for a, b in ab]
        self.absA = [abs(arb(fmpq(a.numerator, a.denominator))) for a, b in ab]
        self.ps = primes_upto(n)
        self.K = len(self.ps)
        self.V = [[vp(k, p) for p in self.ps] for k in range(1, n + 1)]
        self.L = [arb(k).log() for k in range(1, n + 1)]

    def balls(self, xc, hx, wc, hw):
        """Centre values (acb) and radii (arb) of P, N = -P', and their partial derivatives (x first, then w)."""
        K1 = 1 + self.K
        P_c, N_c = acb(0), acb(0)
        P_r, N_r = arb(0), arb(0)
        dP_c = [acb(0)] * K1
        dN_c = [acb(0)] * K1
        dP_r = [arb(0)] * K1
        dN_r = [arb(0)] * K1
        two = arb(2)
        for k in range(self.n):
            L = self.L[k]
            Vk = self.V[k]
            Ec = self.beta[k] * (-2 * xc * L).exp()
            phi = arb(0)
            dphi = arb(0)
            for j in range(self.K):
                if Vk[j]:
                    phi += Vk[j] * wc[j]
                    dphi += Vk[j] * hw[j]
            gx = (2 * hx * L).exp()
            mn = dphi if dphi < two else two
            nu = (gx - 1) + gx * mn
            t = Ec * acb(0, -phi).exp()
            Ab = Ec * nu
            t2 = 2 * L * t
            A2 = 2 * L * Ab
            P_c += t
            P_r += Ab
            N_c += t2
            N_r += A2
            dP_c[0] = dP_c[0] - t2
            dP_r[0] = dP_r[0] + A2
            dN_c[0] = dN_c[0] - 2 * L * t2
            dN_r[0] = dN_r[0] + 2 * L * A2
            for j in range(self.K):
                v = Vk[j]
                if v:
                    dP_c[j + 1] = dP_c[j + 1] + acb(0, -v) * t
                    dP_r[j + 1] = dP_r[j + 1] + v * Ab
                    dN_c[j + 1] = dN_c[j + 1] + acb(0, -v) * t2
                    dN_r[j + 1] = dN_r[j + 1] + v * A2
        return P_c, P_r, N_c, N_r, dP_c, dP_r, dN_c, dN_r

    def bound(self, kind, c, h):
        """Rigorous bound over the box (floats c, h): lower bound of |P|^2 (m), upper bound of Re(N/P) (M) or of
        |N/P|^2 (l). Returns an exact arb endpoint, or None when the enclosure does not apply (P may vanish)."""
        xc, hx = exact(c[0]), exact(h[0])
        wc = [exact(t) for t in c[1:]]
        hw = [exact(t) for t in h[1:]]
        hs = [hx] + hw
        P_c, P_r, N_c, N_r, dP_c, dP_r, dN_c, dN_r = self.balls(xc, hx, wc, hw)
        aP = abs(P_c)
        J = range(len(hs))
        if kind == "m":
            F = aP * aP
            s = arb(0)
            for j in J:
                cc = P_c.conjugate() * dP_c[j]
                rad = aP * dP_r[j] + abs(dP_c[j]) * P_r + P_r * dP_r[j]
                s += 2 * (abs(cc.real) + rad) * hs[j]
            lb = (F - s).lower()
            d = aP - P_r
            if d > 0:
                lb = amax(lb, (d * d).lower())
            return lb
        d = aP - P_r
        a2 = aP * aP
        r2 = 2 * aP * P_r + P_r * P_r
        if not (d > 0 and a2 - r2 > 0):
            return None
        c2 = P_c * P_c
        R_c = N_c / P_c
        R_r = (N_r * aP + abs(N_c) * P_r) / (aP * d)
        z_c, z_r = [], []
        for j in J:
            num_c = dN_c[j] * P_c - N_c * dP_c[j]
            num_r = (abs(dN_c[j]) * P_r + aP * dN_r[j] + dN_r[j] * P_r
                     + abs(N_c) * dP_r[j] + abs(dP_c[j]) * N_r + N_r * dP_r[j])
            z_c.append(num_c / c2)
            z_r.append((num_r * a2 + abs(num_c) * r2) / (a2 * (a2 - r2)))
        if kind == "M":
            F = R_c.real
            s = arb(0)
            for j in J:
                s += (abs(z_c[j].real) + z_r[j]) * hs[j]
            return amin_up((F + s).upper(), (R_c.real + R_r).upper())
        aR = abs(R_c)
        F = aR * aR
        s = arb(0)
        for j in J:
            cz = R_c.conjugate() * z_c[j]
            rad = aR * z_r[j] + abs(z_c[j]) * R_r + R_r * z_r[j]
            s += 2 * (abs(cz.real) + rad) * hs[j]
        fo = aR + R_r
        return amin_up((F + s).upper(), (fo * fo).upper())

    def certify(self, kind, target, c, h, y=None):
        """True when the Arb enclosure proves the target on the box (exact decimal target)."""
        if kind == "q":
            ub = self.bound_q(c, h, y)
            if ub is None:
                return False, None
            yy = arb(y)
            rhs = (yy / arb.pi()).log() - 1 / yy
            return bool(ub < rhs), float(ub.mid())
        b = self.bound(kind, c, h)
        if b is None:
            return False, None
        t = dec(target)
        if kind == "m":
            return bool(b > t * t), float(b.mid())
        if kind == "M":
            return bool(b < t), float(b.mid())
        return bool(b < t * t), float(b.mid())

    # ---- kind q (second condition) ----
    def r_ball(self, xc, hx):
        r_c = r_lo = r_hi = rp_hi = arb(0)
        for k in range(self.n):
            L, A = self.L[k], self.absA[k]
            r_c += A * (-2 * xc * L).exp()
            e_lo = (-2 * (xc + hx) * L).exp()
            e_hi = (-2 * (xc - hx) * L).exp()
            r_lo += A * e_lo
            r_hi += A * e_hi
            rp_hi += 2 * L * A * e_hi
        return r_c, r_lo, r_hi, rp_hi

    @staticmethod
    def dabs(Pc, Pr, Dc, Dr):
        """Upper endpoint of sup over the box of |d|P|| (or of |dP| where P may vanish)."""
        alt = (abs(Dc) + Dr).upper()
        aPc = abs(Pc)
        Alo = aPc - Pr
        if Alo > 0:
            cen = (Pc.conjugate() * Dc).real / aPc
            rad_num = aPc * Dr + abs(Dc) * Pr + Pr * Dr
            var = rad_num / Alo + abs(Dc) * (aPc / Alo - 1)
            return amin_up((abs(cen) + var).upper(), alt)
        return alt

    def bound_q(self, c, h, y):
        """Upper endpoint of q over the box (c[0] = sigma centre), or None if the enclosure does not apply."""
        yy = arb(y)
        sc, hs_ = exact(c[0]), exact(h[0])
        wc = [exact(t) for t in c[1:]]
        hw = [exact(t) for t in h[1:]]
        x1c = arb(3) / 2 - sc
        x2c = 1 + sc
        P1, P1r, N1, N1r, dP1, dP1r, _, _ = self.balls(x1c, hs_, wc, hw)
        P2, P2r, N2, N2r, dP2, dP2r, _, _ = self.balls(x2c, hs_, wc, hw)
        r1c, r1lo, r1hi, rp1hi = self.r_ball(x1c, hs_)
        r2c, r2lo, r2hi, rp2hi = self.r_ball(x2c, hs_)
        A1c, A2c = abs(P1), abs(P2)
        a1c = A1c + r1c / yy
        a2c = A2c - r2c / yy
        den_c = 2 * sc - arb(1) / 2
        if not (a2c > 0 and den_c > 0):
            return None
        q_c = (a1c.log() - a2c.log()) / den_c
        A1lo = amax((A1c - P1r).lower(), arb(0))
        a1lo = (A1lo + r1lo / yy).lower()
        a2lo = (A2c - P2r - r2hi / yy).lower()
        den_lo = (2 * (sc - hs_) - arb(1) / 2).lower()
        if not (a2lo > 0 and den_lo > 0):
            return None
        a1hi = (A1c + P1r + r1hi / yy).upper()
        a2hi = (A2c + P2r - r2lo / yy).upper()
        num_hi = (a1hi.log() - a2lo.log()).upper()
        num_lo = (a1lo.log() - a2hi.log()).lower()
        num_abs = amax(abs(num_hi), abs(num_lo))
        # x derivative of P is -N (radius N_r); w derivatives dP[1:]
        s1x = self.dabs(P1, P1r, -N1, N1r)
        s2x = self.dabs(P2, P2r, -N2, N2r)
        tot = arb(0)
        gs = ((s1x + rp1hi / yy) / a1lo + (s2x + rp2hi / yy) / a2lo) / den_lo + 2 * num_abs / (den_lo * den_lo)
        tot += gs * hs_
        for j in range(self.K):
            s1w = self.dabs(P1, P1r, dP1[j + 1], dP1r[j + 1])
            s2w = self.dabs(P2, P2r, dP2[j + 1], dP2r[j + 1])
            tot += (s1w / a1lo + s2w / a2lo) / den_lo * hw[j]
        return (q_c + tot).upper()


# ---- naive Arb check for kind q, independent of the hand made radius formulas ----
class ArbQNaive:
    def __init__(self, n):
        from arb_leafcheck import ArbTorus
        self.AT = ArbTorus(n)
        ab = exact_coeffs(n)
        self.absA = [abs(arb(fmpq(a.numerator, a.denominator))) for a, b in ab]
        self.L = [arb(k).log() for k in range(1, n + 1)]
        self.n = n
        self.K = self.AT.K

    def rr(self, x):
        r = rp = arb(0)
        for k in range(self.n):
            e = self.absA[k] * (-2 * x * self.L[k]).exp()
            r += e
            rp += 2 * self.L[k] * e
        return r, rp

    def bound_contrib(self, c, h, y):
        yy = arb(y)
        sc = exact(c[0])
        wc = [exact(t) for t in c[1:]]
        P1, _, _, _ = self.AT.eval(arb(3) / 2 - sc, wc, grad=False)
        P2, _, _, _ = self.AT.eval(1 + sc, wc, grad=False)
        r1, _ = self.rr(arb(3) / 2 - sc)
        r2, _ = self.rr(1 + sc)
        a1c = abs(P1) + r1 / yy
        a2c = abs(P2) - r2 / yy
        if not (a2c > 0):
            return None, None
        q_c = (a1c.log() - a2c.log()) / (2 * sc - arb(1) / 2)
        sb = arb(c[0], h[0])  # exact centre and radius (floats)
        wb = [arb(t, r) for t, r in zip(c[1:], h[1:])]
        x1b, x2b = arb(3) / 2 - sb, 1 + sb
        P1b, N1b, dP1b, _ = self.AT.eval(x1b, wb, grad=True)
        P2b, N2b, dP2b, _ = self.AT.eval(x2b, wb, grad=True)
        r1b, rp1b = self.rr(x1b)
        r2b, rp2b = self.rr(x2b)
        A1b, A2b = abs(P1b), abs(P2b)
        a1b = A1b + r1b / yy
        a2b = A2b - r2b / yy
        denb = 2 * sb - arb(1) / 2
        if not (a2b > 0 and denb > 0 and a1b > 0):
            return None, None
        numb = a1b.log() - a2b.log()

        def dabs_box(Pb, Db):
            if abs(Pb) > 0:
                return abs((Pb.conjugate() * Db).real / abs(Pb)).upper()
            return abs(Db).upper()
        contrib = []
        # sigma: d a1/d sigma = -d_x|P1| + r'(x1)/y, d a2/d sigma = d_x|P2| + r'(x2)/y (x derivative of P is -N)
        T1 = (dabs_box(P1b, -N1b) + rp1b / yy) / a1b
        T2 = (dabs_box(P2b, -N2b) + rp2b / yy) / a2b
        gs = (abs(T1) + abs(T2)) / denb + 2 * abs(numb) / (denb * denb)
        contrib.append((gs * exact(h[0])).upper())
        for j in range(self.K):
            g = (dabs_box(P1b, dP1b[j + 1]) / a1b + dabs_box(P2b, dP2b[j + 1]) / a2b) / denb
            contrib.append((abs(g) * exact(h[j + 1])).upper())
        s = sum(contrib, arb(0))
        return (q_c + s).upper(), [float(t) for t in contrib]

    def check2(self, c, h, y, budget=2000):
        yy = arb(y)
        rhs = (yy / arb.pi()).log() - 1 / yy
        stack = [(list(c), list(h))]
        used = 0
        while stack:
            cc, hh = stack.pop()
            used += 1
            if used > budget:
                return False, used
            b, contrib = self.bound_contrib(cc, hh, y)
            if b is not None and b < rhs:
                continue
            if contrib is not None:
                j = max(range(len(contrib)), key=lambda i: contrib[i])
            else:
                sc_ = [hh[0] * 30.0] + list(hh[1:])
                j = max(range(len(sc_)), key=lambda i: sc_[i])
            h2 = list(hh)
            h2[j] = hh[j] / 2
            c1, c2 = list(cc), list(cc)
            c1[j] -= h2[j]
            c2[j] += h2[j]
            stack.append((c2, h2))
            stack.append((c1, list(h2)))
        return True, used


def _bound_ratio_mlow(self, kind, m_low, c, h):
    """Arb version of bnb_m.EncM.bound (kinds M, l) with the certified lower bound m_low of |P| on the band."""
    xc, hx = exact(c[0]), exact(h[0])
    wc = [exact(t) for t in c[1:]]
    hw = [exact(t) for t in h[1:]]
    hs = [hx] + hw
    P_c, P_r, N_c, N_r, dP_c, dP_r, dN_c, dN_r = self.balls(xc, hx, wc, hw)
    aP = abs(P_c)
    ml = dec(m_low)
    d = (aP - P_r).lower()
    A_lo = d if d > ml else ml
    a2 = aP * aP
    r2 = 2 * aP * P_r + P_r * P_r
    c2 = P_c * P_c
    R_c = N_c / P_c
    R_r = (N_r * aP + abs(N_c) * P_r) / (A_lo * aP)
    z_c, z_r = [], []
    for j in range(len(hs)):
        num_c = dN_c[j] * P_c - N_c * dP_c[j]
        num_r = (abs(dN_c[j]) * P_r + aP * dN_r[j] + dN_r[j] * P_r
                 + abs(N_c) * dP_r[j] + abs(dP_c[j]) * N_r + N_r * dP_r[j])
        z_c.append(num_c / c2)
        z_r.append((num_r * a2 + abs(num_c) * r2) / (A_lo * A_lo * a2))
    if kind == "M":
        s = arb(0)
        for j in range(len(hs)):
            s += (abs(z_c[j].real) + z_r[j]) * hs[j]
        return amin_up((R_c.real + s).upper(), (R_c.real + R_r).upper())
    aR = abs(R_c)
    s = arb(0)
    for j in range(len(hs)):
        cz = R_c.conjugate() * z_c[j]
        s += 2 * (abs(cz.real) + aR * z_r[j] + abs(z_c[j]) * R_r + R_r * z_r[j]) * hs[j]
    fo = aR + R_r
    return amin_up((aR * aR + s).upper(), (fo * fo).upper())


def _certify_mlow(self, kind, target, m_low, c, h):
    b = self.bound_ratio_mlow(kind, m_low, c, h)
    t = dec(target)
    if kind == "M":
        return bool(b < t), float(b.mid())
    return bool(b < t * t), float(b.mid())


ArbDisc.bound_ratio_mlow = _bound_ratio_mlow
ArbDisc.certify_mlow = _certify_mlow


def _bound_q2(self, c, h, y):
    """Arb version of bnb_q2.EncQ2: the smaller of bound_q (mean value form) and the first order bound."""
    b1 = self.bound_q(c, h, y)
    yy = arb(y)
    sc, hs_ = exact(c[0]), exact(h[0])
    wc = [exact(t) for t in c[1:]]
    hw = [exact(t) for t in h[1:]]
    x1c = arb(3) / 2 - sc
    x2c = 1 + sc
    P1, P1r = self.balls(x1c, hs_, wc, hw)[:2]
    P2, P2r = self.balls(x2c, hs_, wc, hw)[:2]
    _, _, r1hi, _ = self.r_ball(x1c, hs_)
    _, _, r2hi, _ = self.r_ball(x2c, hs_)
    a1hi = (abs(P1) + P1r + r1hi / yy).upper()
    a2lo = (abs(P2) - P2r - r2hi / yy).lower()
    den_lo = (2 * (sc - hs_) - arb(1) / 2).lower()
    den_hi = (2 * (sc + hs_) - arb(1) / 2).upper()
    b2 = None
    if a2lo > 0 and den_lo > 0:
        num_hi = (a1hi.log() - a2lo.log()).upper()
        b2 = (num_hi / den_lo).upper() if num_hi >= 0 else (num_hi / den_hi).upper()
    if b1 is None:
        return b2
    if b2 is None:
        return b1
    return amin_up(b1, b2)


def _certify_q2(self, c, h, y):
    ub = self.bound_q2(c, h, y)
    if ub is None:
        return False, None
    yy = arb(y)
    rhs = (yy / arb.pi()).log() - 1 / yy
    return bool(ub < rhs), float(ub.mid())


ArbDisc.bound_q2 = _bound_q2
ArbDisc.certify_q2 = _certify_q2


def _naive_bound_contrib_fo(self, c, h, y):
    """ArbQNaive.bound_contrib with the first order bound from the interval evaluation of the numerator over the box."""
    b, contrib = self.bound_contrib(c, h, y)
    yy = arb(y)
    sb = arb(c[0], h[0])
    wb = [arb(t, r) for t, r in zip(c[1:], h[1:])]
    x1b, x2b = arb(3) / 2 - sb, 1 + sb
    P1b = self.AT.eval(x1b, wb, grad=False)[0]
    P2b = self.AT.eval(x2b, wb, grad=False)[0]
    r1b, _ = self.rr(x1b)
    r2b, _ = self.rr(x2b)
    a1b = abs(P1b) + r1b / yy
    a2b = abs(P2b) - r2b / yy
    denb = 2 * sb - arb(1) / 2
    if a2b > 0 and a1b > 0 and denb > 0:
        num_hi = (a1b.log() - a2b.log()).upper()
        fo = (num_hi / denb.lower()).upper() if num_hi >= 0 else (num_hi / denb.upper()).upper()
        if b is None or fo < b:
            b = fo
    if contrib is None:
        contrib = [float(h[0]) * 30.0] + [float(t) for t in h[1:]]
    return b, contrib


def _naive_check2_fo(self, c, h, y, budget=2000):
    yy = arb(y)
    rhs = (yy / arb.pi()).log() - 1 / yy
    stack = [(list(c), list(h))]
    used = 0
    while stack:
        cc, hh = stack.pop()
        used += 1
        if used > budget:
            return False, used
        b, contrib = self.bound_contrib_fo(cc, hh, y)
        if b is not None and b < rhs:
            continue
        j = max(range(len(contrib)), key=lambda i: contrib[i])
        h2 = list(hh)
        h2[j] = hh[j] / 2
        c1, c2 = list(cc), list(cc)
        c1[j] -= h2[j]
        c2[j] += h2[j]
        stack.append((c2, h2))
        stack.append((c1, list(h2)))
    return True, used


ArbQNaive.bound_contrib_fo = _naive_bound_contrib_fo
ArbQNaive.check2_fo = _naive_check2_fo
