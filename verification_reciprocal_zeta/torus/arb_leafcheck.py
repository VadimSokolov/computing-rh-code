"""Independent Arb check of branch and bound leaves (python-flint ball arithmetic, no hand made radius formulas).

For a box with centre c and half widths h (x first, then the phases), the mean value form
  F(box) subset F(c) + sum_j [-1,1] sup_box |d_j F| h_j
is evaluated with F(c) computed at the exact centre and d_j F computed by Arb over the balls x = c_0 +/- h_0,
w_p = c_p +/- h_p. F is |P|^2 (kind m, lower bound), Re(N/P) with N = -P' (kind M, upper bound) or |N/P|^2
(kind l, upper bound). If Arb cannot decide a box, it is bisected along its widest scaled coordinate, up to a depth.
"""
from fractions import Fraction as Fr
from flint import arb, acb, fmpq, ctx
from torus_common import exact_coeffs, primes_upto, vp

ctx.prec = 64


class ArbTorus:
    def __init__(self, n):
        self.n = n
        ab = exact_coeffs(n)
        self.beta = [arb(fmpq(b.numerator, b.denominator)) for a, b in ab]
        self.ps = primes_upto(n)
        self.V = [[vp(k, p) for p in self.ps] for k in range(1, n + 1)]
        self.logk = [arb(k).log() for k in range(1, n + 1)]
        self.K = len(self.ps)

    def eval(self, x, w, grad=True):
        """x arb, w list of arb. Returns P, N, dP[j], dN[j] as acb (j = 0 for x, 1.. for phases)."""
        P = acb(0)
        N = acb(0)
        K1 = 1 + self.K
        dP = [acb(0)] * K1
        dN = [acb(0)] * K1
        for k in range(self.n):
            phi = arb(0)
            Vk = self.V[k]
            for j in range(self.K):
                if Vk[j]:
                    phi += Vk[j] * w[j]
            L = self.logk[k]
            t = self.beta[k] * (-2 * x * L).exp() * acb(0, -phi).exp()
            tn = 2 * L * t
            P += t
            N += tn
            if grad:
                dP[0] = dP[0] - tn
                dN[0] = dN[0] - 2 * L * tn
                for j in range(self.K):
                    if Vk[j]:
                        dP[j + 1] = dP[j + 1] + acb(0, -Vk[j]) * t
                        dN[j + 1] = dN[j + 1] + acb(0, -Vk[j]) * tn
        return P, N, dP, dN


def bound_box(AT, kind, c, h):
    """Arb bound over the box: lower bound of |P|^2 (m) or upper bound of F (M, l), as an arb upper/lower."""
    xc = arb(c[0])
    wc = [arb(t) for t in c[1:]]
    Pc, Nc, _, _ = AT.eval(xc, wc, grad=False)
    xb = arb(c[0], h[0])
    wb = [arb(t, r) for t, r in zip(c[1:], h[1:])]
    Pb, Nb, dPb, dNb = AT.eval(xb, wb, grad=True)
    hs = [arb(t) for t in h]
    if kind == "m":
        F = abs(Pc) ** 2
        s = arb(0)
        for j in range(len(hs)):
            g = 2 * (Pb.conjugate() * dPb[j]).real
            s += abs(g).upper() * hs[j]
        return F.lower() - s.upper()
    # ratio forms need P != 0 on the box
    if not (abs(Pb) > 0):
        return None
    if kind == "M":
        F = (Nc / Pc).real
        s = arb(0)
        for j in range(len(hs)):
            Z = (dNb[j] * Pb - Nb * dPb[j]) / (Pb * Pb)
            s += abs(Z.real).upper() * hs[j]
        return F.upper() + s.upper()
    R = Nb / Pb
    F = abs(Nc / Pc) ** 2
    s = arb(0)
    for j in range(len(hs)):
        Z = (dNb[j] * Pb - Nb * dPb[j]) / (Pb * Pb)
        g = 2 * (R.conjugate() * Z).real
        s += abs(g).upper() * hs[j]
    return F.upper() + s.upper()


def check(AT, kind, target, c, h, depth=0, maxdepth=8):
    """True if Arb proves the target on the box (bisecting if needed). Returns (ok, boxes_used)."""
    tgt = arb(target) ** 2 if kind in ("m", "l") else arb(target)
    b = bound_box(AT, kind, c, h)
    if b is not None:
        if kind == "m" and b > tgt:
            return True, 1
        if kind != "m" and b < tgt:
            return True, 1
    if depth >= maxdepth:
        return False, 1
    # bisect along the coordinate with the largest scaled half width
    sc = [h[0] * 30.0] + list(h[1:])
    j = max(range(len(sc)), key=lambda i: sc[i])
    h2 = list(h)
    h2[j] = h[j] / 2
    c1 = list(c)
    c2 = list(c)
    c1[j] -= h2[j]
    c2[j] += h2[j]
    ok1, n1 = check(AT, kind, target, c1, h2, depth + 1, maxdepth)
    if not ok1:
        return False, n1
    ok2, n2 = check(AT, kind, target, c2, h2, depth + 1, maxdepth)
    return ok2, n1 + n2


def bound_box_contrib(AT, kind, c, h):
    """As bound_box, but also returns the per coordinate contributions sup|d_j F| h_j (floats, upper bounds)."""
    xc = arb(c[0])
    wc = [arb(t) for t in c[1:]]
    Pc, Nc, _, _ = AT.eval(xc, wc, grad=False)
    xb = arb(c[0], h[0])
    wb = [arb(t, r) for t, r in zip(c[1:], h[1:])]
    Pb, Nb, dPb, dNb = AT.eval(xb, wb, grad=True)
    hs = [arb(t) for t in h]
    contrib = []
    if kind == "m":
        F = abs(Pc) ** 2
        for j in range(len(hs)):
            g = 2 * (Pb.conjugate() * dPb[j]).real
            contrib.append(abs(g).upper() * hs[j])
        s = sum(contrib, arb(0))
        return F.lower() - s.upper(), [float(t.upper()) for t in contrib]
    if not (abs(Pb) > 0):
        return None, None
    if kind == "M":
        F = (Nc / Pc).real
        for j in range(len(hs)):
            Z = (dNb[j] * Pb - Nb * dPb[j]) / (Pb * Pb)
            contrib.append(abs(Z.real).upper() * hs[j])
        s = sum(contrib, arb(0))
        return F.upper() + s.upper(), [float(t.upper()) for t in contrib]
    R = Nb / Pb
    F = abs(Nc / Pc) ** 2
    for j in range(len(hs)):
        Z = (dNb[j] * Pb - Nb * dPb[j]) / (Pb * Pb)
        g = 2 * (R.conjugate() * Z).real
        contrib.append(abs(g).upper() * hs[j])
    s = sum(contrib, arb(0))
    return F.upper() + s.upper(), [float(t.upper()) for t in contrib]


def check2(AT, kind, target, c, h, budget=4000):
    """Arb proof of the target on the box, bisecting along the coordinate with the largest Arb contribution.
    Returns (ok, boxes_used). Depth first with a box budget."""
    tgt = arb(target) ** 2 if kind in ("m", "l") else arb(target)
    stack = [(list(c), list(h))]
    used = 0
    while stack:
        cc, hh = stack.pop()
        used += 1
        if used > budget:
            return False, used
        b, contrib = bound_box_contrib(AT, kind, cc, hh)
        if b is not None:
            if kind == "m" and b > tgt:
                continue
            if kind != "m" and b < tgt:
                continue
            j = max(range(len(contrib)), key=lambda i: contrib[i])
        else:
            sc = [hh[0] * 30.0] + list(hh[1:])
            j = max(range(len(sc)), key=lambda i: sc[i])
        h2 = list(hh)
        h2[j] = hh[j] / 2
        c1 = list(cc)
        c2 = list(cc)
        c1[j] -= h2[j]
        c2[j] += h2[j]
        stack.append((c2, h2))
        stack.append((c1, list(h2)))
    return True, used


def bound_box_contrib_mlow(AT, kind, c, h, m_low):
    """Naive Arb bound for kinds M, l that also uses a certified lower bound m_low of |P| on the band: where the Arb
    ball of P over the box contains 0, |1/P| <= 1/m_low is used instead of dividing by the ball (cruder, valid)."""
    xc = arb(c[0])
    wc = [arb(t) for t in c[1:]]
    Pc, Nc, _, _ = AT.eval(xc, wc, grad=False)
    xb = arb(c[0], h[0])
    wb = [arb(t, r) for t, r in zip(c[1:], h[1:])]
    Pb, Nb, dPb, dNb = AT.eval(xb, wb, grad=True)
    hs = [arb(t) for t in h]
    ml = arb(m_low)
    lo = abs(Pb).lower()
    A = lo if lo > ml else ml          # |P| >= A on the box
    contrib = []
    for j in range(len(hs)):
        num = dNb[j] * Pb - Nb * dPb[j]
        if kind == "M":
            g = abs(num).upper() / (A * A)   # |d_j Re(N/P)| <= |num| / |P|^2
        else:
            g = 2 * (abs(Nb).upper() / A) * abs(num).upper() / (A * A)
        contrib.append((g * hs[j]).upper())
    s = sum(contrib, arb(0))
    if kind == "M":
        F = (Nc / Pc).real
        return (F + s).upper(), [float(t) for t in contrib]
    F = abs(Nc / Pc) ** 2
    return (F + s).upper(), [float(t) for t in contrib]


def check2_mlow(AT, kind, target, c, h, m_low, budget=4000):
    """As check2, with bound_box_contrib_mlow."""
    tgt = arb(target) ** 2 if kind == "l" else arb(target)
    stack = [(list(c), list(h))]
    used = 0
    while stack:
        cc, hh = stack.pop()
        used += 1
        if used > budget:
            return False, used
        b, contrib = bound_box_contrib_mlow(AT, kind, cc, hh, m_low)
        if b < tgt:
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
