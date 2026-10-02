# Chapter 13, Tables 13.6 and 13.7: certified enclosures of the total variation distances TV(q, p) = (1/2) int |q - p|
# of the approximations q of the Polya density p = Phi/Xi(0), in the book's units,
#   Phi(x) = 2 sum_n (2 pi^2 n^4 e^{9x/2} - 3 pi n^2 e^{5x/2}) exp(-pi n^2 e^{2x}).
# Every density is even, so TV = int_0^infty |q - p|. On [0, A], A = 3, the integral is a sum over the intervals
# [c - r, c + r], r = 2^-12: with d = q - p, alpha = d(c), beta = d'(c), gamma = d''(c)/2 in ball arithmetic
# (python-flint arb),
#   | int |d| - int_{-r}^{r} |alpha + beta h| dh |  <=  |gamma| 2r^3/3 + M3 r^4/12,   M3 >= sup_[0,A] |d'''|,
# and int |alpha + beta h| dh is exact: 2r|alpha| if |alpha| >= |beta| r, else (alpha^2 + beta^2 r^2)/|beta|.
# M3 comes from interval evaluation of the third derivative on a cover of [0, A] (Taylor series in arb), except for the
# laws f_K, where |f_K'''| <= (1/pi) int_0^infty theta^3 |phi_K(theta)| d theta, bounded by interval evaluation of
# |phi_K| on a cover of [0, Theta] and an explicit tail. Beyond A every density is bounded by an explicit function whose
# logarithm decreases at a known rate, and the tail integrals are bounded by value over rate. The normalising constants
# int w p of the tilts and of Polya's kernels come from the same quadrature, with its error bounds.
# The series for Phi is summed to n = 8; the rest is bounded, with its first three derivatives, for x >= 0 by
# 4 (2 pi^2 9^4 + 3 pi 81)(2 pi 81 + 6.5)^3 e^{-81 pi}.
# The laws f_K = g_v(x) sum_j e_j v^-j He_2j(x/sqrt v), v = 2 beta_K, use the certified inputs of fm_qcert.py
# (fm_zeros_hp.txt: Var X and the zeros from Arb), so the distances are those of the true laws.
# Usage: python3 fm_tvcert.py      reads fm_zeros_hp.txt; writes fm_tvcert.json
import json, math, os, sys, time
from multiprocessing import Pool
from flint import arb, acb, arb_series, ctx

A = 3
R = arb(1)/4096                 # half width of the quadrature intervals (exact)
NI = 6144                       # A/(2R)
KS = [0, 1, 2, 3, 5, 10, 16, 17, 20, 50, 100, 150]
CPUS = int(os.environ.get('SLURM_CPUS_PER_TASK', '8'))


def F(x):
    return float(x.mid().str(20, radius=False))


def lo(x):
    """A float below every point of the ball x."""
    v = float(x.lower().str(20, radius=False))
    return v - abs(v)*1e-12 - 1e-300


def hi(x):
    v = float(x.upper().str(20, radius=False))
    return v + abs(v)*1e-12 + 1e-300


def ub(x):
    """Exact arb upper bound of |x|."""
    return arb(abs(x).upper())


def inputs(W):
    ctx.prec = W
    L = [l.strip() for l in open('fm_zeros_hp.txt') if l.strip() and not l.startswith('#')]
    return arb(L[0]), [arb(s) for s in L[1:200]]


def xi_half():
    s = arb(1)/2
    return s*(s - 1)/2*arb.pi()**(-s/2)*(s/2).gamma()*acb(s).zeta().real


def ser(x0):
    return arb_series([x0, 1], prec=4)


def coeffs4(s):
    c = s.coeffs()
    return c + [arb(0)]*(4 - len(c))


def Phi_ser(x0):
    """Taylor coefficients [Phi, Phi', Phi''/2, Phi'''/6] at x0 >= 0 (a point or a ball)."""
    X = ser(x0); PI = arb.pi()
    e9, e5, e2 = (X*arb('4.5')).exp(), (X*arb('2.5')).exp(), (X*2).exp()
    s = arb_series([0], prec=4)
    for n in range(1, 9):
        s += ((e9*(2*PI**2*n**4) - e5*(3*PI*n**2))*(e2*(-PI*n**2)).exp())*2
    return [ci + arb(0, TB) for ci in coeffs4(s)[:4]]


def cosh_ser(X, a):
    return ((X*a).exp() + (X*(-a)).exp())/2


def kernel_ser(name, x0, xpow=0):
    """Taylor coefficients of x^xpow times the unnormalised density, at x0 >= 0."""
    X = ser(x0); PI = arb.pi()
    if name == 'polya_star':
        s = cosh_ser(X, arb('4.5'))*(8*PI**2)*(cosh_ser(X, arb(2))*(-2*PI)).exp()
    elif name == 'psi':
        s = (cosh_ser(X, arb('4.5'))*(8*PI**2) - cosh_ser(X, arb('2.5'))*(12*PI))*(cosh_ser(X, arb(2))*(-2*PI)).exp()
    elif name == 'gauss':
        s = (X*X*(-1/(2*VAR))).exp()/(2*PI*VAR).sqrt()
    else:
        P = arb_series(Phi_ser(x0), prec=4)
        if name == 'cosh_half':
            w = cosh_ser(X, arb(1)/2)
        elif name == 'cosh_sqrt01':
            w = cosh_ser(X, arb('0.1').sqrt())
        elif name == 'gauss_lambda_half':
            w = (X*X/8).exp()
        elif name == 'dbn_0.2':
            w = (X*X/20).exp()
        elif name == 'uniform_sum_quarter':      # sinh(b x)/(b x), b = sqrt(3)/2, by its Taylor series to x^118
            b2 = arb(3)/4
            w = arb_series([1], prec=4); term = arb_series([1], prec=4)
            for k in range(1, 60):
                term = term*(X*X)*b2/((2*k)*(2*k + 1))
                w += term
            # the rest, sum_{k>=60} (bx)^2k/(2k+1)!, and its first three derivatives are at most
            # 2 (120)^3 (bA)^120/121! on [0, A] (A >= 1, ratio of terms below 1/2)
            rem = 2*arb(120)**3*(b2*A*A)**60/arb(122).gamma()
            w = arb_series([wc + arb(0, rem.upper()) for wc in coeffs4(w)], prec=4)
        else:
            raise ValueError(name)
        s = w*P
    if xpow:
        s = s*X**xpow
    return coeffs4(s)[:4]


def M3_cover(fn):
    """Exact arb upper bound of sup_[0,A] |f'''| for f with Taylor coefficients fn(x0), by interval evaluation."""
    m = arb(0)
    n = 768
    for i in range(n):
        x0 = arb(arb(2*i + 1)*A/(2*n), arb(A)/(2*n))
        c = fn(x0)
        u = ub(6*c[3])
        if u > m:
            m = u
    return m


def absint(al, be, r):
    """Ball containing int_{-r}^{r} |al + be h| dh."""
    a, b = abs(al), abs(be)
    if a >= b*r:
        return 2*r*a
    if a < b*r:
        return (al*al + be*be*r*r)/b
    l, u = 2*r*a, 2*r*a + b*r*r
    return arb((l.lower() + u.upper())/2, (u.upper() - l.lower())/2)


def quad(coef, M3):
    """Ball containing int_0^A g, for g with Taylor coefficients coef(i) at the centre of interval i."""
    s = arb(0)
    e = M3*R**4/12
    for i in range(NI):
        c = coef(i)
        s += 2*R*c[0] + 2*R**3/3*c[2] + arb(0, e.upper())
    return s


def tv_sum(coef_d, M3):
    """Ball containing int_0^A |d|."""
    s = arb(0)
    e3 = M3*R**4/12
    for i in range(NI):
        c = coef_d(i)
        e = ub(c[2])*2*R**3/3 + e3
        s += absint(c[0], c[1], R) + arb(0, e.upper())
    return s


def centre(i):
    return arb(2*i + 1)*R


# ------------------------------------------------------------------ the laws f_K
class FK:
    def __init__(self, K):
        self.K = K
        if K == 0:
            self.W = 128
        else:
            g = GF[:K]; S = sum(x**-2 for x in g); v = VARF - 2*S
            self.W = int(1.12*3.33*(30 + K*max(0.0, math.log10(32*S/v)))) + 256
        ctx.prec = self.W
        V, G = inputs(self.W)
        w = [1/(gg*gg) for gg in G[:K]]
        e = [arb(1)] + [arb(0)]*K
        for wi in w:
            for j in range(K, 0, -1):
                e[j] += wi*e[j - 1]
        self.S = sum(w, arb(0)); self.v = V - 2*self.S; self.beta = self.v/2
        self.c = [e[j]/self.v**j for j in range(K + 1)]
        self.sv = self.v.sqrt()
        self.G = G[:K]

    def coef(self, x0):
        """[f, f', f''/2] at the point x0, from f^(m)(x) = g_v(x) (-1)^m v^(-m/2) sum_j c_j He_{2j+m}(x/sqrt v)."""
        ctx.prec = self.W
        y = x0/self.sv
        K = self.K
        H = [arb(1), y]
        for n in range(1, 2*K + 2):
            H.append(y*H[n] - n*H[n - 1])
        g = (-y*y/2).exp()/(2*arb.pi()*self.v).sqrt()
        s = [sum((self.c[j]*H[2*j + m] for j in range(K + 1)), arb(0)) for m in range(3)]
        return [g*s[0], -g*s[1]/self.sv, g*s[2]/self.v/2]

    def M3(self):
        """Bound of (1/pi) int_0^infty theta^3 |phi_K(theta)|, phi_K = exp(-beta th^2) prod (1 - th^2/gamma_k^2)."""
        ctx.prec = 64
        beta = arb(self.beta.lower())
        K = self.K
        G = [arb(gg) for gg in self.G]
        m = 2*K + 3
        Th = max(1.5*(F(G[-1]) if K else 0.0), math.sqrt(2*m/F(beta)) + 1, 10.0)
        n = int(Th*64) + 1                      # cover [0, n/64] by intervals of width 1/64
        tot = arb(0)
        for i in range(n):
            a = arb(i)/64; b = arb(i + 1)/64
            th = arb(arb(2*i + 1)/128, arb(1)/128)
            p = (-beta*a*a).exp()
            for gk in G:
                p *= abs(1 - th*th/(gk*gk))
            tot += arb(p.upper())*b**3/64
        # tail: |phi_K| <= exp(-beta th^2) th^2K / prod gamma_k^2 for th >= T >= sqrt 2 gamma_K, and
        # int_T^inf th^m e^{-beta th^2} <= T^m e^{-beta T^2} / (2 beta T - m/T) when 2 beta T^2 > m.
        T = arb(n)/64
        pg = arb(1)
        for gk in G:
            pg *= gk*gk
        if K:
            assert T*T >= 2*G[-1]**2
        rate = 2*beta*T - m/T
        assert rate > 0
        tail = T**m*(-beta*T*T).exp()/rate/pg
        return arb(((tot + tail)/arb.pi()).upper()), hi(tail)

    def tail(self):
        """int_A^inf f_K <= g_v(A) Pi(A/sqrt v)/lam, since 0 < Q(y) <= Pi(y) = prod (1 + y^2/(v gamma_k^2)) for
        y >= sqrt(8K+2), and the log of g_v(x) Pi(x/sqrt v) decreases at rate >= A/v - 2K/A for x >= A."""
        ctx.prec = 128
        y = arb(A)/self.sv
        if self.K == 0:
            return (-y*y/2).exp()*self.v/A/(2*arb.pi()*self.v).sqrt()
        assert y*y >= 8*self.K + 2
        Pi = arb(1)
        for gk in self.G:
            Pi *= 1 + y*y/(self.v*gk*gk)
        lam = A/self.v - arb(2*self.K)/A
        assert lam > 0
        return (-y*y/2).exp()/(2*arb.pi()*self.v).sqrt()*Pi/lam


def tail_Phi(extra=0):
    """Bound of int_A^inf e^{x + x^2/8} x^extra Phi(x) dx (covers every weight used, and p itself):
    Phi(x) <= 4 (2 pi^2 + 3 pi) e^{9x/2} exp(-pi e^{2x}) for x >= 0, x^extra <= e^{extra x}, and the log of the bound
    decreases at rate >= 2 pi e^{2A} - 9/2 - 1 - A/4 - extra for x >= A."""
    PI = arb.pi()
    b = 4*(2*PI**2 + 3*PI)*((arb('5.5') + extra)*A + arb(A)**2/8).exp()*(-PI*(arb(2)*A).exp()).exp()
    rate = 2*PI*(arb(2)*A).exp() - arb('5.5') - arb(A)/4 - extra
    assert rate > 0
    return b/rate


def tail_kernel(extra=0):
    """Bound of int_A^inf x^extra |Polya kernel|: (8 pi^2 cosh 9x/2 + 12 pi cosh 5x/2) exp(-2 pi cosh 2x)
    <= 32 pi^2 e^{9x/2} exp(-pi e^{2x})."""
    PI = arb.pi()
    b = 32*PI**2*((arb('4.5') + extra)*A).exp()*(-PI*(arb(2)*A).exp()).exp()
    rate = 2*PI*(arb(2)*A).exp() - arb('4.5') - extra
    return b/rate


def run_fk(K):
    t0 = time.time()
    f = FK(K)
    M3f, th_tail = f.M3()
    M3 = M3f + M3P
    def cd(i):
        a = f.coef(centre(i))
        ctx.prec = 128
        pp = Pc[i]
        return [a[0] - pp[0], a[1] - pp[1], a[2] - pp[2]]
    s = tv_sum(cd, M3)
    ctx.prec = 128
    tl = f.tail() + TP
    tv = s + arb(0, tl.upper())
    out = dict(name='f_%d' % K, K=K, W=f.W, TV=[lo(tv), hi(tv)], M3=hi(M3), M3_fK=hi(M3f), fourier_tail=th_tail,
               tails=hi(tl), two_beta=f.v.str(12), secs=round(time.time() - t0, 1))
    print(json.dumps(out), flush=True)
    return out


def run_kernel(name):
    t0 = time.time()
    ctx.prec = 128
    fn = lambda x0: kernel_ser(name, x0)
    M3k = M3_cover(fn)
    cs = [fn(centre(i)) for i in range(NI)]
    if name in ('polya_star', 'psi'):
        tk = tail_kernel()
    elif name == 'gauss':
        tk = (-arb(A)**2/(2*VAR)).exp()*VAR/A/(2*arb.pi()*VAR).sqrt()
    else:
        tk = tail_Phi()
    m = quad(lambda i: cs[i], M3k)
    half = m + arb(0, tk.upper())                 # int_0^inf of the unnormalised density
    C2 = arb(1) if name == 'gauss' else 2*half    # its integral over the line
    M3 = M3k/C2 + M3P
    def cd(i):
        a = cs[i]; pp = Pc[i]
        return [a[0]/C2 - pp[0], a[1]/C2 - pp[1], a[2]/C2 - pp[2]]
    s = tv_sum(cd, M3)
    tl = tk/C2 + TP
    tv = s + arb(0, tl.upper())
    out = dict(name=name, mass=[lo(C2), hi(C2)], TV=[lo(tv), hi(tv)], M3=hi(M3), tails=hi(tl), secs=round(time.time() - t0, 1))
    if name in ('polya_star', 'psi'):
        fn2 = lambda x0: kernel_ser(name, x0, 2)
        v2 = quad(lambda i: fn2(centre(i)), M3_cover(fn2)) + arb(0, tail_kernel(2).upper())
        var = 2*v2/C2
        out['var'] = [lo(var), hi(var)]
        q0 = kernel_ser(name, arb(0))[0]/C2
        out['at0'] = [lo(q0), hi(q0)]
    print(json.dumps(out), flush=True)
    return out


def setup():
    global VAR, VARF, GF, XI0, TB, Pc, M3P, TP
    ctx.prec = 128
    VAR, G = inputs(128)
    VARF = F(VAR); GF = [F(g) for g in G]
    XI0 = xi_half()
    PI = arb.pi()
    TB = (4*(2*PI**2*9**4 + 3*PI*81)*(2*PI*81 + arb('6.5'))**3*(-81*PI).exp()).upper()
    Pc = []
    for i in range(NI):
        c = Phi_ser(centre(i))
        Pc.append([c[0]/XI0, c[1]/XI0, c[2]/XI0])
    M3P = M3_cover(lambda x0: [ci/XI0 for ci in Phi_ser(x0)])
    TP = tail_Phi()/XI0


if __name__ == '__main__':
    t0 = time.time()
    setup()
    print('Xi(0) = %s, Var X = %s, sup|p\'\'\'| <= %.4g, tail of p beyond %d <= %.3g, series tail %.3g'
          % (XI0.str(15), VAR.str(15), hi(M3P), A, hi(TP), F(arb(TB))), flush=True)
    kernels = ['gauss', 'polya_star', 'psi', 'cosh_half', 'uniform_sum_quarter', 'gauss_lambda_half', 'dbn_0.2',
               'cosh_sqrt01']
    with Pool(CPUS) as pool:
        rf = pool.map_async(run_fk, sorted(KS, reverse=True), chunksize=1)
        rk = pool.map_async(run_kernel, kernels, chunksize=1)
        res = rk.get() + rf.get()
    ctx.prec = 128
    p0 = Phi_ser(arb(0))[0]/XI0
    # checks: int cosh(x/2) Phi = xi(1) = 1/2, int cosh(sqrt(0.1) x) Phi = xi(1/2 + sqrt 0.1)
    s = arb(1)/2 + arb('0.1').sqrt()
    xi_s = s*(s - 1)/2*arb.pi()**(-s/2)*(s/2).gamma()*acb(s).zeta().real
    out = dict(A=A, R=F(R), NI=NI, Xi0=XI0.str(20), VarX=VAR.str(20), p0=p0.str(10), M3_p=hi(M3P), tail_p=hi(TP),
               xi_check=dict(cosh_half=0.5, cosh_sqrt01=xi_s.str(15)), results=res, secs=round(time.time() - t0, 1))
    json.dump(out, open('fm_tvcert.json', 'w'), indent=1)
    for r in res:
        print('%-22s TV in [%.6e, %.6e]%s' % (r['name'], r['TV'][0], r['TV'][1],
                                             '  mass in [%.12f, %.12f]' % tuple(r['mass']) if 'mass' in r else ''))
    print('xi(1/2 + sqrt 0.1) =', xi_s.str(15), '(the mass for cosh_sqrt01); 1/2 for cosh_half')
    print('done in %.0f s' % (time.time() - t0))
