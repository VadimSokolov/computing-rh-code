# Chapter 13, Section ch:fmexp, Proposition fm:prop:hcm: the HCM threshold of e^{beta X}, a computer assisted proof in
# ball arithmetic (python-flint arb and acb).
# The computation runs in the variable u = t/2 of Rodgers and Tao, with the kernel Phi_N(u) = Phi(2u)/2
#   Phi_N(u) = sum_n (2 pi^2 n^4 e^{9u} - 3 pi n^2 e^{5u}) exp(-pi n^2 e^{4u}),  analytic in |Im u| < pi/8,
# whose law is that of X/2. So e^{kappa X/2} = e^{beta X} with beta = kappa/2, and the printed book thresholds are half
# of the thresholds in kappa.
# By the elasticity criterion (Theorem bs:thm:crit, as in Proposition bs:prop:Wcrit), e^{kappa X/2} has an HCM density
# iff Phi_N has no zeros in |Im u| < pi/kappa and h(x, v) = Im(-Phi_N'/Phi_N)(x + iv) >= 0 for 0 < v < pi/kappa.
#  (Z) for x >= 0, 0 <= v <= Vz: |Phi_N/phi_1 - 1| <= S(Vz) = sum_{n>=2} n^2 (2 pi n^2 + 3)/(2 pi - 3) exp(-pi (n^2-1) cos 4Vz),
#      and S(Vz) < 1 gives Phi_N != 0; by Phi_N(-x+iv) = conj Phi_N(x+iv) the same holds for x <= 0.
#  (T) top line: h(x, V) > 0 for 0 <= x <= X = 1, by adaptive bisection with rigorous enclosures (the tail n > 30 of the
#      series is bounded explicitly); h is even in x.
#  (S) side and outer region x >= X: h >= 4 pi e^{4x} sin(4v) - |Im q'|, with q = log(1 + eps), eps = (Phi_N - phi_1)/phi_1,
#      |Im q'(x+iv)| <= v max|q''| <= v 4 S_eps/r^2 (Cauchy on discs of radius r), and sin 4v >= 8v/pi.
#  Then h = 0 on v = 0, h is harmonic on the strip (no zeros), h >= 0 on the boundary of [-X, X] x [0, V], so h >= 0 inside
#  (minimum principle), and h > 0 for |x| >= X by (S). Hence e^{kappa X/2} is HCM for every kappa >= pi/V.
#  (A) necessity: F(v) = Phi_N(iv) is real and F'(v) = -Im Phi_N'(iv); at a point v_b with F(v_b) > 0 > F'(v_b), h(0, v_b) < 0,
#      so the condition fails on the imaginary axis and e^{kappa X/2} is not HCM for kappa < pi/v_b. The bisection brackets
#      a local maximum v_m of F, with certified signs F' > 0 at the left end and F' < 0 at the right end v_b; that v_m is the
#      first critical point is not certified, and the necessity bound uses only v_b.
# Usage: python3 fm_hcm.py   Writes fm_hcm.json.
import sys, time, json
from flint import arb, acb, ctx

ctx.prec = 256
PI = arb.pi()
N = 30
t0 = time.time()
out = {}

def fl(a):
    try:
        return float(a)
    except Exception:
        return float(a.mid().str(20, radius=False))

def tail_bounds(X, Vmax):
    """Rigorous bounds for sum_{n>N} |phi_n(u)| and |phi_n'(u)| on -0.001 <= Re u <= X, |Im u| <= Vmax."""
    c = (4*Vmax).cos().lower()*arb(-0.004).exp()          # lower bound for Re e^{4u}
    e9 = (9*X).exp(); e5 = (5*X).exp(); e4 = (4*X).exp()
    def T(n):
        A = 2*PI**2*n**4*e9 + 3*PI*n**2*e5
        Ad = 18*PI**2*n**4*e9 + 15*PI*n**2*e5 + A*4*PI*n**2*e4
        E = (-PI*c*n**2).exp()
        return A*E, Ad*E
    n = N + 1
    ratio = arb(n + 1)**6/arb(n)**6*(-PI*c*(2*n + 1)).exp()
    assert ratio < arb(0.5), 'tail ratio too large'
    tP, tD = T(n)
    return (2*tP).upper(), (2*tD).upper()

def PhiD(u, tP, tD):
    e4 = (4*u).exp(); e5 = (5*u).exp(); e9 = (9*u).exp()
    P = acb(0); D = acb(0)
    for n in range(1, N + 1):
        n2 = n*n
        A = 2*PI**2*(n2*n2)*e9 - 3*PI*n2*e5
        Ad = 18*PI**2*(n2*n2)*e9 - 15*PI*n2*e5
        E = (-PI*n2*e4).exp()
        P += A*E
        D += (Ad - A*4*PI*n2*e4)*E
    U = arb(0, 1)
    P += acb(tP*U, tP*U)
    D += acb(tD*U, tD*U)
    return P, D

def zero_free_S(Vz):
    """S(Vz) = sum_{n>=2} n^2 (2 pi n^2+3)/(2 pi-3) exp(-pi(n^2-1) cos 4Vz), with a tail bound."""
    c = (4*Vz).cos().lower()
    s = arb(0)
    for n in range(2, 41):
        s += arb(n)**2*(2*PI*n**2 + 3)/(2*PI - 3)*(-PI*(n**2 - 1)*c).exp()
    n = 41
    s += 2*arb(n)**2*(2*PI*n**2 + 3)/(2*PI - 3)*(-PI*(n**2 - 1)*c).exp()   # ratio of terms < 1/2 beyond 40
    return s

def side_bound(X, V, r):
    """For x >= X, 0 < v <= V: returns (lower bound for (h/v) from the phi_1 term, bound for max|q''|)."""
    y0 = (4*(X - r)).exp()                   # |e^{4u}| >= y0 on the discs
    c = (4*(V + r)).cos().lower()            # Re e^{4u}/|e^{4u}| >= c on the discs (needs V + r < pi/8)
    Se = arb(0)
    for n in range(2, 41):
        Se += arb(n)**2*(2*PI*n**2*y0 + 3)/(2*PI*y0 - 3)*(-PI*(n**2 - 1)*y0*c).exp()
    Se = 2*Se                                 # crude allowance for n > 40 (terms decay superexponentially)
    assert Se < arb(0.5)
    q2 = 4*Se/r**2                           # |log(1+eps)| <= 2|eps|, Cauchy: |q''| <= 2 max|q| / r^2
    main = 32*(4*X).exp()                    # 4 pi e^{4x} sin 4v >= 4 pi e^{4X} (8 v / pi)
    return main, q2, Se

def F_and_dF(v, tP, tD):
    P, D = PhiD(acb(arb(0), v), tP, tD)
    return P.real, -D.imag, P.imag, D.real          # F, F', and the parts that must contain 0

def certify_line(V, X, tP, tD, init=256, minw=2.0**-34, maxeval=400000):
    stack = [(X*k/init, X*(k + 1)/init) for k in range(init - 1, -1, -1)]
    ok = 0; fails = []; minlow = None; where = None; nev = 0
    while stack and nev < maxeval:
        a, b = stack.pop()
        u = acb(arb((a + b)/2, (b - a)/2), V)
        P, D = PhiD(u, tP, tD); nev += 1
        h = (-D/P).imag
        if h > 0:
            ok += 1
            lo = fl(h.lower())
            if minlow is None or lo < minlow:
                minlow, where = lo, (a, b)
        elif b - a < minw:
            fails.append((a, b, h.str(6)))
        else:
            m = (a + b)/2
            stack += [(m, b), (a, m)]
    return dict(intervals=ok, failures=fails[:5], nfail=len(fails), unfinished=len(stack), evaluations=nev,
                min_lower_bound=minlow, at=where)


ctx.prec = 256
PI = arb.pi()
t0 = time.time()
X = 1.0
tP, tD = tail_bounds(arb('1.01'), arb('0.33'))
out = {}
S = zero_free_S(arb('0.2817'))
print('zero free |Im u| <= 0.2817:', S < 1, S.str(6))
# bracket for v_m with certified signs of F'
lo, hi = arb('0.28'), arb('0.283')
for it in range(48):
    mid = (lo + hi)/2
    F, dF, _, _ = F_and_dF(mid, tP, tD)
    if dF > 0: lo = mid
    elif dF < 0: hi = mid
    else:
        print('undecided at', mid.str(20)); break
Fl, dFl, _, _ = F_and_dF(lo, tP, tD)
Fh, dFh, _, _ = F_and_dF(hi, tP, tD)
print('v_lo =', lo.str(17, radius=False), ' F\'(v_lo) =', dFl.str(5), ' F(v_lo) =', Fl.str(12))
print('v_hi =', hi.str(17, radius=False), ' F\'(v_hi) =', dFh.str(5), ' F(v_hi) =', Fh.str(12),
      ' necessity point certified:', bool(Fh > 0 and dFh < 0))
res = {}
for name, V in [('pi/12', PI/12), ('0.281548', arb('0.281548'))]:
    r = certify_line(V, X, tP, tD, minw=2.0**-40)
    print('top line V = %s: %s' % (name, r), flush=True)
    res[name] = r
for Vs in ['0.2618', '0.281548']:
    main, q2, Se = side_bound(arb(X), arb(Vs), arb('0.05'))
    print('side x >= 1, V = %s: h/v >= %s - %s, positive %s' % (Vs, main.str(6), q2.str(3), bool(main - q2 > 0)))
kl = PI/hi; ku = PI/arb('0.281548')
print('RESULT: not HCM for |kappa| < %s; HCM for |kappa| >= %s (variable u = t/2).'
      % (kl.lower().str(12, radius=False), ku.upper().str(12, radius=False)))
print('        book units e^{beta X}: not HCM for |beta| < %s; HCM for |beta| >= %s.'
      % ((kl/2).lower().str(12, radius=False), (ku/2).upper().str(12, radius=False)))
print('        kappa = 12 (beta = 6) is HCM:', res['pi/12']['nfail'] == 0 and res['pi/12']['unfinished'] == 0)
out = dict(v_lo=lo.str(17, radius=False), v_hi=hi.str(17, radius=False), kappa_not_below=kl.lower().str(12, radius=False),
           kappa_hcm_from=ku.upper().str(12, radius=False), beta_not_below=(kl/2).lower().str(12, radius=False),
           beta_hcm_from=(ku/2).upper().str(12, radius=False), v_m_book=[(2*lo).str(12, radius=False), (2*hi).str(12, radius=False)],
           zero_free_S_02817=S.str(8), top=res)
json.dump(out, open('fm_hcm.json', 'w'), indent=1, default=str)
print('total %.0f s' % (time.time() - t0))
