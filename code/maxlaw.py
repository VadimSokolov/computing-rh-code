# Maxima of exponentials and the law of Ramanujan's Delta (Section 8.5.4). M* is the maximum of independent exponential
# variables with rates 2 pi n, each rate taken 24 times, so P(M* <= y) = F(y) = prod (1 - e^{-2 pi n y})^24 = e^{2 pi y} Delta(iy);
# M_1 is the maximum with each rate taken once, P(M_1 <= y) = F_1(y) = prod (1 - e^{-2 pi n y}) = e^{pi y/12} eta(iy).
# Checks, in the normalisation xi(s, Delta) = (2 pi)^{-s-11/2} Gamma(s + 11/2) L(s, Delta) of Section 8.5.2:
#   1. xi(s, Delta) = (2 pi)^{-s-11/2} E Gamma(s + 11/2, 2 pi M*) and the symmetric form xi(s) = A(s) + A(1 - s),
#      A(s) = (2 pi)^{-s-11/2} E Gamma(s + 11/2, 2 pi max(M*, 1)), against Hecke's formula; the first zero from the symmetric form;
#   2. L(s, chi_12) = E Q(s/2, pi M_1/12) at s = 1;
#   3. the Laplace transform of M_1: closed form, product and quadrature; the factorisation M_1 + Z = S into two GGCs;
#      the characteristic function at u = 5 in three forms; the Levy density by Poisson summation; the Bondesson function N_1;
#   4. the Laplace transform psi* of M*: its poles, its real zero between -12 pi and -10 pi, the sign of psi* along the real
#      axis, and a count of its zeros in a rectangle by the argument principle;
#   5. the first zero of Polya's factor K_{it}(2 pi), and the two halves of Hecke's cut at y = 1 off the critical line.
# mpmath only. Writes maxlaw.json.
import json, time
import mpmath as mp

mp.mp.dps = 40
t0 = time.time()
out = {}
def put(key, val):
    out[key] = val
    print(key, val, round(time.time()-t0, 1), flush=True)

# tau(n) from (n - 1) tau(n) = -24 sum_{m < n} sigma(m) tau(n - m), which is q d/dq log Delta = E_2
NT = 400
sig = [0]*(NT+1)
for d in range(1, NT+1):
    for m in range(d, NT+1, d): sig[m] += d
tau = [0]*(NT+1); tau[1] = 1
for n in range(2, NT+1):
    tau[n] = -24*sum(sig[m]*tau[n-m] for m in range(1, n))//(n-1)
assert tau[2] == -24 and tau[3] == 252 and tau[11] == 534612
put('tau_2_to_8', tau[2:9])

H = mp.mpf(11)/2
TINY = mp.mpf(10)**(-mp.mp.dps-8)

# 1. xi(s, Delta)
def A_hecke(s, N=30):
    a = s + H
    return mp.fsum(tau[n]*(2*mp.pi*n)**(-a)*mp.gammainc(a, 2*mp.pi*n) for n in range(1, N+1))
def xi_hecke(s):
    return A_hecke(s) + A_hecke(1-s)
def F(y):
    y = mp.mpf(y)
    if y < 1:
        return y**(-12)*mp.exp(2*mp.pi*y - 2*mp.pi/y)*F(1/y)
    p = mp.mpf(1); n = 1
    while True:
        q = mp.exp(-2*mp.pi*n*y)
        if q < TINY: break
        p *= (1-q)**24; n += 1
    return p
def dlogF(y):
    # (log F)'(y) = 48 pi sum sigma(m) e^{-2 pi m y} for y >= 1; for y < 1 differentiate F(y) = y^-12 e^{2 pi y - 2 pi/y} F(1/y)
    if y < 1:
        return -12/y + 2*mp.pi + 2*mp.pi/y**2 - dlogF(1/y)/y**2
    return 48*mp.pi*mp.fsum(sig[m]*mp.exp(-2*mp.pi*m*y) for m in range(1, 80))
def dens(y):
    return F(y)*dlogF(y)
BRK = [0, 0.15, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24, 32]
def xi_max(s):
    a = s + H
    return (2*mp.pi)**(-a)*mp.quad(lambda y: mp.gammainc(a, 2*mp.pi*y)*dens(y), BRK)
def A_max(s):
    return mp.quad(lambda y: mp.exp(-2*mp.pi*y)*F(y)*y**(s + H - 1), [1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24, 32])
put('density_integral', mp.nstr(mp.quad(dens, BRK), 25))
chk = {}
for s in [mp.mpf(1)/2, mp.mpc(0.5, 5), mp.mpf(3)]:
    chk[mp.nstr(s, 4)] = {'Hecke': mp.nstr(xi_hecke(s), 25), 'E Gamma(s+11/2, 2 pi M*)': mp.nstr(xi_max(s), 25),
                          'A(s)+A(1-s), A from the law of max(M*,1)': mp.nstr(A_max(s) + A_max(1-s), 25)}
put('xi_Delta', chk)
# L(s, Delta) = E Q(s + 11/2, 2 pi M*) at s = 17/2, against the absolutely convergent Dirichlet series
a = mp.mpf(14)
put('L(17/2, Delta)', {'E Q(14, 2 pi M*)': mp.nstr(mp.quad(lambda y: mp.gammainc(a, 2*mp.pi*y, regularized=True)*dens(y), BRK), 25),
                       'Hecke': mp.nstr(xi_hecke(mp.mpf(17)/2)*(2*mp.pi)**a/mp.gamma(a), 25),
                       'sum tau(n) n^-14, n <= 400 (tail below 1e-19)': mp.nstr(mp.fsum(mp.mpf(tau[n])/mp.mpf(n)**14 for n in range(1, NT+1)), 18)})
z_sym = mp.findroot(lambda t: mp.re(A_max(mp.mpc(0.5, t))), mp.mpf('9.22'))
z_hec = mp.findroot(lambda t: mp.re(A_hecke(mp.mpc(0.5, t))), mp.mpf('9.22'))
put('first_zero', {'symmetric form from the law of max(M*,1)': mp.nstr(z_sym, 20), 'Hecke': mp.nstr(z_hec, 20)})

# 2. eta and L(s, chi_12)
def F1(y):
    y = mp.mpf(y)
    if y < 1:
        return y**(-mp.mpf(1)/2)*mp.exp(mp.pi*y/12 - mp.pi/(12*y))*F1(1/y)
    p = mp.mpf(1); n = 1
    while True:
        q = mp.exp(-2*mp.pi*n*y)
        if q < TINY: break
        p *= 1-q; n += 1
    return p
def F1m1(y):
    # F_1(y) - 1 without cancellation for large y
    y = mp.mpf(y)
    if y < 1:
        return F1(y) - 1
    s = mp.mpf(0); n = 1
    while True:
        q = mp.exp(-2*mp.pi*n*y)
        if q < TINY: break
        s += mp.log1p(-q); n += 1
    return mp.expm1(s)
def dlogF1(y):
    if y < 1:
        return -1/(2*y) + mp.pi/12 + mp.pi/(12*y**2) - dlogF1(1/y)/y**2
    return 2*mp.pi*mp.fsum(sig[m]*mp.exp(-2*mp.pi*m*y) for m in range(1, 80))
BRK1 = [0, 0.15, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 100]
L1 = mp.quad(lambda y: mp.erfc(mp.sqrt(mp.pi*y/12))*F1(y)*dlogF1(y), BRK1)
put('L(1, chi_12)', {'E erfc(sqrt(pi M_1/12))': mp.nstr(L1, 25), 'log(2 + sqrt 3)/sqrt 3': mp.nstr(mp.log(2+mp.sqrt(3))/mp.sqrt(3), 25)})

# 3. the eta law M_1
b = lambda n: mp.pi*(n*n-1)/12
def wof(th):
    return mp.sqrt(1 - 12*mp.mpc(th)/mp.pi)
def lt1_closed(th):
    w = wof(th)
    return mp.pi*(1-w**2)*mp.sin(mp.pi*w/3)/(2*mp.sqrt(3)*w*mp.cos(mp.pi*w/2))
def lt1_prod(th):
    return mp.nprod(lambda k: (1+th/b(6*k))/((1+th/b(6*k-1))*(1+th/b(6*k+1))), [1, mp.inf])
def lt1_quad(th):
    pts = [mp.mpf(x) for x in BRK1[:6]] + [mp.mpf(x)/2 for x in range(3, 201)]
    return th*mp.quad(lambda y: mp.exp(-th*y)*F1(y), pts)
def ltS_closed(th):
    w = wof(th)
    return mp.pi/4*(1-w**2)/mp.cos(mp.pi*w/2)
def ltS_prod(th):
    return mp.nprod(lambda j: 1/(1+th/b(2*j+1)), [1, mp.inf])
def ltZ_closed(th):
    w = wof(th)
    return mp.sqrt(3)*w/2/mp.sin(mp.pi*w/3)
def ltZ_prod(th):
    return mp.nprod(lambda m: 1/(1+th/b(3*m)), [1, mp.inf])
lt = {}
for th in [mp.mpf(2), mp.mpc(1, 3)]:
    lt[mp.nstr(th, 3)] = {'M_1 closed form': mp.nstr(lt1_closed(th), 25), 'M_1 product': mp.nstr(lt1_prod(th), 25),
                          'M_1 quadrature': mp.nstr(lt1_quad(th), 25),
                          'S closed form': mp.nstr(ltS_closed(th), 25), 'S product': mp.nstr(ltS_prod(th), 25),
                          'Z closed form': mp.nstr(ltZ_closed(th), 25), 'Z product': mp.nstr(ltZ_prod(th), 25),
                          '(M_1 times Z) / S - 1': mp.nstr(lt1_closed(th)*ltZ_closed(th)/ltS_closed(th) - 1, 5)}
put('laplace_M1', lt)
put('laplace_exponent_M1_at_2', mp.nstr(-mp.log(mp.re(lt1_closed(2))), 25))
# zeros of the transform at theta = -b_{6k}: the closed form there, and its poles at -b_n, n = +-1 mod 6, n >= 5
put('closed_form_at_-b6_-b12', [mp.nstr(abs(lt1_closed(-b(6))), 5), mp.nstr(abs(lt1_closed(-b(12))), 5)])
# characteristic function at u = 5
u = mp.mpf(5)
phi_closed = lt1_closed(mp.mpc(0, -u))
phi_quad = 1 - 1j*u*mp.quad(lambda y: mp.exp(1j*u*y)*F1m1(y), [mp.mpf(x)/4 for x in range(0, 81)])
thorin_sum = mp.nsum(lambda k: mp.log(1+u**2/b(6*k-1)**2) + mp.log(1+u**2/b(6*k+1)**2) - mp.log(1+u**2/b(6*k)**2), [1, mp.inf])
# Bondesson function N_1: +1 at b_n for n = +-1 mod 6, n >= 5, -1 at b_{6k}; integral of N_1(z) 2u^2/(z(z^2+u^2)) dz
KB = 20000
atoms = sorted([(b(n), 1) for n in range(5, 6*KB+2) if n % 6 in (1, 5)] + [(b(6*k), -1) for k in range(1, KB+1)])
G = lambda z: mp.log(z**2/(z**2+u**2))
Nval = 0; acc = mp.mpf(0); Nmin = 0
for j in range(len(atoms)-1):
    Nval += atoms[j][1]; Nmin = min(Nmin, Nval)
    acc += Nval*(G(atoms[j+1][0]) - G(atoms[j][0]))
Nval += atoms[-1][1]
Z0 = atoms[-1][0]
# beyond the last atom, N_1(z) = sqrt(3 z/pi)/3 + O(1)
acc += mp.quad(lambda z: mp.sqrt(3*z/mp.pi)/3*2*u**2/(z*(z**2+u**2)), [Z0, mp.inf])
put('char_fn_M1_u5', {'-log|phi|^2, closed form': mp.nstr(-2*mp.log(abs(phi_closed)), 20),
                      '-log|phi|^2, quadrature of the law': mp.nstr(-2*mp.log(abs(phi_quad)), 20),
                      'sum over the Thorin atoms': mp.nstr(thorin_sum, 20),
                      'integral of the Bondesson function': mp.nstr(acc, 20),
                      '|phi(5)|^-2': mp.nstr(1/abs(phi_closed)**2, 20)})
put('bondesson_N1', {'min over the first %d atoms' % len(atoms): Nmin, 'N_1 at the last atom': Nval,
                     'asymptote sqrt(3z/pi)/3 there': mp.nstr(mp.sqrt(3*Z0/mp.pi)/3, 12)})
# Levy density: x nu_1(x) = sum_k (e^{-b_{6k-1}x} + e^{-b_{6k+1}x} - e^{-b_{6k}x}), and the Poisson summation form
def xnu_series(x):
    return mp.fsum(mp.exp(-b(6*k-1)*x) + mp.exp(-b(6*k+1)*x) - mp.exp(-b(6*k)*x) for k in range(1, 201))
def xnu_theta(x):
    th = lambda al: mp.fsum(mp.exp(-12*mp.pi*(m+al)**2/x) for m in range(-40, 41))
    return mp.exp(mp.pi*x/12)/2*(1 + mp.sqrt(12/x)*(th(0)/6 - th(mp.mpf(1)/2)/2 - 2*th(mp.mpf(1)/3)/3)) - 1
put('levy_density_M1', {mp.nstr(x, 3): {'series': mp.nstr(xnu_series(x), 20), 'Poisson form': mp.nstr(xnu_theta(x), 20),
                                         '2 sqrt3 x^{3/2} nu_1(x)': mp.nstr(2*mp.sqrt(3)*mp.sqrt(x)*xnu_series(x), 10)}
                        for x in [mp.mpf('0.01'), mp.mpf('0.1'), mp.mpf(1), mp.mpf(3)]})

# 4. the Laplace transform of M*. With c = 0.3, psi*(theta) = theta [ int_0^c e^{-theta y} F(y) dy
#    + sum_n tau(n) e^{-(theta + 2 pi (n-1)) c} / (theta + 2 pi (n-1)) ], meromorphic, with simple poles at -2 pi m.
#    The integral is taken over [0.02, c] (F(y) < 1e-110 below 0.02) by Gauss Legendre, 8 panels of 96 nodes.
CUT = mp.mpf('0.3'); Y0 = mp.mpf('0.02')
gl = mp.calculus.quadrature.GaussLegendre(mp.mp)
base = gl.calc_nodes(6, mp.mp.prec)
PSI_NODES = []
hp = (CUT - Y0)/8
for k in range(8):
    lo = Y0 + k*hp
    for x, w in base:
        y = lo + hp/2*(x+1)
        PSI_NODES.append((y, hp/2*w*F(y)))
put('F(0.02)', mp.nstr(F(Y0), 5))
def psi_star(th):
    th = mp.mpc(th)
    E = mp.fsum(w*mp.exp(-th*y) for y, w in PSI_NODES)
    S = mp.fsum(tau[n]*mp.exp(-(th + 2*mp.pi*(n-1))*CUT)/(th + 2*mp.pi*(n-1)) for n in range(1, 121))
    return th*(E + S)
two = mp.mpf(2)
put('psi_star_at_2', {'cut form': mp.nstr(mp.re(psi_star(two)), 25),
                      'quadrature 2 int e^{-2y} F(y) dy': mp.nstr(2*mp.quad(lambda y: mp.exp(-2*y)*F(y), BRK + [48, 64]), 25)})
put('residues_at_-2pi_m_over_pi', {m: int(-2*m*tau[m+1]) for m in range(1, 9)})
th0 = mp.findroot(lambda x: mp.re(psi_star(x)), -11*mp.pi)
put('zero_of_psi_star', {'theta': mp.nstr(th0, 25), 'theta/(2 pi)': mp.nstr(th0/(2*mp.pi), 15), '|psi*(theta)|': mp.nstr(abs(psi_star(th0)), 3)})
# sign changes of psi* on the real axis between its poles, on (-15 pi, -pi/10); psi* > 0 on (-2 pi, oo)
ends = [-15*mp.pi] + [-2*mp.pi*m for m in range(7, 0, -1)] + [-mp.pi/10]
sc = {}
for lo, hi in zip(ends[:-1], ends[1:]):
    xs = [lo + (hi-lo)*(j + mp.mpf(1)/2)/400 for j in range(400)]
    vals = [mp.re(psi_star(x)) for x in xs]
    sc['(%s pi, %s pi)' % (mp.nstr(lo/mp.pi, 4), mp.nstr(hi/mp.pi, 4))] = sum(1 for v0, v1 in zip(vals[:-1], vals[1:]) if v0*v1 < 0)
put('real_sign_changes', sc)
# argument principle on the rectangle |Re theta| <= 15 pi, |Im theta| <= HT, which contains the 7 poles -2 pi m, m <= 7:
# zeros = winding number + 7. The change of arg psi* along each edge is summed over steps of at most h0, halving a step
# until the change over it is below 0.3.
def darg(f, z0, z1, f0, f1, depth=0):
    d = mp.im(mp.log(f1/f0))
    if abs(d) < mp.mpf('0.3') or depth > 30:
        return d, depth > 30
    zm = (z0+z1)/2; fm = f(zm)
    d0, e0 = darg(f, z0, zm, f0, fm, depth+1); d1, e1 = darg(f, zm, z1, fm, f1, depth+1)
    return d0 + d1, e0 or e1
def winding(corners, h0):
    tot = mp.mpf(0); bad = False
    for p, q in zip(corners, corners[1:] + corners[:1]):
        n = int(mp.ceil(abs(q-p)/h0)); zs = [p + (q-p)*j/n for j in range(n+1)]; fs = [psi_star(z) for z in zs]
        for j in range(n):
            d, e = darg(psi_star, zs[j], zs[j+1], fs[j], fs[j+1]); tot += d; bad = bad or e
    return tot/(2*mp.pi), bad
for HT, h0s in [(20, ['0.1', '0.05']), (60, ['0.1'])]:
    c = [mp.mpc(-15*mp.pi, -HT), mp.mpc(15*mp.pi, -HT), mp.mpc(15*mp.pi, HT), mp.mpc(-15*mp.pi, HT)]
    res = {}
    for h0 in h0s:
        wnd, bad = winding(c, mp.mpf(h0))
        res['step ' + h0] = {'winding': mp.nstr(wnd, 12), 'zeros = winding + 7': int(mp.nint(wnd)) + 7, 'depth limit hit': bad}
    put('argument_principle_|Re|<=15pi_|Im|<=%d' % HT, res)

# 5. Polya's factor and the halves of Hecke's cut
mp.mp.dps = 30
kb = lambda t: mp.re(mp.besselk(mp.mpc(0, t), 2*mp.pi))
put('first_zero_K_it_2pi', mp.nstr(mp.findroot(kb, mp.mpf('9.77')), 20))
ratio = lambda s: abs(A_hecke(1-s))/abs(A_hecke(s))
put('hecke_halves_ratio_t0', {mp.nstr(d, 3): mp.nstr(ratio(mp.mpf(1)/2 + d), 8) for d in [mp.mpf('0.02'), mp.mpf(4)]})
cross = {}
for d in [mp.mpf('0.1'), mp.mpf('0.5'), mp.mpf(1), mp.mpf(2), mp.mpf(4)]:
    r = lambda t: mp.log(ratio(mp.mpc(mp.mpf(1)/2 + d, t)))
    ts = [mp.mpf(j)/20 for j in range(0, 501)]
    vs = [r(t) for t in ts]
    xs = []
    for t1, t2, v1, v2 in zip(ts[:-1], ts[1:], vs[:-1], vs[1:]):
        if v1*v2 < 0:
            lo, hi, vlo = t1, t2, v1
            for _ in range(30):
                mid = (lo+hi)/2; vm = r(mid)
                if vm*vlo < 0: hi = mid
                else: lo, vlo = mid, vm
            xs.append(mp.nstr((lo+hi)/2, 6))
    cross[mp.nstr(d, 3)] = {'sign at t = 0': int(mp.sign(vs[0])), 'crossings': xs}
put('hecke_halves_equal_modulus_ordinates', cross)

out['seconds'] = time.time()-t0
json.dump(out, open('maxlaw.json', 'w'), indent=1)
print('done', round(time.time()-t0, 1))
