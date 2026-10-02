# Chapter 13, the Polya law as a ferromagnetic law: moments, cumulants, Newton ratios, the kernel's convexity,
# the uniform mixture, and the distances of the simple approximations, all in the book's units:
#   Phi(t) = 2 sum_n (2 pi^2 n^4 e^{9t/2} - 3 pi n^2 e^{5t/2}) exp(-pi n^2 e^{2t}),  p = Phi / Xi(0),
#   E e^{i theta X} = Xi(theta)/Xi(0), whose real zeros are the ordinates gamma_k.
# Reads zeros_1700.txt (the first 1700 ordinates to 70 digits); writes fm_polya.json.
import json, sys, time
import numpy as np, mpmath as mp

t0 = time.time()
R = {}
pi = np.pi


def Phi_np(t):
    t = np.abs(np.asarray(t, float)); s = 0.0
    for n in range(1, 8):
        s = s + 2*(2*pi**2*n**4*np.exp(4.5*t) - 3*pi*n**2*np.exp(2.5*t))*np.exp(-pi*n**2*np.exp(2*t))
    return s


mp.mp.dps = 40


def Phi_mp(t):
    t = abs(mp.mpf(t))
    return 2*mp.fsum((2*mp.pi**2*n**4*mp.exp(mp.mpf(9)/2*t) - 3*mp.pi*n**2*mp.exp(mp.mpf(5)/2*t))
                     * mp.exp(-mp.pi*n**2*mp.exp(2*t)) for n in range(1, 12))


def xi(s):
    return s*(s - 1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)


X0 = xi(mp.mpf(1)/2)
cuts = [0, mp.mpf(1)/4, mp.mpf(1)/2, 1, mp.mpf(3)/2, 2, 3, 4]
mass = 2*mp.quad(Phi_mp, cuts)
R['Xi0'] = mp.nstr(X0, 20); R['int_Phi'] = mp.nstr(mass, 20)
print('Xi(0) =', mp.nstr(X0, 20), ' int Phi =', mp.nstr(mass, 20), flush=True)

# ---- moments m_2j and Newton ratios of b_j = m_2j/(2j)!  (j = 0..25)
m = [2*mp.quad(lambda t, j=j: t**(2*j)*Phi_mp(t), cuts)/X0 for j in range(26)]
b = [m[j]/mp.factorial(2*j) for j in range(26)]
newton = []
for j in range(1, 25):
    r = b[j]**2/(b[j - 1]*b[j + 1])
    newton.append(dict(j=j, b=mp.nstr(b[j], 8), ratio=mp.nstr(r, 8), bound=mp.nstr(1 + mp.mpf(1)/j, 8),
                       margin=mp.nstr(r - 1 - mp.mpf(1)/j, 6)))
R['newton'] = newton
print('Newton ratios b_j^2/(b_{j-1} b_{j+1}) against 1+1/j:')
for d in newton:
    print('  ', d['j'], d['b'], d['ratio'], d['bound'], d['margin'])
k2 = m[1]; k4 = m[2] - 3*m[1]**2; k6 = m[3] - 15*m[2]*m[1] + 30*m[1]**3
R['kappa_quad'] = dict(k2=mp.nstr(k2, 12), k4=mp.nstr(k4, 12), k6=mp.nstr(k6, 12))
print('cumulants by quadrature: k2', mp.nstr(k2, 12), ' k4', mp.nstr(k4, 12), ' k6', mp.nstr(k6, 12), flush=True)

# ---- the same cumulants from the Taylor coefficients of log xi(1/2+s) (no hypothesis), and from the zeros under RH
mp.mp.dps = 50
lx = lambda s: mp.log(xi(mp.mpf(1)/2 + s))
tc = mp.taylor(lx, 0, 6)
kt = [tc[2*j]*mp.factorial(2*j) for j in range(4)]
R['kappa_taylor'] = dict(k2=mp.nstr(kt[1], 12), k4=mp.nstr(kt[2], 12), k6=mp.nstr(kt[3], 12))
print('cumulants from log xi:   k2', mp.nstr(kt[1], 12), ' k4', mp.nstr(kt[2], 12), ' k6', mp.nstr(kt[3], 12))
mp.mp.dps = 40
gam = [mp.mpf(l.strip()) for l in open('zeros_1700.txt') if l.strip()]
G = gam[-1]; NG = len(gam)
Ns = lambda t: t/(2*mp.pi)*mp.log(t/(2*mp.pi)) - t/(2*mp.pi) + mp.mpf(7)/8
pz = {}
for mm in (1, 2, 3):
    head = mp.fsum(g**(-2*mm) for g in gam)
    tail = mp.quad(lambda t: (Ns(t) - NG)*2*mm*t**(-2*mm - 1), [G, 2*G, 10*G, mp.inf])
    pz[mm] = head + tail
kz = [None, 2*pz[1], -12*pz[2], 240*pz[3]]
R['kappa_zeros'] = dict(k2=mp.nstr(kz[1], 12), k4=mp.nstr(kz[2], 12), k6=mp.nstr(kz[3], 12),
                        p1=mp.nstr(pz[1], 12), p2=mp.nstr(pz[2], 12), p3=mp.nstr(pz[3], 12))
print('cumulants from 1700 zeros and a smooth tail: k2', mp.nstr(kz[1], 12), ' k4', mp.nstr(kz[2], 12),
      ' k6', mp.nstr(kz[3], 12), flush=True)
var = float(k2)

# ---- Maxwell range of the zero completed laws: Var X >= 3 sum_{k<=K} gamma_k^-2
S = mp.mpf(0); Kmax = None
for K in range(1, 40):
    S += gam[K - 1]**(-2)
    if 3*S > k2 and Kmax is None:
        Kmax = K - 1
        R['maxwell_range'] = dict(K=Kmax, three_sum_K=mp.nstr(3*(S - gam[K - 1]**(-2)), 8),
                                  three_sum_K1=mp.nstr(3*S, 8), var=mp.nstr(k2, 8))
print('Maxwell construction possible for K <=', Kmax, R['maxwell_range'], flush=True)

# ---- the kernel's convexity: F = -log Phi, F'' and F'''' at 0 and the sign of F''' on (0, 4]
F = lambda t: -mp.log(Phi_mp(t))
R['F2_0'] = mp.nstr(mp.diff(F, 0, 2), 10); R['F4_0'] = mp.nstr(mp.diff(F, 0, 4), 10)
f3 = [(float(u), float(mp.diff(F, u, 3))) for u in [mp.mpf(i)/40 for i in range(1, 161)]]
R['F3_min_grid'] = min(f3, key=lambda r: r[1]); R['F3_at'] = {str(u): v for u, v in f3 if u in (0.25, 0.5, 1.0, 2.0, 4.0)}
print("F''(0) =", R['F2_0'], " F''''(0) =", R['F4_0'], " min F''' on (0,4]:", R['F3_min_grid'], flush=True)

# ---- grids for distances
x = np.linspace(-2.4, 2.4, 480001); dx = x[1] - x[0]
p = Phi_np(x)/float(X0)
R['p0'] = float(p[len(x)//2]); R['mass_grid'] = float(np.trapezoid(p, x))
def tv(q): return 0.5*np.trapezoid(np.abs(q - p), x)
def ks(q): return float(np.max(np.abs(np.cumsum(q - p)*dx)))
def norm(q): return q/np.trapezoid(q, x)
def var_of(q): return float(np.trapezoid(x*x*q, x))
gauss = np.exp(-x**2/(2*var))/np.sqrt(2*pi*var)
D = {}
D['gauss'] = dict(tv=tv(gauss), ks=ks(gauss), var=var, at0=float(gauss[len(x)//2]))
Ps = norm(8*pi**2*np.cosh(4.5*x)*np.exp(-2*pi*np.cosh(2*x)))
Pt = norm((8*pi**2*np.cosh(4.5*x) - 12*pi*np.cosh(2.5*x))*np.exp(-2*pi*np.cosh(2*x)))
D['polya_star'] = dict(tv=tv(Ps), ks=ks(Ps), var=var_of(Ps), at0=float(Ps[len(x)//2]))
D['psi'] = dict(tv=tv(Pt), ks=ks(Pt), var=var_of(Pt), at0=float(Pt[len(x)//2]), min=float(Pt.min()))
for lam in (0.2, 0.22):
    q = norm(np.exp(lam*x**2/4)*p); D[f'dbn_{lam}'] = dict(tv=tv(q), ks=ks(q))
def shc(z):                      # sinh(z)/z, equal to 1 at z = 0
    out = np.ones_like(z); nz = z != 0; out[nz] = np.sinh(z[nz])/z[nz]; return out
mult = {'cosh_half': np.cosh(x/2), 'uniform_sum_quarter': shc(np.sqrt(3)/2*x),
        'gauss_lambda_half': np.exp(x**2/8), 'cosh_sqrt01': np.cosh(np.sqrt(0.1)*x),
        'uniform_sum_twelfth': shc(x/2)}
for k, w in mult.items():
    q = norm(w*p); D[k] = dict(tv=tv(q), ks=ks(q))
R['distances'] = D
R['p0_check'] = float(Phi_mp(0)/X0)
for k, v in D.items():
    print('  ', k, {a: float('%.6g' % b) for a, b in v.items()}, flush=True)

# ---- Khinchine: X = U Y with h(y) = -2 y p'(y); E Y = 2 E|X|, E Y^2 = 3 E X^2
def dPhi_np(t):
    t = np.asarray(t, float); s = 0.0
    for n in range(1, 8):
        A = 2*pi**2*n**4*np.exp(4.5*t) - 3*pi*n**2*np.exp(2.5*t)
        dA = 9*pi**2*n**4*np.exp(4.5*t) - 7.5*pi*n**2*np.exp(2.5*t)
        E = np.exp(-pi*n**2*np.exp(2*t))
        s = s + 2*(dA - A*2*pi*n**2*np.exp(2*t))*E
    return s
y = np.linspace(0, 2.4, 240001)
h = -2*y*dPhi_np(y)/float(X0)
R['khinchine'] = dict(h_min=float(h.min()), mass=float(np.trapezoid(h, y)), EY=float(np.trapezoid(y*h, y)),
                      EY2_over_3=float(np.trapezoid(y*y*h, y)/3))
EabsX = 2*mp.quad(lambda t: t*Phi_mp(t), cuts)/X0
R['khinchine']['E_absX'] = float(EabsX); R['khinchine']['EY_exact'] = float(2*EabsX)
print('Khinchine mixing density:', R['khinchine'], flush=True)

# ---- Polya's Phi* and two term Psi: zeros of the transforms on (0.25, 62.5) and in the box |Im| <= 2
mp.mp.dps = 20
Kb = lambda nu: mp.besselk(nu, 2*mp.pi)
Hs = lambda th: Kb((9 + 2j*th)/4) + Kb((9 - 2j*th)/4)
Hp = lambda th: 4*mp.pi**2*(Kb((9 + 2j*th)/4) + Kb((9 - 2j*th)/4)) - 6*mp.pi*(Kb((5 + 2j*th)/4) + Kb((5 - 2j*th)/4))
def real_zeros(f, a, Z, h):
    zz = np.arange(a, Z, h); v = [float(mp.re(f(z))) for z in zz]; out = []
    for i in range(len(zz) - 1):
        if v[i]*v[i + 1] < 0:
            out.append(float(mp.findroot(lambda z: mp.re(f(z)), (zz[i], zz[i + 1]), solver='anderson')))
    return out
def winding(f, a, Z, H, n):
    pts = [complex(a, -H) + (Z - a)*s for s in np.linspace(0, 1, n)]
    pts += [complex(Z, -H + 2*H*s) for s in np.linspace(0, 1, n)[1:]]
    pts += [complex(Z - (Z - a)*s, H) for s in np.linspace(0, 1, n)[1:]]
    pts += [complex(a, H - 2*H*s) for s in np.linspace(0, 1, n)[1:]]
    ang = np.unwrap(np.angle([complex(f(z)) for z in pts]))
    return int(round((ang[-1] - ang[0])/(2*pi)))
Z = 62.5
gz = [float(g) for g in gam if g < Z]
for name, f in (('polya_star', Hs), ('psi', Hp)):
    zs = real_zeros(f, 0.25, Z, 0.025)
    w1 = winding(f, 0.25, Z, 2.0, 1500); w2 = winding(f, 0.25, Z, 2.0, 3000)
    R[name + '_zeros'] = dict(real=zs, n_real=len(zs), box=[w1, w2])
    print(name, 'real zeros', len(zs), [round(z, 4) for z in zs[:6]], ' box counts', w1, w2, flush=True)
R['gamma_below_62.5'] = gz

# ---- finite uniform mixtures: J atoms at the quantiles of Y
Hc = np.cumsum(h); Hc /= Hc[-1]
for J in (5, 20):
    yj = np.interp((np.arange(J) + 0.5)/J, Hc, y)
    f = lambda z, yj=yj: mp.fsum([mp.sin(yy*z)/(yy*z) for yy in yj])/J
    zs = real_zeros(f, 0.25, 100.0, 0.005)
    w1 = winding(f, 0.25, 100.0, 30.0, 3000); w2 = winding(f, 0.25, 100.0, 30.0, 6000)
    R[f'mix{J}'] = dict(y_max=float(yj.max()), n_real=len(zs), box=[w1, w2])
    print(f'J={J}: real zeros in (0.25,100): {len(zs)}; zeros in (0.25,100)x(-30,30): {w1} {w2}', flush=True)

R['seconds'] = time.time() - t0
json.dump(R, open('fm_polya.json', 'w'), indent=1, default=str)
print('done in %.0f s' % R['seconds'])
