# Chapter 13, Section ch:fmexp: the size biased laws B_q of B_0 = (pi/2) S_2 and their Fourati phases on the cut.
# B_0 has density f_0(y) = sum_m (4 pi^2 m^4 y - 6 pi m^2) e^{-pi m^2 y} = 4 y omega''(y) + 6 omega'(y), omega(y) = sum_m e^{-pi m^2 y},
# and E B_0^s = 2 xi(2s) (bp:eq:bpy). For 0 < q < 1, B_q has density y^q f_0(y)/E B_0^q, so E B_q^s = xi(2s + 2q)/xi(2q) and
# B_{1/4} = e^{2X}.
# Lambda_q(lam) = int_0^inf y^q f_0(y) e^{-lam y} dy is continued to lam = -r + i0, 0 < r < pi (M+1)^2, by splitting f_0 into its
# first M terms, whose transforms are gamma functions, and a remainder R_M = O(e^{-pi (M+1)^2 y}), bounded near 0, whose transform
# is a real convergent integral; M = 5 covers r <= 90 < 36 pi.
# The Fourati phase theta_q(r) = -arg Lambda_q(-r + i0)/pi, continued in r, rises by 2 + q at each branch point pi n^2;
# np.unwrap sees only q of that rise, hence the + 2 N(r), N(r) = #{n : pi n^2 < r}.
# Checks: (a) theta_q is 2n + q just after pi n^2 and 2n just before pi (n+1)^2, and nonincreasing in between (fm:thm:sbnotggc);
# (b) theta_G = (2+q) N - theta_q is nondecreasing and equal to (n-1) q at pi n^2 (fm:thm:factor);
# (c) k_q(x) = (2+q) omega(x) - x int_0^inf e^{-tx} theta_G(t) dt is positive and decreasing on [1e-3, 3] (fm:thm:sd), with
#     theta_G(t) = q (sqrt(t/pi) - 1) beyond the scan, the value that a mass q on each gap gives at the points pi n^2;
# (d) for q = 1/4: -log|E e^{iyB}|^2 three ways, the direct Fourier integral of the density, the Thorin form (the atoms 9/4 at
#     pi n^2 summed in closed form, minus the smooth part from theta_G) and the theta Levy form int 2(1 - cos yx) k(x) dx/x;
#     and E G_{1/4} = int theta_G(t) t^{-2} dt against 3 pi/8 - xi(5/2)/xi(1/2).
# Usage: python3 fm_sizebias.py   Writes fm_sizebias.json.
import json, os, time
from multiprocessing import Pool
import numpy as np, mpmath as mp
from scipy.integrate import quad
from scipy.special import gammaincc, gamma as Gam

M = 5
QS = [0.1, 0.25, 0.5, 0.75, 0.9]
RS = np.linspace(0.05, 90, 1800)
t0 = time.time()


def S(y, lo, hi):
    return mp.fsum((4*mp.pi**2*m**4*y - 6*mp.pi*m**2)*mp.exp(-mp.pi*m*m*y) for m in range(lo, hi))


def f0(y):
    if y >= 1: return S(y, 1, 30)
    return y**(-mp.mpf(5)/2)*S(1/y, 1, 30)            # f_0(y) = y^{-5/2} f_0(1/y)


def RM(y):
    if y > mp.mpf('0.3'): return S(y, M + 1, 60)
    return f0(y) - S(y, 1, M + 1)


def Lam(args):
    q, r = args
    mp.mp.dps = 30 + int(r/4)
    q = mp.mpf(q); r = mp.mpf(r)
    lam = mp.mpc(-r, mp.mpf(10)**(-mp.mp.dps + 5))
    fin = mp.fsum(4*mp.pi**2*m**4*mp.gamma(q + 2)*(mp.pi*m*m + lam)**(-q - 2)
                  - 6*mp.pi*m**2*mp.gamma(q + 1)*(mp.pi*m*m + lam)**(-q - 1) for m in range(1, M + 1))
    rem = mp.quad(lambda y: y**q*RM(y)*mp.exp(r*y), [0, 0.05, 0.3, 1, 3, mp.inf])
    return complex(fin + rem)


def omega(x):
    if x >= 1: return float(sum(np.exp(-np.pi*n*n*x) for n in range(1, 12)))
    return 0.5*(x**-0.5*(1 + 2*omega(1/x)) - 1)       # theta(x) = x^{-1/2} theta(1/x)


def analyse(q, v):
    N = np.floor(np.sqrt(RS/np.pi))
    th = -np.unwrap(np.angle(v))/np.pi + 2*N
    thG = (2 + q)*N - th
    gaps = []
    for n in range(1, 5):
        i0 = np.searchsorted(RS, np.pi*n*n); i1 = np.searchsorted(RS, np.pi*(n + 1)**2) - 1
        seg = th[i0:i1 + 1]
        gaps.append(dict(n=n, after=float(th[i0]), before=float(th[i1]), max_increment=float(np.max(np.diff(seg))),
                         thetaG_before_next=float(thG[i1])))
    dG = np.diff(thG[RS > np.pi])
    T = RS[-1]
    tail = lambda t: q*(np.sqrt(t/np.pi) - 1)
    def lap_thG(x):                                   # int_0^inf e^{-tx} theta_G(t) dt
        head = np.trapezoid(np.exp(-RS*x)*thG, RS)
        # int_T^inf e^{-tx} q (sqrt(t/pi) - 1) dt = q [pi^{-1/2} x^{-3/2} Gamma(3/2, Tx) - e^{-Tx}/x]
        tl = q*(np.pi**-0.5*x**-1.5*Gam(1.5)*gammaincc(1.5, T*x) - np.exp(-T*x)/x)
        return head + tl
    xs = np.geomspace(1e-3, 3, 300)
    k = np.array([(2 + q)*omega(x) - x*lap_thG(x) for x in xs])
    EG = np.trapezoid(thG/RS**2, RS) + q*(2/np.sqrt(np.pi*T) - 1/T)
    res = dict(q=q, gaps=gaps, thetaG_min_increment=float(dG.min()), thetaG_at_end=float(thG[-1]),
               k_min=float(k.min()), k_decreasing=bool(np.all(np.diff(k) < 0)),
               subtracted_fraction_at_1e3=float(1 - k[0]/((2 + q)*omega(xs[0]))), q_over_2_plus_q=q/(2 + q),
               EG_thorin=float(EG))
    return res, th, thG, lap_thG


if __name__ == '__main__':
    tasks = [(q, r) for q in QS for r in RS]
    with Pool(int(os.environ.get('SLURM_CPUS_PER_TASK', '8'))) as P:
        vals = P.map(Lam, tasks, chunksize=8)
    V = np.array(vals).reshape(len(QS), len(RS))
    np.save('fm_sizebias_phase.npy', V)
    print('phases computed, %.0f s' % (time.time() - t0), flush=True)
    R = dict(r_grid=[float(RS[0]), float(RS[-1]), len(RS)], M=M, analysis=[])
    for i, q in enumerate(QS):
        res, th, thG, lap = analyse(q, V[i])
        R['analysis'].append(res); print(res, flush=True)
        if q == 0.25:
            mp.mp.dps = 25
            Z = mp.quad(lambda w: w**0.25*f0(w), [0, 0.2, 0.5, 1, 2, 5, mp.inf])
            xi = lambda s: s*(s - 1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
            exact_EG = 3*mp.pi/8 - xi(mp.mpf(5)/2)/xi(mp.mpf(1)/2)
            T = RS[-1]
            tail = lambda t: q*(np.sqrt(t/np.pi) - 1)
            xg = np.geomspace(1e-7, 12, 20000)
            kg = np.array([(2 + q)*omega(x) - x*lap(x) for x in xg])
            rows = []
            for y in [0.5, 1, 2, 4, 8]:
                re = mp.quad(lambda w: w**0.25*f0(w)*mp.cos(y*w), [0, 0.2, 0.5, 1, 2, 3, 5, 8, mp.inf])/Z
                im = mp.quad(lambda w: w**0.25*f0(w)*mp.sin(y*w), [0, 0.2, 0.5, 1, 2, 3, 5, 8, mp.inf])/Z
                direct = float(-mp.log(re*re + im*im))
                s = np.sqrt(2*np.pi*y)
                atoms = 2.25*np.log((np.cosh(s) - np.cos(s))/(2*np.pi*y))
                smooth = np.trapezoid(thG*2*y*y/(RS*(RS**2 + y*y)), RS) + quad(lambda t: tail(t)*2*y*y/(t*(t*t + y*y)), T, np.inf, limit=200)[0]
                levy = np.trapezoid(2*(1 - np.cos(y*xg))*kg/xg, xg)
                rows.append(dict(y=y, direct=direct, thorin=float(atoms - smooth), theta_levy=float(levy), atoms_alone=float(atoms)))
                print(rows[-1], flush=True)
            R['q_quarter'] = dict(normaliser=mp.nstr(Z, 12), two_xi_half=mp.nstr(2*xi(mp.mpf(1)/2), 12), sin2_table=rows,
                                  EG_exact=mp.nstr(exact_EG, 10), EG_thorin=res['EG_thorin'])
    R['seconds'] = time.time() - t0
    json.dump(R, open('fm_sizebias.json', 'w'), indent=1)
