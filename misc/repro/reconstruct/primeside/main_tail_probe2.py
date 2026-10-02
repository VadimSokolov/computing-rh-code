"""RECONSTRUCTION AID, NOT THE AUTHORS' CODE.

Second probe of the fields tv and sup of data/results_main.json. main_tail_probe.py showed
that with the tail of code/analyse.py (smooth density from the Gram point T*, quad with
epsabs=1e-16, epsrel=1e-12) the residual is a near constant offset and tv/eps = 1.05e-9,
against 1.73e-9 archived. Here the same tail integrals are evaluated with scipy.integrate.quad
at its default tolerances (epsabs = epsrel = 1.49e-8), which let a single Gauss Kronrod
panel through, and with the tail started at gamma_1700 plus the boundary term
-K(theta, gamma_1700) (N(gamma_1700) - theta(gamma_1700)/pi - 1), N(gamma_1700) = 1700.
Usage: python main_tail_probe2.py BOOK_results_main.json   (reads main_EPS.npz)
"""
import json, sys
import numpy as np, mpmath as mp
from scipy import integrate
from scipy.interpolate import CubicSpline

book = json.load(open(sys.argv[1]))
rows = {r['eps']: r for r in book['rows']}
trap = integrate.trapezoid
dens = {'log': lambda t: np.log(t / (2 * np.pi)) / (2 * np.pi),
        'log48': lambda t: (0.5 * np.log(t / (2 * np.pi)) + 1 / (48 * t * t)) / np.pi}
mp.mp.dps = 30
out = {}
for e, r in rows.items():
    d = np.load(f'main_{e}.npz'); th = d['th']; rho = d['F'] / np.pi; R = d['Rld']
    Ts, g1700 = float(d['Tstar']), float(d['g1700'])
    Nsm = float(mp.siegeltheta(g1700) / mp.pi + 1)
    c = np.linspace(0, 58, 233)
    print(f"eps {e}: archived tv/eps {r['tv'] / e:.4e}, sup {r['sup']:.4e}")
    res = {}
    for dn, n in dens.items():
        for tol in ['default', 'tight']:
            kw = {} if tol == 'default' else dict(epsabs=1e-16, epsrel=1e-12, limit=200)
            for start in ['Tstar', 'g1700_boundary']:
                T0 = Ts if start == 'Tstar' else g1700
                v = []
                for t in c:
                    K = lambda g: (e / (e ** 2 + (t - g) ** 2) + e / (e ** 2 + (t + g) ** 2)) / np.pi
                    q, _ = integrate.quad(lambda g: K(g) * n(g), T0, np.inf, **kw)
                    if start == 'g1700_boundary':
                        q -= K(g1700) * (1700 - Nsm)
                    v.append(q)
                D = rho - R - CubicSpline(c, v)(th)
                name = f'{start}_{dn}_{tol}'
                res[name] = dict(tv_over_eps=float(trap(np.abs(D), th) / e), sup=float(np.abs(D).max()))
                print(f"   {name:30s} tv/eps {res[name]['tv_over_eps']:.4e}  sup {res[name]['sup']:.4e}")
    out[e] = res
json.dump(out, open('main_tail_probe2.json', 'w'), indent=1)
