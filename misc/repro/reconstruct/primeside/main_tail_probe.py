"""RECONSTRUCTION AID, NOT THE AUTHORS' CODE.

Probe of the two fields of data/results_main.json that main_analyse.py does not reproduce,
tv and sup: the total variation and the largest difference on [0, 58] between the Polya
density and the Poisson sum over 1700 ordinates plus a smooth tail. Both are floors of
order 1e-9 eps and 1e-11, set by the residual of the tail model. The probe prints, for each
tail model of main_compute.py, tv/eps and sup next to the archived values, and then asks
which start T0 of the tail (density log(t/2pi)/(2pi)) would reproduce the archived tv at
each resolution, using tail(T0) = tail(T*) - int_{T*}^{T0} K(theta, t) n(t) dt evaluated by quad.
Usage: python main_tail_probe.py BOOK_results_main.json    (reads main_EPS.npz)
"""
import json, sys
import numpy as np
from scipy import integrate, optimize

book = json.load(open(sys.argv[1]))
rows = {r['eps']: r for r in book['rows']}
trap = integrate.trapezoid
out = {}
for e, r in rows.items():
    d = np.load(f'main_{e}.npz'); th = d['th']; rho = d['F'] / np.pi; R = d['Rld']
    print(f"eps {e}: archived tv/eps {r['tv'] / e:.4e}, sup {r['sup']:.4e}")
    res = {}
    for tn in sorted(k for k in d.files if k.startswith('tail_')):
        D = rho - R - d[tn]
        res[tn] = dict(tv_over_eps=float(trap(np.abs(D), th) / e), sup=float(np.abs(D).max()), mean_D=float(D.mean()))
        print(f"   {tn:22s} tv/eps {res[tn]['tv_over_eps']:.4e}  sup {res[tn]['sup']:.4e}  mean D {res[tn]['mean_D']:+.3e}")
    Ts = float(d['Tstar']); n = lambda t: np.log(t / (2 * np.pi)) / (2 * np.pi)
    base = rho - R - d['tail_Tstar_log']
    c = np.linspace(0, 58, 59)

    def tv_at(T0):
        v = []
        for t in c:
            q, _ = integrate.quad(lambda g: (e / (e ** 2 + (t - g) ** 2) + e / (e ** 2 + (t + g) ** 2)) / np.pi * n(g),
                                  Ts, T0, epsabs=1e-20, epsrel=1e-13)
            v.append(q)
        extra = np.interp(th, c, v)
        return float(trap(np.abs(base + extra), th))   # tail(T0) = tail(T*) - extra

    grid = Ts + np.linspace(-2e-3, 2e-3, 81)
    tvs = np.array([tv_at(T0) for T0 in grid])
    k = int(np.argmin(np.abs(tvs - r['tv'])))
    j = int(np.argmin(tvs))
    print(f"   start T0 - T* giving the archived tv: {grid[k] - Ts:+.2e} (tv/eps there {tvs[k] / e:.4e}); "
          f"minimum of tv over T0: {tvs[j] / e:.4e} at T0 - T* = {grid[j] - Ts:+.2e}")
    out[e] = dict(models=res, T0_minus_Tstar_matching=float(grid[k] - Ts), tv_min_over_eps=float(tvs[j] / e),
                  T0_minus_Tstar_at_min=float(grid[j] - Ts))
json.dump(out, open('main_tail_probe.json', 'w'), indent=1)
