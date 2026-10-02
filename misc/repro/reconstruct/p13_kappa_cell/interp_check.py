# NEW SCRIPT, not in the archive: a diagnosis of the column "cell TV" of Table tab:tv (ch/ch07.tex), which
# code/analyse.py computes as tv_cells in data/analysis.json by linear interpolation (np.interp) of the
# cumulative mass M_eps = arg xi(1/2+eps+i theta)/pi of code/part2.py (files eps_*.npz) at the cell boundaries.
# kappa_cell.py (this folder) found the exact d_cell(eps) below the archived tv_cells by 2.5e-3 at eps = 1 and
# 5.0e-4 at eps = 0.5. This script evaluates M_eps exactly (Arb, 128 bits) at the two grid points on either
# side of every cell boundary and at the boundary itself, M_eps(theta) = N(theta) + phi(theta)/pi with
# phi(theta) = arg xi(1/2+eps+i theta) - arg xi(1/2+i theta) continued along the horizontal segment (the proof
# of prop:mass), and separates the error of the grid values from the error of the interpolation. It also tries
# cubic Hermite interpolation with the density F/pi that part2.py saves beside the argument (a one line fix).
# Inputs (read only): eps_*.npz and part1.json of the earlier rerun in /scratch/vsokolov/rh_book_repro/run,
# data/analysis.json in /scratch/vsokolov/rh_book_repro/book_data. Output: interp_check.json.
import json, os, time
import numpy as np
from multiprocessing import Pool
from scipy.interpolate import CubicHermiteSpline
from flint import arb, acb, ctx

RUN = '/scratch/vsokolov/rh_book_repro/run/'
BOOK = '/scratch/vsokolov/rh_book_repro/book_data/'
EPS = [1.0, 0.5, 0.25, 0.1, 0.05, 0.02, 0.01, 0.005]
T0 = time.time()


def phi_at(args):
    """phi(theta, m/1000) by steps of 0.001 in sigma (principal argument of each step of zeta)."""
    bf, m = args
    ctx.prec = 128
    b = arb(bf)
    s0 = acb(arb(1) / 2, b)
    zprev, acc = s0.zeta(), arb(0)
    for j in range(1, m + 1):
        s = acb(arb(1) / 2 + arb(j) / 1000, b)
        z = s.zeta()
        acc += (z / zprev).arg()
        zprev = z
    v = (s.log().imag + (s - 1).log().imag + (s / 2).lgamma().imag
         - s0.log().imag - (s0 - 1).log().imag - (s0 / 2).lgamma().imag + acc)
    return float((v / arb.pi()).mid())


if __name__ == '__main__':
    P1 = json.load(open(RUN + 'part1.json'))
    AN = json.load(open(BOOK + 'analysis.json'))
    near_hp = np.array([np.longdouble(s) for s in P1['ref']])[:80]
    zf = near_hp.astype(float)
    cells = np.concatenate([[0.0], [float((near_hp[k] + near_hp[k + 1]) / 2) for k in range(78)], [200.0]])
    jobs, meta = [], []
    for eps in EPS:
        d = np.load(RUN + f'eps_{eps}.npz'); th = d['theta']
        m = int(round(eps * 1000))
        for k in range(1, 80):
            i = int(np.searchsorted(th, cells[k], side='right')) - 1
            pts = [float(cells[k])] + ([float(th[i]), float(th[i + 1])] if i + 1 < len(th) else [float(th[i])])
            for p in pts:
                jobs.append((p, m)); meta.append((eps, k, p, i))
    with Pool(int(os.environ.get('SLURM_CPUS_PER_TASK', '8'))) as pool:
        vals = pool.map(phi_at, jobs, chunksize=4)
    ex = {(e, float(p)): v for (e, k, p, i), v in zip(meta, vals)}
    out = dict(note='exact M_eps at the cell boundaries and at the grid points around them', rows=[])
    for eps in EPS:
        d = np.load(RUN + f'eps_{eps}.npz'); th = d['theta']; M = d['arg'] / np.pi
        Mc_interp = np.interp(cells, th, M)
        # candidate one line fix for code/analyse.py: cubic Hermite interpolation with the saved derivative
        # M' = rho_eps = F/pi (part2.py saves F = pi rho_eps = d arg/d theta beside arg)
        Mc_herm = CubicHermiteSpline(th, M, d['F'] / np.pi)(cells)
        Mc_exact = np.zeros(80); grid_err = 0.0
        for k in range(1, 80):
            Mc_exact[k] = np.searchsorted(zf, cells[k]) + ex[(eps, float(cells[k]))]
            i = int(np.searchsorted(th, cells[k], side='right')) - 1
            for j in (i, i + 1):
                if j < len(th):
                    ex_j = np.searchsorted(zf, th[j]) + ex[(eps, float(th[j]))]
                    grid_err = max(grid_err, abs(M[j] - ex_j))
        tv_interp = float(np.sum(np.abs(np.diff(Mc_interp) - 1)))
        tv_exact = float(np.sum(np.abs(np.diff(Mc_exact) - 1)))
        tv_herm = float(np.sum(np.abs(np.diff(Mc_herm) - 1)))
        a = [r for r in AN['rows'] if r['eps'] == eps][0]
        row = dict(eps=eps, grid_step=float(th[1] - th[0]), npts=len(th),
                   max_abs_grid_value_error=float(grid_err),
                   max_abs_interp_error_at_cells=float(np.max(np.abs(Mc_interp - Mc_exact))),
                   tv_cells_interp_rerun=tv_interp, tv_cells_exact=tv_exact, tv_cells_archived=a['tv_cells'],
                   tv_cells_hermite=tv_herm, max_abs_hermite_error_at_cells=float(np.max(np.abs(Mc_herm - Mc_exact))),
                   cell_min_hermite=float(np.min(np.diff(Mc_herm))), cell_max_hermite=float(np.max(np.diff(Mc_herm))),
                   cell_min_exact=float(1 + np.min(np.diff(Mc_exact) - 1)),
                   cell_max_exact=float(1 + np.max(np.diff(Mc_exact) - 1)))
        out['rows'].append(row)
        print('eps %-6g step %.5f: grid values err %.1e, interpolation err at cells %.1e; tv_cells interp %.10f,'
              ' archived %.10f, exact %.10f (archived-exact %+.2e); Hermite %.10f (err at cells %.1e, tv diff %+.1e)' % (
                  eps, row['grid_step'], grid_err, row['max_abs_interp_error_at_cells'], tv_interp, a['tv_cells'],
                  tv_exact, a['tv_cells'] - tv_exact, tv_herm, row['max_abs_hermite_error_at_cells'],
                  tv_herm - tv_exact), flush=True)
    json.dump(out, open('interp_check.json', 'w'), indent=1)
    print('done in %.1fs' % (time.time() - T0))
