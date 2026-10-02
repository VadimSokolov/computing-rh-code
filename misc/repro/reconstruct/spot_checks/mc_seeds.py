#!/usr/bin/env python3
# NOT an author script. Audit helper for the three Monte Carlo lines that
# verification_reciprocal_zeta/spot_checks.py prints for Proposition prop:hmfalse
# (they are in the regenerated spot_checks_report.md, not in the archived one).
#
# It repeats the script's sampling code verbatim (random.seed, two random.expovariate(1.0)
# per draw, N = 200000, H = 2xy/(x+y), P(H>h) counted against an mpf threshold) for the
# script's seed 20260920 and for 400 other seeds, and compares with the closed form
# P(H>h) = h e^{-h} K_1(h) that the script also prints. It answers two questions:
# does the seeded run reproduce exactly, and how far from the closed form should a
# run with another seed or another generator be expected to land.
import json
import os
import random
from multiprocessing import Pool

from mpmath import mp, mpf, besselk, exp, nstr, sqrt

mp.dps = 25
N = 200000
HS = ['0.1', '1', '2']
SEED = 20260920
OTHER = list(range(1, 401))


def closed(h):
    h = mpf(h)
    return h * exp(-h) * besselk(1, h)


def run(seed):
    random.seed(seed)
    draws = sorted(2 * x * y / (x + y) for x, y in
                   ((random.expovariate(1.0), random.expovariate(1.0)) for _ in range(N)))
    out = {}
    for h in HS:
        hv = mpf(h)
        out[h] = mpf(sum(1 for v in draws if v > hv)) / N
    return seed, out


if __name__ == '__main__':
    ncpu = int(os.environ.get('SLURM_CPUS_PER_TASK', '8'))
    with Pool(ncpu) as pool:
        results = dict(pool.map(run, [SEED] + OTHER))
    res = {'N': N, 'seed': SEED, 'other_seeds': f'{OTHER[0]}..{OTHER[-1]}', 'h': {}}
    for h in HS:
        p = closed(h)
        se = sqrt(p * (1 - p) / N)
        est = [results[s][h] for s in OTHER]
        mean = sum(est) / len(est)
        sd = sqrt(sum((e - mean) ** 2 for e in est) / (len(est) - 1))
        z0 = (results[SEED][h] - p) / se
        zs = [abs((e - p) / se) for e in est]
        res['h'][h] = {
            'closed_form': nstr(p, 10),
            'script_seed_estimate_printed_as_in_report': nstr(results[SEED][h], 6),
            'script_seed_z': nstr(z0, 4),
            'standard_error_binomial': nstr(se, 4),
            'other_seeds_mean': nstr(mean, 8),
            'other_seeds_sd': nstr(sd, 4),
            'mean_minus_closed_in_se_of_mean': nstr((mean - p) / (sd / sqrt(len(est))), 3),
            'fraction_other_seeds_with_abs_z_at_least_script_seed': nstr(mpf(sum(1 for z in zs if z >= abs(z0))) / len(zs), 4),
            'max_abs_diff_other_seeds': nstr(max(abs(e - p) for e in est), 4),
            'script_tolerance_0.004_in_se': nstr(mpf('0.004') / se, 3),
        }
    with open('mc_seeds.json', 'w') as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1))
