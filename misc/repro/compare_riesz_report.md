# Reproducibility: archived outputs against a fresh run on Hopper (riesz)

## Script runs

- almost: {'rc': 0, 'seconds': 3.3, 'host': 'hop070.orc.gmu.edu'}
- check_large: {'rc': 0, 'seconds': 53.8, 'host': 'hop066.orc.gmu.edu'}
- figs_riesz: {'rc': 0, 'seconds': 6.6, 'host': 'hop066.orc.gmu.edu'}
- figs_strategy: {'rc': 0, 'seconds': 7.2, 'host': 'hop066.orc.gmu.edu'}
- ggc_coef: {'rc': 0, 'seconds': 21.3, 'host': 'hop050.orc.gmu.edu'}
- high: {'rc': 0, 'seconds': 97.6, 'host': 'hop069.orc.gmu.edu'}
- killing: {'rc': 0, 'seconds': 5.2, 'host': 'hop066.orc.gmu.edu'}
- killing2: {'rc': 0, 'seconds': 12.0, 'host': 'hop066.orc.gmu.edu'}
- ks: {'rc': 0, 'seconds': 114.4, 'host': 'hop056.orc.gmu.edu'}
- more2: {'rc': 0, 'seconds': 82.3, 'host': 'hop070.orc.gmu.edu'}
- riesz_large_tw: {'rc': 0, 'seconds': 16.5, 'host': 'hop050.orc.gmu.edu'}
- riesz_thorin: {'rc': 0, 'seconds': 6.5, 'host': 'hop050.orc.gmu.edu'}
- salem_riesz: {'rc': 0, 'seconds': 19.1, 'host': 'hop069.orc.gmu.edu'}
- sens: {'rc': 0, 'seconds': 16.7, 'host': 'hop070.orc.gmu.edu'}
- twisted: {'rc': 0, 'seconds': 66.8, 'host': 'hop056.orc.gmu.edu'}

## Outputs

For each archived file: the number of values compared, the largest absolute and relative differences, and the worst entries (relative difference above 1e-9).

### code/almost.json: identical
37948 values; max abs diff 0.000e+00; max rel diff 0.000e+00

### code/check_large.json: DIFFERS
20 values; max abs diff 5.792e+00; max rel diff 1.656e-01
- [1]/sec: abs 7.126e-01, rel 1.656e-01 (1 values)
- [2]/sec: abs 2.126e+00, rel 1.569e-01 (1 values)
- [3]/sec: abs 5.792e+00, rel 1.371e-01 (1 values)
- [0]/sec: abs 3.076e-03, rel 8.493e-03 (1 values)

### code/ggc_coef.json: identical
45 values; max abs diff 0.000e+00; max rel diff 0.000e+00

### code/high.json: identical
3260 values; max abs diff 0.000e+00; max rel diff 0.000e+00

### code/killing.json: agrees to 1e-9
146 values; max abs diff 3.638e-12; max rel diff 2.254e-16

### code/ks.json: agrees to 1e-9
53 values; max abs diff 9.095e-13; max rel diff 1.729e-16

### code/more.json: identical
117 values; max abs diff 0.000e+00; max rel diff 0.000e+00

### code/riesz_explicit.json: identical
9 values; max abs diff 0.000e+00; max rel diff 0.000e+00

### code/riesz_thorin.json: identical
244 values; max abs diff 0.000e+00; max rel diff 0.000e+00

### code/salem_riesz.json: identical
28912 values; max abs diff 0.000e+00; max rel diff 0.000e+00

### code/sens.json: agrees to 1e-9
289 values; max abs diff 1.289e-13; max rel diff 2.200e-13

### code/twisted.json: identical
1103 values; max abs diff 0.000e+00; max rel diff 0.000e+00


## Constants hard-coded in ks.py

ks.py hard-codes the arithmetic factors a_1, a_2, a_3 of the Keating and Snaith moments (its function ak, which computes them, is never called). misc/repro/ks_ak.py recomputes them on Hopper from the same Euler product over the primes below 2e5 at 30 digits (job 1229074, output misc/repro/ks_ak.json): the hard-coded values agree to relative 2e-14. The truncated product for a_2 differs from the exact value 6/pi^2 by a relative 3.8e-7, so it agrees with 6/pi^2 to six digits, not seven as the note says.

## Rerun of v1min.py (2026-09-26, 20:45 EDT)

riesz_program/code/v1min.py, the computation behind Proposition bp:prop:v1min, was rerun on an Intel node of Hopper (job 1240390, /scratch/vsokolov/rh_book_repro/v1check/rerun): its output agrees with riesz_program/data/v1min.json in all 38 fields, exactly. (Note of 28 September 2026: that v1min.json has since been replaced by the output of the rewritten v1min.py, in which every inequality of the proof is a directed comparison of Arb balls, job 1313130; misc/gpt-review/fixes/basepoint.md.) An independent check with mpmath (misc/repro/v1_grid_mp.py, zeros_mp.py, v1_collect.py; arrays 1240387 and 1240388) is reported below when it finishes.

The independent check finished at 20:55 (collector job 1240389, output misc/repro/v1_check.json). On the grid of step 0.005 on (0, 2100], 420000 points evaluated with mpmath at 25 digits, v_1 exceeds b_1 = 0.0230957090 at every point; the smallest excess is 2.8e-9 at theta = 0.005, as the second order vanishing at the origin predicts (v_1(1) - b_1 = 1.12e-4). On [13.5, 2000] the smallest value is 0.1132554 at theta = 17.495, against 0.1132564 at 17.5 on the coarser grid of v1min.py. The first 1850 zeros from mpmath.zetazero agree with those from Arb to double precision; the largest gap below 2000 is 6.887314, between gamma_1 and gamma_2, so every theta in [13.5, 2000] lies within 3.4437 of an ordinate, below the threshold 4.6259 of the proof.

Note of 28 September 2026: run_riesz.py now also runs riesz_program/code/v1min.py, in the rewritten version whose every inequality is a directed comparison of Arb balls, and explicit_bounds.py, added for the book (the certificate for Theorem ag:thm:explicit and part (3) of Theorem sm:thm:strip). Both run in the rest step. Reruns on Hopper in fresh folders reproduce the archived outputs: v1min.json (job 1323300) in every field but the run time, host and job number, and explicit_bounds.json (job 1318565) in every field but the run time. compare6.py will report those fields as differences.
