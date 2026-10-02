# Part VI reproducibility: archived outputs against a fresh run on Hopper

## Script runs

- eta:arith: {'rc': 0, 'seconds': 3.3, 'host': 'hop070.orc.gmu.edu'}
- eta:base: {'rc': 0, 'seconds': 39.9, 'host': 'hop070.orc.gmu.edu'}
- eta:figs: {'rc': 0, 'seconds': 6.0, 'host': 'hop050.orc.gmu.edu'}
- eta:hit: {'rc': 0, 'seconds': 203.5, 'host': 'hop070.orc.gmu.edu'}
- eta:hit2: {'rc': 0, 'seconds': 12.4, 'host': 'hop070.orc.gmu.edu'}
- eta:zeros12: {'rc': 0, 'seconds': 229.5, 'host': 'hop070.orc.gmu.edu'}
- pm:bnd_alpha_0.6: {'rc': 0, 'seconds': 96.5, 'host': 'hop050.orc.gmu.edu'}
- pm:bnd_alpha_0.75: {'rc': 0, 'seconds': 90.2, 'host': 'hop056.orc.gmu.edu'}
- pm:comp: {'rc': 0, 'seconds': 98.8, 'host': 'hop050.orc.gmu.edu'}
- pm:comp2: {'rc': 0, 'seconds': 94.7, 'host': 'hop050.orc.gmu.edu'}
- pm:figs: {'rc': 0, 'seconds': 1.8, 'host': 'hop050.orc.gmu.edu'}
- pm:ou: {'rc': 0, 'seconds': 1.2, 'host': 'hop066.orc.gmu.edu'}
- sb:comp: {'rc': 0, 'seconds': 4.8, 'host': 'hop050.orc.gmu.edu'}
- sb:figs: {'rc': 0, 'seconds': 1.7, 'host': 'hop050.orc.gmu.edu'}
- sb:weil: {'rc': 0, 'seconds': 9.9, 'host': 'hop066.orc.gmu.edu'}
- sb:weil_mp: {'rc': 0, 'seconds': 10.8, 'host': 'hop069.orc.gmu.edu'}
- sb:weil_mp30: {'rc': 0, 'seconds': 16.7, 'host': 'hop066.orc.gmu.edu'}
- sm:comp: {'rc': 0, 'seconds': 39.6, 'host': 'hop056.orc.gmu.edu'}
- sm:dens2: {'rc': 0, 'seconds': 28.3, 'host': 'hop066.orc.gmu.edu'}
- sm:figs2: {'rc': 0, 'seconds': 9.6, 'host': 'hop050.orc.gmu.edu'}
- sm:improve: {'rc': 0, 'seconds': 12.6, 'host': 'hop050.orc.gmu.edu'}
- sm:strip: {'rc': 0, 'seconds': 112.3, 'host': 'hop066.orc.gmu.edu'}
- stage: {'rc': 0, 'seconds': 0.1, 'host': 'hop050.orc.gmu.edu'}
- wr:salem_riesz: {'rc': 0, 'seconds': 16.2, 'host': 'hop056.orc.gmu.edu'}
- wr:wiener: {'rc': 0, 'seconds': 46.2, 'host': 'hop050.orc.gmu.edu'}
- xi:below: {'rc': 0, 'seconds': 2.4, 'host': 'hop070.orc.gmu.edu'}
- xi:compute: {'rc': 0, 'seconds': 10.2, 'host': 'hop066.orc.gmu.edu'}
- xi:errformula: {'rc': 0, 'seconds': 0.3, 'host': 'hop070.orc.gmu.edu'}
- xi:figs: {'rc': 0, 'seconds': 5.4, 'host': 'hop050.orc.gmu.edu'}
- xi:figs2: {'rc': 0, 'seconds': 4.5, 'host': 'hop050.orc.gmu.edu'}
- xi:figs3: {'rc': 0, 'seconds': 4.4, 'host': 'hop050.orc.gmu.edu'}
- xi:figs4: {'rc': 0, 'seconds': 4.5, 'host': 'hop050.orc.gmu.edu'}
- xi:figs5: {'rc': 0, 'seconds': 4.4, 'host': 'hop050.orc.gmu.edu'}
- xi:improve: {'rc': 0, 'seconds': 24.1, 'host': 'hop050.orc.gmu.edu'}
- xi:levy: {'rc': 0, 'seconds': 2.9, 'host': 'hop056.orc.gmu.edu'}
- xi:mc: {'rc': 0, 'seconds': 3.4, 'host': 'hop066.orc.gmu.edu'}
- xi:rs1: {'rc': 0, 'seconds': 110.8, 'host': 'hop069.orc.gmu.edu'}
- xi:rs1c: {'rc': 0, 'seconds': 5.0, 'host': 'hop069.orc.gmu.edu'}
- xi:rs1d: {'rc': 0, 'seconds': 11.0, 'host': 'hop069.orc.gmu.edu'}

## Outputs

For each archived file: the number of values compared, the largest absolute and relative differences, and the worst entries (relative difference above 1e-9).

### basepoint_one/eta/arith.json: identical
366 values; max abs diff 0.000e+00; max rel diff 0.000e+00

### basepoint_one/eta/base.json: identical
3609 values; max abs diff 0.000e+00; max rel diff 0.000e+00

### basepoint_one/eta/g12.npy: identical
297 values; max abs diff 0.000e+00; max rel diff 0.000e+00

### basepoint_one/eta/hit.json: identical
9242 values; max abs diff 0.000e+00; max rel diff 0.000e+00

### basepoint_one/process_map/alpha1.json: identical
6116 values; max abs diff 0.000e+00; max rel diff 0.000e+00

### basepoint_one/process_map/bnd_0.6.json: agrees to 1e-9
4801 values; max abs diff 1.772e-10; max rel diff 5.309e-10

### basepoint_one/process_map/bnd_0.75.json: agrees to 1e-9
4801 values; max abs diff 6.092e-11; max rel diff 2.386e-10

### basepoint_one/process_map/comp.json: identical
56 values; max abs diff 0.000e+00; max rel diff 0.000e+00

### basepoint_one/process_map/comp2.json: identical
1064 values; max abs diff 0.000e+00; max rel diff 0.000e+00

### basepoint_one/process_map/ou.json: identical
523 values; max abs diff 0.000e+00; max rel diff 0.000e+00

### basepoint_one/process_map/rs1.json: agrees to 1e-9
7997 values; max abs diff 3.883e-11; max rel diff 2.014e-10

### basepoint_one/process_map/volterra.json: agrees to 1e-9
5805 values; max abs diff 5.578e-13; max rel diff 2.248e-13

### basepoint_one/scale_mixture/comp.json: identical
2454 values; max abs diff 0.000e+00; max rel diff 0.000e+00

### basepoint_one/scale_mixture/improve.json: DIFFERS
981 values; max abs diff 2.855e-08; max rel diff 1.000e+00
- /dens/0.5/f: abs 1.502e-08, rel 1.000e+00 (161 values)
- /dens/0.75/f: abs 2.176e-08, rel 1.743e-01 (161 values)
- /dens/1.0/f: abs 2.855e-08, rel 1.741e-01 (161 values)

### basepoint_one/scale_mixture/strip.json: identical
765 values; max abs diff 0.000e+00; max rel diff 0.000e+00

### basepoint_one/wiener/riesz_explicit.json: NOT PRODUCED by the run

### basepoint_one/wiener/salem_riesz.json: identical
28912 values; max abs diff 0.000e+00; max rel diff 0.000e+00

### basepoint_one/wiener/wiener.json: DIFFERS
121 values; max abs diff 2.168e-19; max rel diff 1.353e-16
- /amp: missing in run
- /dxi: missing in run

### basepoint_one/xi/alpha1.json: identical
6116 values; max abs diff 0.000e+00; max rel diff 0.000e+00

### basepoint_one/xi/below.json: identical
400 values; max abs diff 0.000e+00; max rel diff 0.000e+00

### basepoint_one/xi/errformula.json: agrees to 1e-9
240 values; max abs diff 1.694e-21; max rel diff 1.163e-16

### basepoint_one/xi/improve.json: identical
4933 values; max abs diff 0.000e+00; max rel diff 0.000e+00

### basepoint_one/xi/levy.json: identical
107 values; max abs diff 0.000e+00; max rel diff 0.000e+00

### basepoint_one/xi/mc.json: identical
34 values; max abs diff 0.000e+00; max rel diff 0.000e+00

### basepoint_one/xi/rs1.json: agrees to 1e-9
7997 values; max abs diff 3.883e-11; max rel diff 2.014e-10

### basepoint_one/xi/rs1d.json: agrees to 1e-9
874 values; max abs diff 1.645e-15; max rel diff 1.039e-13

### basepoint_one/xi/surv.json: identical
120 values; max abs diff 0.000e+00; max rel diff 0.000e+00

### sato_bondesson/comp.json: agrees to 1e-9
2275 values; max abs diff 8.882e-16; max rel diff 2.215e-16

### sato_bondesson/improve.json: identical
4933 values; max abs diff 0.000e+00; max rel diff 0.000e+00

### sato_bondesson/rs1.json: agrees to 1e-9
7997 values; max abs diff 3.883e-11; max rel diff 2.014e-10


## Follow-up (26 September 2026)

- basepoint_one/wiener/wiener.json: the archived file holds `amp` (the values of |xi(alpha+i gamma_1)| at alpha = 0.5, 0.6, 0.75, 1) and `dxi` (|xi'(1/2+i gamma_1)|), which Section wr:sec:wiener quotes but the archived wiener.py did not compute. wiener.py now computes them; a rerun on Hopper (job 1226949) gives 0.00013831870802157972, 0.00034641150135522304, 0.0006972274929721557 and 0.001382719089216254, equal to the archived values to about 1e-15.
- basepoint_one/scale_mixture/improve.json: the densities agree to 2.9e-8 in absolute terms (relative 1.5e-5 where the density exceeds 1e-3); the large relative differences are in the far tail, where the density is below 1e-7. The printed values are unaffected.
- basepoint_one/wiener/riesz_explicit.json has no script in the archive (the README says the computation was inline); the same file ships with riesz_program/data.

## Rerun of wiener.py after its edit (2026-09-26, 20:45 EDT)

basepoint_one/wiener/wiener.py was edited at 12:12, after Chapter ch:strip was written, to compute the fields amp and dxi that the archived wiener.json already held. The edited script was rerun on an Intel node of Hopper (job 1240391, /scratch/vsokolov/rh_book_repro/part6/wiener_rerun): its wiener.json has the same 126 fields as the archived one and agrees with it to a relative difference of 1.1e-15.
