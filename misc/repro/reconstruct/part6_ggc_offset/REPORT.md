# Reproducibility report: ggc_offset

Written by the main session from the structured result of workflow wf_4047bf71-f6a, because the harness does not let subagents write report files. Status: done.

## Findings

Computed on Hopper (job 1243029, Intel) and archived: riesz_program/code/ggc_offset.py, output riesz_program/data/ggc_offset.json. Rerunning ggc_coef.py reproduces ggc_coef.json bit for bit. Without quadrature error (Arb v_1, Gauss Legendre; 24 and 12 nodes agree to 2.6e-23), the Weyl tail error is 2.985e-6 to 2.989e-6 times (2k-1)^2, so that clause stands. The quadrature offset of ggc_coef.py for k=1..6 is -2.980e-5, -3.001e-5, -3.040e-5, -3.099e-5, -3.176e-5, -3.271e-5. The part that is the same for every k, -2.98e-5 = h(2-log 2pi)b_1/pi, comes from the trapezoid rule (step 0.025) starting one step from the log singularity at theta=0. The rest comes mostly from the step 0.1 beyond 60. The error range and the k=6 sentence stand; only the quadrature sentence on line 164 needs replacing, as given. The harness refused REPORT.md, so its content is in these fields.

## Files added

* /Users/vsokolov/Dropbox/papers/computing_rh_book/riesz_program/code/ggc_offset.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/riesz_program/data/ggc_offset.json
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/part6_ggc_offset/slurm-1243029.log
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/part6_ggc_offset/job-1243029.slurm

## Hopper jobs

* 1242479 (hop047, intel): environment probe, completed
* 1242819 (hop048, intel): first run of ggc_offset.py; computation finished and wrote the JSON, then the final summary print failed on an mpf format; superseded
* 1243029 (hop051, intel): archived run of ggc_offset.py, 16 cores, 15 s, exit 0; output riesz_program/data/ggc_offset.json, log misc/repro/reconstruct/part6_ggc_offset/slurm-1243029.log

## Changes the book needs

* ch/riesz.tex line 164 (paragraph before Table ri:tab:ggccoef). OLD: They consist of the error of the Weyl tail, about $3.0\times10^{-6}(2k-1)^2$, and a quadrature error of about $-3\times10^{-5}$ on $[0,200]$, the same for every $k$. NEW: They consist of the error of the Weyl tail, about $3.0\times10^{-6}(2k-1)^2$, and the error of the trapezoid rule on $[0,200]$, between $-2.98\times10^{-5}$ at $k=1$ and $-3.27\times10^{-5}$ at $k=6$ (\texttt{ggc\_offset.py}, output \texttt{ggc\_offset.json}). Of the latter, $-2.98\times10^{-5}$ is the same for every $k$: the rule of step $h=0.025$ starts one step from the logarithmic singularity of the integrand at $\theta=0$, and by Stirling's formula this makes the integral too large by $h(2-\log2\pi)\,b_1/\pi=2.98\times10^{-5}$.
* Optional, ch/riesz.tex line 164 (the present sentence is correct, since without quadrature error the k=6 value is 1.0001150). OLD: At $k=6$ the truncated value exceeds one, which $1/\zeta(12)$ cannot, and the excess is the error of the Weyl tail. NEW: At $k=6$ the truncated value exceeds one, which $1/\zeta(12)$ cannot: the error of the Weyl tail, $3.6\times10^{-4}$, exceeds $1-1/\zeta(12)=2.5\times10^{-4}$.
* No change to Table ri:tab:ggccoef. Its Thorin form and GGC factor columns are what ggc_coef.py prints, and the Hopper rerun reproduces them exactly. The range 3e-6 to 3.3e-4 (largest at k=6) and the coefficient 3.0e-6 (2k-1)^2 stand. If Part VI lists its scripts and outputs, add ggc_offset.py and ggc_offset.json (riesz_program).

## For the authors

* REPORT.md was not written. The harness refused /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/part6_ggc_offset/REPORT.md because subagents may not write report files. Its content is in findings, chapter_changes and these flags. The evidence is riesz_program/data/ggc_offset.json, slurm-1243029.log and job-1243029.slurm in misc/repro/reconstruct/part6_ggc_offset/.
* What was missing: the -3e-5 quadrature offset in the paragraph on line 164 came from the reproduction in check M3 (misc/reviews/riesz/ggc-math.md), and no archived script or output supported it.
* Per row k=1..6, relative errors (value/reference minus 1). Table vs 1/zeta(2k): -2.6811e-5, -3.1179e-6, +4.4289e-5, +1.1537e-4, +2.1010e-4, +3.2843e-4. Weyl tail (truncated integral without quadrature error vs exact): +2.9886e-6, +2.6895e-5, +7.4696e-5, +1.4637e-4, +2.4187e-4, +3.6116e-4; divided by (2k-1)^2 these are 2.9886, 2.9883, 2.9878, 2.9871, 2.9861, 2.9848 times 1e-6. Quadrature (ggc_coef.json vs the truncated integral without quadrature error): -2.9799e-5, -3.0012e-5, -3.0404e-5, -3.0988e-5, -3.1759e-5, -3.2714e-5. The product (1+Weyl)(1+quad) equals 1+table to 2e-16. The GGC factor column has the same errors.
* The words 'the same for every k' in the chapter and in M3 are only approximate: the offset spreads by 10 percent. The part that does not depend on k is 2.97965e-5 in the exponent, h(2 - log 2pi) b_1/pi with h = 0.025 (Stirling). It matches the trapezoid error on [0.025, 60] to within 2.2e-8 (k=1) to 3.0e-7 (k=6). The rest, -3.4e-9 (k=1) to -2.92e-6 (k=6), comes from the step 0.1 on [60, 200], which adds +2.7e-8 to +3.22e-6 to the exponent. The constant density on [0, 0.025] adds about -2e-9 and the authors' Weyl tail quadrature below 1e-14.
* Values without quadrature error, in case the authors would rather print them (not recommended, since the table would then no longer be ggc_coef.py output). Thorin form: 0.6079289, 0.9239633, 0.9830260, 0.9960850, 0.9992480, 1.0001150. GGC factor: 0.9772440, 0.8144496, 0.5701841, 0.3394447, 0.1742778, 0.0782578, equal to the M3 checker's reproduction to all 7 digits.
* Checks. (1) The rerun of ggc_coef.py matches ggc_coef.json with relative difference 0.0 for both columns at every k, and v_1/pi from mp.diff and from Arb agree exactly in double precision on the authors' grid. (2) Gauss Legendre with 24 and 12 nodes per panel agree to 2.6e-23 on [0,200], the largest Arb radius of v_1 is 2.6e-37, and the two rules agree to 2.3e-27 on [1900, 2000]. (3) The closed form Weyl tail agrees with mpmath quad to 1e-17. (4) Integrating rho on to 500, 1000 and 2000 changes the tail error per unit w from 2.99e-6 to 9.72e-7, -1.511e-8 and -1.606e-8. The first order prediction -arg zeta(1+iL)/(pi L^2) + (11/12)/(3 pi L^3) gives 9.73e-7, -1.519e-8 and -1.611e-8, within 1.2e-9 (5e-11 at 2000).
* Optional explanation of the Weyl tail coefficient: to first order the error is -(2k-1)^2 arg zeta(1+200i)/(200^2 pi) + (11/12)(2k-1)^2/(3 pi 200^3) = 3.08e-6 (2k-1)^2, with arg zeta(1+200i) = -0.38522. This is within 3 percent of the computed 2.99e-6 (2k-1)^2, and is the basepoint one analogue of the weighted average of S in sec:ch18:kent.
