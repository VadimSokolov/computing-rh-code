# Reproducibility report: riesz_explicit

Written by the main session from the structured result of workflow wf_36eed84b-56b, because the harness does not let subagents write report files. Status: reproduced.

## Findings

All four fields reproduced bit for bit. Definitions: maxdev = max |R(x)/x^(1/4) − E(x)| over the 36 salem_riesz.json grid points in [1e3,1e10], E being the explicit formula (60 zeros, 3 trivial terms) written as the pred lambda of figs_riesz.py and check_large.py at mpmath 20 digits; signal = max |R(x)/x^(1/4)| over the same points (Möbius values, attained at x = 1e3); amp1 = 2|Γ(1−ρ1/2)/(2ζ'(ρ1))| with γ1 from g_1_400.npy; triv = k!/(2ζ'(−2k)). The reconstructed riesz_explicit.py (job 1241435) writes a file byte identical to both archived copies, which are identical: 9 of 9 numbers bit identical. maxdev sits at x = 1e10 and is the truncation error of the Möbius sum at N = 1e7 (job 1241800), not an explicit formula error; ch/strip.tex line 302 wrongly blames double precision cancellation. REPORT.md was refused by the harness (subagents may not write report files); details are in author_flags.

## Scripts

* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/riesz_explicit/riesz_explicit.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/riesz_explicit/explore.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/riesz_explicit/compare_riesz_explicit.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/riesz_explicit/truncation_check.py

## Hopper jobs

* 1241214
* 1241435
* 1241595
* 1241800

## Changes the book needs

* ch/strip.tex line 302 (Section wr:sec:riesz): the clause 'beyond $x\approx10^{10}$ double precision cancellation limits the accuracy' gives the wrong cause. Suggested text: 'beyond $x\approx10^{10}$ the truncation at $N=10^7$ limits the accuracy: the omitted terms are about $-x^2\sum_{n>N}\mu(n)n^{-4}\approx9\times10^{-26}x^2$, that is $3\times10^{-8}$ in $R(x)/x^{1/4}$ at $x=10^{10}$ and $9\times10^{-5}$ at $x=10^{12}$'. Evidence (truncation_check.py, job 1241800): sum over 1e7 < n <= 1e8 of mu(n)/n^4 = -9.11e-26 predicts -2.88e-8, -1.62e-6, -9.11e-5 in R(x)/x^(1/4) at 1e10, 1e11, 1e12, and the observed deviations of salem_riesz.json from the explicit formula are -3.11e-8, -1.59e-6, -9.07e-5; with N = 1e8 in the same double precision the largest deviation falls to 6.5e-9 on [1e3,1e10] and 1.5e-7 on [1e3,1e12]. No number quoted in the book changes: 3.1e-8, 3e-3 (2.96e-3), 7.8e-5 (7.775e-5), -16.4 and -16.42, 125.25, -508.50 (-508.495) and -1.6e-4 (-1.642e-4 at x = 1e4) all agree with the file.
* ch/strip.tex line 302 (optional clarification): after 'reproduces the computed values to $3.1\times10^{-8}$ on $[10^3,10^{10}]$' say that the largest difference occurs at $x=10^{10}$ and is this truncation error, while below $10^9$ the agreement is $6\times10^{-10}$ (per point deviations in rerun.log).

## For the authors

* Add the missing script. basepoint_one/wiener/riesz_explicit.json had no script (README.md line 5: 'script inline in the paper's description'). misc/repro/reconstruct/riesz_explicit/riesz_explicit.py reproduces it byte for byte (182 bytes). The script is marked at its top as a reconstruction and not the authors' script. It uses numpy, mpmath at 20 digits and the coef, triv and pred lines of riesz_program/code/figs_riesz.py and check_large.py, reads salem_riesz.json and g_1_400.npy, and writes riesz_explicit.json. Please add it, or the original if it survives, as basepoint_one/wiener/riesz_explicit.py and riesz_program/code/riesz_explicit.py. Record the command lines 'python3 salem_riesz.py; python3 riesz_explicit.py' (in riesz_program, run them from code/ after ./setup.sh). Replace the phrase on README.md line 5 with the script name, and list the script in the code section of riesz_program/README.md.
* Exact definitions (explore.py, job 1241214, explore1.log). maxdev uses the ratio form |R/x^(1/4) - pred| of figs_riesz.py. The other algebraic forms give 3.1103171913005e-08 or 3.11031719130027e-08, within 1.1e-13 relative: taking the difference before dividing by x^(1/4), or the complex form of twisted.py summed over both signs of gamma. signal comes from the Möbius values; max |E| would give 0.0029587439723876663 (relative 8.7e-8). amp1 used gamma_1 = g_1_400.npy[0] = 14.134725141734695; with mp.zetazero(1) it would be 7.775062764452562e-05 (relative 7e-16). triv comes out as the same doubles at every precision from 15 to 50 digits. maxdev is attained at x = 1e10 and signal at x = 1e3.
* Comparison (compare_riesz_explicit.py, job 1241595, compare.log). The rerun (riesz_explicit_rerun.json, job 1241435) is byte identical to basepoint_one/wiener/riesz_explicit.json and to riesz_program/data/riesz_explicit.json, and those two are byte identical to each other. So are the two copies of salem_riesz.json, salem_riesz.py and g_1_400.npy. 9 of 9 numbers are bit identical and the worst relative difference is 0. The computation is deterministic, so bit identity is the agreement to expect.
* What maxdev measures (truncation_check.py, job 1241800, truncation_check.log/json). The deviation is 2.56e-10 at x = 1e3, which is the omitted fourth trivial term 4!/(2 zeta'(-8)) = 1442.97 (it gives 2.57e-10 there). It stays below 6.3e-10 on [1e4,1e9], then grows like x^(7/4) to 3.11e-8 at 1e10. The N = 1e7 recomputation reproduces all 56 values of salem_riesz.json bit for bit. The Mertens values are M(1e7) = 1037 and M(1e8) = 1928. riesz_large.py (check_large.json, K = 200) agrees with the explicit formula to 1.9e-9 at 1e10. So maxdev is the truncation error of salem_riesz.py at the single end point x = 1e10, not an error of the explicit formula.
* ch/strip.tex line 276 (caption of Figure wr:fig:salem) names salem_riesz.py for salem_riesz.json but gives no script for riesz_explicit.json. Once the script is added, write 'and in riesz_explicit.json, written by riesz_explicit.py'.
* The reproduction drivers misc/repro/run6.py and misc/repro/run_riesz.py should run riesz_explicit.py after salem_riesz.py. At present run_riesz.py line 22 copies the archived JSON in as an input of figs_riesz.py.
* The authors' source notes give the same wrong cause for the accuracy limit (double precision cancellation): riesz_program/notes/riesz_thorin.tex line 65 and basepoint_one/merged/basepoint_one_merged.tex line 981.
* REPORT.md was not written because the harness refuses report .md files from subagents. The evidence is in misc/repro/reconstruct/riesz_explicit/: explore1.log, rerun.log, compare.log, truncation_check.log, truncation_check.json and riesz_explicit_rerun.json. The Hopper run folders are under /scratch/vsokolov/rh_book_repro/agents/reconstruct_riesz_explicit/.
