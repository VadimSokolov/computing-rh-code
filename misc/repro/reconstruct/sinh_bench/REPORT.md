# Reproducibility report: sinh_bench

Written by the main session from the structured result of workflow wf_36eed84b-56b, because the harness does not let subagents write report files. Status: reproduced.

## Findings

data/sinh_bench.json holds the sinh clock benchmark of Figure fig:ch3:bench (code/plots3.py lines 117 to 123): cell distance and mass on 0<θ≤50.5 at five ε, as arctangent sums over the lattice cut at |k|≤400000. The reconstructed sinh_bench.py (numpy, cumulative mass over atom and reflection pairs) reproduces all 15 numbers within 1e-9 (worst 1.12e-12, celltv at ε=0.01; 11 of 15 bit for bit). An mpmath check with log Gamma tails shows the archive equals the exact truncated sums to 2.1e-12 and pins K=400000. Against the closed forms the truncation raises celltv by 6e-5 to 8e-5 relative and lowers the mass defect by 1.3 percent, which cannot be seen in the figure. Every value printed in ch03 and ch15 already uses the closed forms, so the text needs no correction. The harness refused REPORT.md; its content is in author_flags.

## Scripts

* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/sinh_bench/sinh_bench.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/sinh_bench/check_sinh_bench.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/sinh_bench/probe_arctan.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/sinh_bench/probe_codings.py

## Hopper jobs

* 1242050: sinh_bench.py (K=400000), writes sinh_bench.json; hop062, folder /scratch/vsokolov/rh_book_repro/agents/reconstruct_sinh_bench_run2
* 1242115: sinh_bench.py closed, writes sinh_bench_closed.json; folder reconstruct_sinh_bench_closed
* 1242241: check_sinh_bench.py on the archived file (book_data copy, same sha256 bc79a1a4...) and both outputs, writes check_sinh_bench.json; folder reconstruct_sinh_bench_check2
* 1242522: probe_arctan.py (numpy X86_V4 arctan dispatch on and off), writes probe_arctan.json; folder reconstruct_sinh_bench_probe
* 1242738: probe_codings.py (twelve codings of the same sums), writes probe_codings.json; folder reconstruct_sinh_bench_codings
* superseded: 1241522 and 1241623 (first coding with one arctan2 per atom and cell, agreement 1.6e-12), 1242373 and 1242431 (probe_arctan runs that failed or did not switch the dispatch)

## For the authors

* Add a script for data/sinh_bench.json to code/, for example misc/repro/reconstruct/sinh_bench/sinh_bench.py with its reconstruction header replaced, and record its command line: python3 sinh_bench.py (K=400000, about 2 s) reproduces the archived file to 1.1e-12; python3 sinh_bench.py closed writes the closed forms of Section sec:ch3:bench.
* Decide which data the figure uses. The archived file is the lattice truncated at |k|<=400000. The closed forms, which the text quotes, differ from it by at most 7.9e-5 relative, and the difference cannot be seen in fig/bench.pdf. To switch, as U02 F4 and R2-02 recommend, use the closed form values in sinh_bench_closed.json: celltv 0.4968485199452337, 0.1250356005955584, 0.06313587666093763, 0.012667666542139005, 0.006334469699604266 and mass 50.00315148005477, 50.000630315782445, 50.000315158200166, 50.000063031659806, 50.000031515830216 at eps 0.5, 0.1, 0.05, 0.01, 0.005. Write them under the name sinh_bench.json, which plots3.py line 118 reads. To keep the truncated file, state the cut |k|<=400000 in the script or in Appendix A.
* ch/appA.tex line 23: remove sinh_bench.json from the data files without a script (Two data files, hankel.json and spacing.json, have no script that writes them). Line 33 (Table tab:app:scripts, row ch:ggc): list sinh_bench.py (sinh_bench.json) beside plots3.py.
* Regenerate fig/bench.pdf without the title 'Irregular spacing costs a factor of three' (code/plots3.py line 122), as U02 F4 and R2-02 already ask. The last clause of the caption at ch/ch03.tex line 145 can then go.
* Comparison detail. All 15 archived numbers agree at relative tolerance 1e-9. Equal bit for bit: the five eps, and celltv and mass at eps 0.5, 0.05 and 0.005. The rest differ as follows: celltv at eps 0.1 by 5.7e-14, celltv at eps 0.01 by 1.12e-12, and mass at eps 0.1 and 0.01 by one and two ulp. Archived against the exact truncated sums (mpmath at 40 digits with log Gamma tails; the formula matches direct sums at K=300 to 3.4e-41): at most 2.1e-12 apart. The K implied by each of the ten numbers lies between 399999.986 and 400000.008, and one unit of K changes celltv at eps 0.005 by 1.55e-10. The four numbers that differ in their last bits do not change with numpy's AVX512 arctan switched off, and none of 23 codings reproduces all ten. This is rounding in double precision, so a rerun on other hardware should agree to about 1e-12.
* All 13 benchmark values the book prints are the correctly rounded closed forms: ch03 lines 134, 140, 149 and 159, and ch15 lines 3 and 37. The archived file would instead give 3.111e-3, 0.4969, 1.26697 and 0.00622 in place of 3.151e-3, 0.4968, 1.26694 and 0.00630.
* The harness refused REPORT.md because subagents may not write report files. The number by number record is in misc/repro/reconstruct/sinh_bench/check_sinh_bench.json, logs/slurm-1242241.log, probe_arctan.json and probe_codings.json.
