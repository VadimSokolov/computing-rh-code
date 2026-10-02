# Reproducibility report: ttest2

Written by the main session from the structured result of workflow wf_36eed84b-56b, because the harness does not let subagents write report files. Status: partly_reproduced.

## Findings

Both scripts failed only because misc/repro/run_all.py:30 passes no arguments (ttest.py:8 ValueError, ttest2.py:7 IndexError). The archived ttest2.py run as `0.10 200 0.004 4 200 500 1000` and `0.02 220 0.001 4.2 1000 2000 5000` reproduces 17 of 18 computed entries of Table tab:tilt: all node counts, the last column, and five of six errors to the printed two digits (worst 1.644e-45 against 1.6e-45). A 116 task sbatch search finds h=0.004 the only step fitting the δ=0.1 rows; for δ=0.02, h=0.001 is the only round fit (h=0.00093 also gives 2.9e-25). The entry <5e-45 (δ=0.02, θ=1000) cannot be reproduced: every usable h prints 5.14e-45. The three entries near 1e-45 measure the 45 digit rounding inside ttest2.py, not the reference; true errors are 1e-55 to 3e-51. Times: 6.8 and 30 ms against 7 and 28.

## Scripts

* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/ttest2/run_archived.slurm
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/ttest2/search_task.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/ttest2/search.slurm
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/ttest2/collect_search.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/ttest2/compare_table.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/ttest2/post.slurm
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/ttest2/floor_check.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/ttest2/horner100.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/ttest2/diag.slurm
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/ttest2/extra_check.py

## Hopper jobs

* 1241681: archived ttest2.py and ttest.py run unchanged with the recovered command lines (hop071)
* 1241682: sbatch array 0 to 115, the search over the step h (114 distinct steps, archived ttest2.py unchanged)
* 1242113: rerun of search tasks 47 and 115 with other cuts (3.9965, 4.1991), same output as cuts 4 and 4.2
* 1241683: floor_check.py and horner100.py (hop073)
* 1241684: first collection of the search
* 1242114: final collection and table comparison
* 1242357: extra_check.py through misc/tools/hopper_run.sh (hop061)
* 1205029: earlier run of misc/repro/ttest2.slurm (hop067), read only; same digits

## Changes the book needs

* ch/ch04.tex:132 (Table tab:tilt, δ=0.02, θ=1000): `$<5\times10^{-45}$` is not what the archived code/ttest2.py prints; it prints 5.14e-45 (5.13636e-45) for every usable step from 0.00083 to 0.00116. Write `$5.1\times10^{-45}$`, or `$<6\times10^{-45}$`.
* ch/ch04.tex:121 ('The reference ... at 50 digits, which limits the comparison at heights 1000 and 2000') and the caption at ch/ch04.tex:137 ('< marks agreement to the accuracy of the reference'): the limit is the rounding of the tilted value to 45 significant digits in tomp of code/ttest2.py (str(45)). The 50 digit reference is good to 1.0e-50 or better, and the tilted values are right to 1.4e-55 (δ=0.1, θ=200), 2.1e-53 (δ=0.02, θ=1000) and 3.4e-51 (δ=0.02, θ=2000). Say so, and mark the θ=200 row at ch/ch04.tex:129 (1.6e-45, also at this floor) with '<' like rows 132 and 133. Alternative: change str(45 to str(60 in tomp; the three entries then become 2.4e-51, 9.9e-51 and 8.8e-51, the sentence at line 121 becomes true as written, and the other three entries do not change.
* ch/ch04.tex:127: the header says 'error of Re(ξ'/ξ)' but code/ttest2.py prints the modulus of the complex difference, which bounds it. Write 'error of ξ'/ξ', or give the real part errors of rows 129 to 134: 2.8e-47, 3.4e-42, 7.0e-21, 4.6e-45, 3.3e-46, 1.4e-25 (first, fourth and fifth at the 45 digit floor). The relative error of Re(ξ'/ξ) at θ=5000 is 1.24e-25, so '25 digits at height 5000' (line 121) holds.
* ch/ch04.tex:173: 'gives balls of radius $10^{-52}$ at the same point' (θ=100) is not reproduced: at θ=100 the direct evaluator gives radius 10^-55.2 for Re(ξ'/ξ) and relative radius 10^-55.5 for ξ (horner100.py, job 1241683); 2.6e-52 is the radius ttest2.py prints at θ=200. Write 'about $10^{-55}$'. The other numbers of that paragraph are reproduced with h=0.004, cut 4: growth 1.3107 per step (median over the last 1500 steps), Horner radius 10^109.1 for P(E), midpoint right to 10^-57.3 relative, coefficient 2.0e-12 at x=2, 1.3107^1500 = 10^176.3.

## For the authors

* Record the command lines in a comment at the top of code/ttest2.py and in ch/appA.tex:11 (which now says only that the runs take δ, the precision, h and the cut from the command line): `python ttest2.py 0.10 200 0.004 4 200 500 1000` (2001 nodes) and `python ttest2.py 0.02 220 0.001 4.2 1000 2000 5000` (8401 nodes). Optionally give h=0.004 with |x|≤4 and h=0.001 with |x|≤4.2 in the caption of tab:tilt, as round 1 (U03 B7d, C18) asked.
* Confirm h=0.001 (cut 4.2) for the δ=0.02 run: the table does not determine it, since h=0.00093 (cut about 3.906) prints 2.94e-25, which also rounds to 2.9e-25; h=0.001 is the only round value that fits, and it reproduces 2.89e-25. The δ=0.1 run is determined: among 47 steps only h=0.004 (any cut in (3.996, 4]) gives both 9.8e-42 and 1.8e-20; its neighbours 0.003999 and 0.004001 give 4.0e-42 and 1.5e-42 at θ=500.
* code/ttest.py: record `python ttest.py 0.1 200 0.004 4`. It prints an error of 2.02e-15 at θ=14 because it evaluates at 1/2+double(0.05) and takes the reference at double(0.55), 4.16e-17 apart (with the exact inputs of ttest2.py its Horner value agrees to 5.8e-40), and nan at θ≥200, where the Horner ball of ξ contains zero. It has no θ=100, so the Horner numbers of Section sec:ch4:balls come from an unarchived run; add 100 to its heights and use acb(arb('0.05'), th) and mp.mpf('0.55').
* misc/repro/run_all.py:30 starts ttest.py and ttest2.py with no arguments, the cause of rc 1; pass them the command lines. misc/reviews/round1/U03-verify.md:41 says the ttest2.slurm run failed; its job 1205029 completed and printed exactly the numbers reproduced here.
* Timings are machine dependent: medians of 20 evaluations are 6.8 to 6.9 ms per point with 2001 nodes (book 7) and 29.8 to 30.3 ms with 8401 nodes (book 28) on an Intel Xeon Gold 6240R at 2.40 GHz; the book does not name its machine. Deterministic outputs agree digit for digit across jobs 1205029, 1241681 and the search tasks on other nodes.
* REPORT.md was not written: the harness refused a report file from this subagent, and it was not written by another route. Its content is in these fields; the numbers are in /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/ttest2/out/ (post_1242114.log for the search and the comparison, compare.json, floor_check.json, horner100.json, extra_check.json, search_results.json, search_raw.json, jobs.txt).
