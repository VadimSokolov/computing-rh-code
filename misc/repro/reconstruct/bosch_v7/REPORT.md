# Reproducibility report: bosch_v7

Written by the main session from the structured result of workflow wf_36eed84b-56b, because the harness does not let subagents write report files. Status: reproduced.

## Findings

Cause: the archived file joins 13 scripts; block 13 execs gaptest.py (lines 390, 392), which is not archived, though block 12 is its text. Fix: comment out line 390, make line 392 out=[]. The fixed copy runs rc 0 in 403 s; its 95 lines match the archive byte for byte except one timing (627 numbers, max rel diff 0), and match the unchanged archive run beside a reconstructed gaptest.py. Reconstructions reproduce output lines 37 to 39 and 109 exactly. In stale lines 93 to 102 the printed code differs on 4 cases, giving the values of lines 103 to 106. At t=0.003, 150 digit series and Talbot both give 3.32007507753e-72 (max term 2.694e20, 91.9 digits cancel). All 146 book checks pass; no corrections needed. The harness blocked REPORT.md; the report is in the final message.

## Scripts

* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/bosch_v7/verify_bosch_xi_v7_fixed.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/bosch_v7/gaptest.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/bosch_v7/supp_lines37_39.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/bosch_v7/supp_line109.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/bosch_v7/supp_gaps_table.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/bosch_v7/residue_hp_task.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/bosch_v7/residue_hp_collect.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/bosch_v7/talbot_hp.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/bosch_v7/compare_bosch_v7.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/bosch_v7/main.slurm
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/bosch_v7/residue.slurm
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/bosch_v7/talbot.slurm
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/bosch_v7/collect.slurm
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/bosch_v7/compare.slurm

## Hopper jobs

* 1242945 (array tasks 0 to 4: fixed copy rc 0 in 403 s; unchanged archive plus reconstructed gaptest.py rc 0 in 401 s; supp_lines37_39; supp_line109; supp_gaps_table)
* 1242946 (array tasks 0 to 49: zeros and residue coefficients at 150 digits, n up to 250)
* 1242947 (array tasks 0 and 1: Talbot inversion at 150 digits, t=0.003 and t=0.005)
* 1243307 (collector residue_hp_collect.py; its comparison step stopped on a regex bug in my script)
* 1243500 (compare_bosch_v7.py, rc 0, 146 of 146 checks pass)

## For the authors

* code_bosch/verify_bosch_xi_v7.py and bosch_hcm_xi/verify_bosch_xi_v7.py: delete line 390 and replace line 392 by out=[] (block 12 already defines P and make_g); the script then runs to the end in about 400 s
* Add the code for output lines 37 to 39 and 109 (supp_lines37_39.py and supp_line109.py reproduce them exactly; the 148 point grid of line 39 is a guess, say which grid was used)
* Extend the case list at script line 414 to the eleven computed rows of Table bs:tab:gaps and drop the stale output lines 93 to 102 (four of them, [2,3], [2,3,5,6], [3,4], [4,5], are superseded)
* In block 6, add t=0.003 at 150 digits with at least 220 zeros (or Talbot at 150 digits) so the output contains the 3.32007507753e-72 that boschB line 341 and appBosch line 28 quote; lines 56 and 60 at t=0.003 are rounding noise
* Record the join: run in an empty folder, then cat the nine block files in block order into verify_bosch_xi_v7_output.txt; stdout is not the printed output
* After the script is fixed and rerun, rewrite appBosch line 28 and the none rows and caption of Table tab:app:bosch (line 24), which describe the present archive
* misc/repro: run the fixed script in run_all.py and have compare.py compare the joined block files (not logs/bosch_v7.out) with the archive, treating timings apart; the zeros dependency is unnecessary
* Minor: the label of output line 66 says [1e-6,100] but the code checks 10^(-5.9) to 10^2; boschB line 265 states it correctly
* boschA numbers not made by this script (integral of g 0.99996, 390 point scan, Table bh:tab:cm, sech product, K_it(2) zero) have no archived script and were not checked
