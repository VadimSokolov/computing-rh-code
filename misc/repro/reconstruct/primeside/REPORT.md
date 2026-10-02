# Reproducibility report: primeside

Written by the main session from the structured result of workflow wf_36eed84b-56b, because the harness does not let subagents write report files. Status: partly_reproduced.

## Findings

primeside.py runs unchanged with a shim polya.py, a wrapper of polya_flint setup(420, 200, 3). The regenerated results_prime.json has 162 numbers within 10⁻⁹, 148 bit for bit; the Pólya column differs by at most 4 ulp, and none of 196 tested conventions reproduces its last bits. plots.py on it regenerates extra.json and sens2.json byte for byte, and all 58 printed numbers of tab:prime and tab:split agree. No script writes or reads results_main.json, and it does not depend on primeside. A reconstruction reproduces 117 of its 135 numbers (not tv and sup, the tail model floors, nor the timings). The earlier "agrees" for results_prime.json compared the file with its own copy. REPORT.md was not written because the harness refuses report files from subagents; its text is in the final message.

## Scripts

* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/primeside/polya.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/primeside/shim_variants.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/primeside/shim_variants2.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/primeside/compare_primeside.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/primeside/archimedean_check.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/primeside/book_tables_check.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/primeside/main_compute.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/primeside/main_analyse.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/primeside/main_tail_probe.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/primeside/main_tail_probe2.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/primeside/pipeline.slurm
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/primeside/variants2.slurm
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/primeside/main_array.slurm
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/primeside/main_analyse.slurm
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/primeside/probe.slurm
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/primeside/probe2.slurm

## Hopper jobs

* 1241557
* 1242889
* 1241999
* 1242243 (array tasks 0 to 5)
* 1242244
* 1242481
* 1242621

## For the authors

* code/primeside.py:3 imports a module polya that is not in the archive. Add code/polya.py (the shim) or use polya_flint setup(420, 200, 3) at lines 3, 28 and 29 as plots.py:102 to 107 do; the Pólya column then moves by at most 5 ulp and no printed number changes.
* data/results_main.json has no writer or reader in the archive. Add its script or remove the file and take γ1700 = 2197.2587 (ch07.tex:25, ch18.tex:91) from the zero list, for instance stored in analysis.json. The tail model behind its tv and sup cannot be recovered.
* Record the inputs of plots.py in appA or the README: plotcurves.npz, analysis.json, hcm.json, results_prime.json and an existing fig/ folder; no arguments.
* Optional: Table tab:prime row α = 3, θ = 50 (ch02.tex:161) with the closed form archimedean integral would read 1.120967015289, 1.120967015289, 1.3×10⁻¹³; the current text at ch02.tex:152 is correct.
* Minor: plots.py:84 at ε = 0.1 evaluates P(0, θ − γ5), which is 0/0 at θ = γ5; Figure fig:weaklimit misses one sample there.
* The compare_report.md lines saying results_prime.json, extra.json and sens2.json agree came from job 1218076, which used a copy of the archived results_prime.json; this reconstruction supersedes them.
