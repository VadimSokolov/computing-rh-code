# Reproducibility report: spacing

Written by the main session from the structured result of workflow wf_36eed84b-56b, because the harness does not let subagents write report files. Status: reproduced.

## Findings

data/spacing.json holds four things from ch15. leaks: the 79 first order cell errors, with cells from b_0=0 to b_79=200 and the zero sum truncated at 1700 ordinates. coef: their absolute sum, 5.415. mean, var, min: statistics of the 1699 spacings unfolded by ϑ(γ)/π; var divides by 1699. The reconstructed spacing.py computes zeros with mpmath.zetazero at 30 digits (identical to Arb in float64), ϑ with mpmath.siegeltheta, and sums in float64 numpy. Its output is byte identical to the archived file: 83 of 83 numbers match, both with fresh zeros and with the earlier rerun's zero lists. A first version agreed to 2.55×10⁻¹²; a search over float64 recipes found the exact one. With 40 bins, the same spacings also reproduce the bars of fig/spacing.pdf. Every ch15 number agrees with the data, so the book needs no correction. REPORT.md was not written because the harness refuses report files from subagents; its content is in these fields.

## Scripts

* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/spacing/spacing.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/spacing/run_spacing.slurm
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/spacing/compare_spacing.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/spacing/variants.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/spacing/variants2.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/spacing/variants3.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/spacing/pdf_bars.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/spacing/fig_check.py

## Hopper jobs

* 1241556 (hop059): first version of spacing.py, run with fresh zeros and with the zero lists of the earlier rerun, then compared with the archive: all 83 numbers within 10⁻⁹, worst 2.55×10⁻¹² (min). Also ran variants.py on unfolding, variance divisor, finite ε and boundaries
* 1242010 (hop053): variants2.py tried 198 recipes for the leaks; the best reproduced 61 of 79 leaks bit for bit, and the other 18 differed by one unit in the last place
* 1242128 (hop062): variants3.py tried orders of the final differences; one order reproduces all 79 leaks and coef bit for bit
* 1242323 (hop055): final spacing.py, run twice and compared: 83 of 83 numbers identical, output byte identical to data/spacing.json (md5 7ad57c0957ef4851617c32735be6d001)
* 1242576 (hop064): pdf_bars.py read the book's fig/spacing.pdf: 40 contiguous bars from the smallest spacing (0.0892529) to the largest (2.2025)
* 1242659 (hop047): pdf_bars.py dumped the drawing operators of fig/spacing.pdf (colours, line widths, alpha, fonts, 6 by 3.4 inches)
* 1242815 (hop060): final run_spacing.slurm: byte identical again. fig_check.py: the 40 bars match the book's figure to 6×10⁻⁹ in edges and 4×10⁻⁹ in heights (data units)

## For the authors

* Add a script that writes data/spacing.json and draws fig/spacing.pdf, for example code/spacing.py. The reconstruction /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/spacing/spacing.py reproduces the data file byte for byte and the 40 bars of the figure to PDF precision. Run it with: cd code && python3 spacing.py. It regenerates the zero lists if they are absent (53 s on 16 cores; 11 s with the lists present)
* ch/appA.tex:23: remove spacing.json from the data files that have no script. ch/appA.tex:47 (Table tab:app:scripts, row ch:spacing): add the script with spacing.json, and add Figure fig:ch15:spacing. ch/appA.tex:21: nine scripts load the zero lists, which becomes ten
* Archive the generator of the zero lists (misc/repro/zeros.py or the authors' own) and state its precision. mpmath.zetazero at 20 or 30 digits and Arb give identical float64 lists; at 15 digits one ordinate differs by one unit in the last place
* ch/ch15.tex:24 (AUTHORS comment): κ_cell from Ξ′/Ξ at the 80 cell boundaries is 5.436650959450707 (key kappa_cell_xi in misc/repro/reconstruct/spacing/out/spacing_extra.json). If the text is to quote 5.43665, have the script write it. At the same boundaries, −Ξ′(200)/(πΞ(200)) = 0.3981063
* Exact recipe of the archived numbers: x_k is float(mpmath.siegeltheta(γ_k)) divided by np.pi, and var is the population variance. Each leak is (A(b_{k−1})−A(b_k))/π + (C(b_{k−1})−C(b_k))/π, where A(c) and C(c) are numpy sums of 1/(c−γ_j) and 1/(c+γ_j) at the float64 midpoints of the float64 ordinates; coef is the numpy sum of the absolute leaks. The figure is a density histogram with 40 bins over the range of the spacings, lightsteelblue bars with white edges, and a crimson surmise and gray dashed Poisson curve of width 1.5
* The file holds first order coefficients: cell masses computed at ε = 0.005 would differ by 8×10⁻⁵. The archived leaks use float64 midpoints, not the long double cells of code/analyse.py; on Hopper those cells differ from the midpoints in 17 of 80 boundaries by one unit in the last place
