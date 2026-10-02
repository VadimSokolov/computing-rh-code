# Reproducibility report: kernel_fig

Written by the main session from the structured result of workflow wf_65af78f3-bcf, because the harness does not let subagents write report files. Status: done.

## Findings

fig/kernel.pdf replaced (old kept as old_kernel.pdf). kernel_eta.py copies the kernel block of code/plots3.py, changing only the legend (η=0.1000, 0.3927, 0.6854, 0.7654) and the right title (|Φ(x+iη)| inside the strip). Run on Hopper with matplotlib 3.10.8, the archived PDF's version, the unchanged code reproduces that PDF exactly (only /CreationDate differs). New against archived: 61 of 66 strokes, 30 of 31 fills, 29 of 34 text blocks identical to 0.0 pt; pixel changes lie only in the title and legend boxes; by eye identical apart from the labels. Curve values match a fresh computation to 4.7e−8 in log10. No number changes. Caption sentence to delete (ch02:47): 'in the figure these distances are labelled $\delta$'. The harness refused REPORT.md; its content is in these fields.

## Files added

* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kernel_fig/kernel_eta.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kernel_fig/mpl3108.slurm
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kernel_fig/pdf_compare.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kernel_fig/vertex_check.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kernel_fig/peak_check.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kernel_fig/old_kernel.pdf
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kernel_fig/old_kernel.png
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kernel_fig/kernel.pdf
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kernel_fig/kernel.png
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kernel_fig/kernel_curves.json
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kernel_fig/kernel_oldlabels.pdf
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kernel_fig/kernel_oldlabels.png
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kernel_fig/kernel_oldlabels_curves.json
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kernel_fig/pdf_compare.json
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kernel_fig/diff_old_new.png
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kernel_fig/diff_old_control.png
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kernel_fig/diff_control_new.png
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kernel_fig/run_compare.log
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kernel_fig/run_peak.log
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kernel_fig/slurm_1243711_mpl3108.log
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kernel_fig/mpl3112/ (the same outputs, comparison and logs from the matplotlib 3.11.2 runs: kernel.pdf, kernel_oldlabels.pdf, their PNG and JSON files, pdf_compare.json, diff_*.png, run_*.log, slurm_1243258_new.log, slurm_1243259_old.log)

## Files replaced

* /Users/vsokolov/Dropbox/papers/computing_rh_book/fig/kernel.pdf

## Hopper jobs

* 1243258: kernel_eta.py with the new labels, matplotlib 3.11.2 of the shared venv (agents/p13_kernel_fig/new; outputs in mpl3112/)
* 1243259: kernel_eta.py old, the archived labels as a control, matplotlib 3.11.2 (agents/p13_kernel_fig/old)
* 1243360: first pdf_compare.py run; stopped with an error (tick marks taken for spines), fixed
* 1243523: pdf_compare.py on the 3.11.2 set with earlier file names (superseded by 1244007, same numbers)
* 1243676: vertex_check.py on the 3.11.2 set (why one curve's compared vertex count differs by one there: path simplification near the clipped ends)
* 1243711: mpl3108.slurm, kernel_eta.py in both modes with matplotlib 3.10.8 unpacked into agents/p13_kernel_fig/mpl3108_pkg; produced the replacement fig/kernel.pdf and its control
* 1243742: pdf_compare.py on the 3.10.8 set with earlier file names (superseded by 1243883, same numbers)
* 1243773 and 1243809: peak_check.py, side check of the ch04 numbers about these curves
* 1243883: pdf_compare.py final on the 3.10.8 set (pdf_compare.json, diff_*.png)
* 1244007: pdf_compare.py final on the 3.11.2 set (mpl3112/pdf_compare.json)

## Changes the book needs

* ch/ch02.tex line 47, caption of fig:ch2:kernel. Old: '(the tilts of Section~\ref{sec:newcomp}); in the figure these distances are labelled $\delta$. As the line approaches' New: '(the tilts of Section~\ref{sec:newcomp}). As the line approaches'. The clause 'in the figure these distances are labelled $\delta$' goes, and the semicolon before it becomes a full stop. The rest of the caption stays as it is.
* No other chapter changes. The figure has no data file and no number changes. A grep of ch/*.tex for 0.6854, 0.3927, 0.7654, kernel.pdf and fig:ch2:kernel finds only ch02:42, 46 to 48, ch04:98 and 99, and appA:32. The ch04 lines are rows of Table tab:ch4:Meta and stay unchanged: their maxima 0.981 and 3.24 match the blue and red peaks, 0.980498 and 3.24183. appA:32 names plots3.py as the script of fig:ch2:kernel, which stays true once code/plots3.py gets the two line change flagged below.

## For the authors

* code/plots3.py lines 20 and 21 (not edited; code/ belongs to the authors) should get the same change, so that the archive redraws the new figure. Line 20 old: ` ax[1].semilogy(xx,np.maximum(v,1e-30),lw=1.4,label=r'$\delta=%s$'%d)` new: ` ax[1].semilogy(xx,np.maximum(v,1e-30),lw=1.4,label=r'$\eta=%.4f$'%(np.pi/4-d))`. Line 21: replace `ax[1].set_title(r'$|\Phi(x+i(\pi/4-\delta))|$ inside the strip')` with `ax[1].set_title(r'$|\Phi(x+i\eta)|$ inside the strip')` and leave the rest of the line as it is. With matplotlib 3.10.8 the edited block then redraws fig/kernel.pdf exactly.
* Matplotlib version: the shared Hopper venv has 3.11.2. With the same code it places the axes 2.16 pt higher at the top and 0.16 pt lower at the bottom, so curve vertices move by up to 4.4 pt (the values are unchanged). The replacement was therefore drawn with 3.10.8, the archived version, unpacked into /scratch/vsokolov/rh_book_repro/agents/p13_kernel_fig/mpl3108_pkg; the shared venv was not touched. The 3.11.2 outputs are kept in mpl3112/.
* Side check outside the item, with no change needed. The grid maxima of Table tab:ch4:Meta (step 0.002, |x| at most 6) recompute as 0.980498, 3.24183, 76.497, 3938.31 and 22207.2 (the table has 0.981, 3.24, 76.5, 3938 and 22207). At δ=0.02, |Φ(iη)| is 0.0121869 (ch04:90 says 0.012). The maximum itself lies at x=(1/2)log 2 and is 3942.36 at δ=0.02 and 22298.6 at δ=0.01. So 'a maximum of 3938' in ch04:90 is the grid value; the authors could write 3942 there if they prefer.
* The R2-02 flags on fig/kernel.pdf (the last list of R2-02-changes.md and item 5 of R2-02-verify.md) are resolved by this item. The logs themselves were not edited.
* REPORT.md was not written because the subagent harness refuses report files. The report content (what was wrong, the jobs, the comparison and the exact changes) is in this structured output, and the evidence files are in the p13_kernel_fig folder.
