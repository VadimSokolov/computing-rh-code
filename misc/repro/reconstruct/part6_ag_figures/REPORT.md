# Reproducibility report: ag_figures

Written by the main session from the structured result of workflow wf_4047bf71-f6a, because the harness does not let subagents write report files. Status: done.

## Findings

No script drew fig/ag_almost.pdf or fig/ag_cancel.pdf. Their vector paths show that the ag_almost density panels used step 0.005 on [284, 300], finer than almost.json (step 0.5). figs_almost.py redraws both figures with v_α/π titles and α in the exponent legend. Its densities equal almost.json bit for bit at all archived points, and it reproduces every stored high.json minimum and negative mass exactly. Every curve vertex of the old and new PDFs lies within 1×10⁻⁷ of the exact data. Limits, colours and widths match, and the PNGs agree by eye. Both fig/ files were replaced. figs_almost.json stores the dip (−5.3873 at θ = 292.40187, α = 0.7) and the background statistics, so "about −5.4" and line 127 stand. Only the two captions change. REPORT.md was not written because the harness refuses report files from subagents; its content is in these fields.

## Files added

* /Users/vsokolov/Dropbox/papers/computing_rh_book/riesz_program/code/figs_almost.py (draws fig/ag_almost.pdf and fig/ag_cancel.pdf from almost.json, more.json, high.json and zeros/g_*.npy, density written v_alpha/pi)
* /Users/vsokolov/Dropbox/papers/computing_rh_book/riesz_program/data/figs_almost.json (exactness checks, dip depths near 292.40 at alpha 0.6 and 0.7, high.json background min and mean)
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/part6_ag_figures/old_ag_almost.pdf and old_ag_cancel.pdf (the replaced files)
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/part6_ag_figures/ag_almost.pdf and ag_cancel.pdf (copies of the new files)
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/part6_ag_figures/pdfpaths.py, mapdata.py, inspect_old.py, compare_figs.py, render_png.py, final_run.py (tools: read the PDF vector paths, map them to data coordinates, check every vertex against the archived data, render PNGs)
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/part6_ag_figures/compare_figs.json and old_ag_almost.pdf.data.json, old_ag_cancel.pdf.data.json, ag_almost.pdf.data.json, ag_cancel.pdf.data.json, figs_almost.json (copy)
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/part6_ag_figures/png/ (14 PNGs: both pages, each panel old above new, close views of the four dips)
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/part6_ag_figures/logs/ (slurm logs of the 7 jobs)

## Files replaced

* /Users/vsokolov/Dropbox/papers/computing_rh_book/fig/ag_almost.pdf (old copy: misc/repro/reconstruct/part6_ag_figures/old_ag_almost.pdf)
* /Users/vsokolov/Dropbox/papers/computing_rh_book/fig/ag_cancel.pdf (old copy: misc/repro/reconstruct/part6_ag_figures/old_ag_cancel.pdf)

## Hopper jobs

* 1242523 parse_old: pdfpaths.py on the old PDFs
* 1242620 map_old: old curves mapped to data coordinates
* 1242773 inspect_old: old density curves lie on a grid of step 0.005 over [284, 300]; marker size 6
* 1242907 figs_almost_run1: first run of figs_almost.py
* 1243072 compare1: first vertex comparison (its pdftoppm step failed: the Hopper poppler module has no pdftoppm, so PyMuPDF 1.28.2 went into a private pylib folder)
* 1243494 render1: PNG renderings
* 1243683 final: figs_almost.py, compare_figs.py and render_png.py with the final scripts; the files in fig/ and riesz_program/data/figs_almost.json come from this job (the JSON is byte identical to 1242907)

## Changes the book needs

* ch/almost.tex line 179 (caption of ag:fig:cancel): old text 'the density $v_{0.6}/\pi$, written $\rho_{0.6}$ in the panels, near the heights' ; new text 'the density $v_{0.6}/\pi$ near the heights'
* ch/almost.tex line 179 (same caption, end): old text 'at the centre of the window (\texttt{high.py}).}' ; new text 'at the centre of the window (\texttt{high.py}). Drawn by \texttt{figs\_almost.py}.}'
* ch/almost.tex line 186 (caption of ag:fig:almost): old text 'the density $v_\alpha/\pi$, written $\rho_\alpha$ in the panels, at $\alpha=0.6$' ; new text 'the density $v_\alpha/\pi$ at $\alpha=0.6$'
* ch/almost.tex line 186 (same caption, end): old text 'with and without the hypothetical zeros (\texttt{almost.py}).}' ; new text 'with and without the hypothetical zeros, evaluated as in \texttt{almost.py} with step $0.005$. Drawn by \texttt{figs\_almost.py}.}'
* Each old text above occurs once in ch/almost.tex. No sentence of the text changes. 'about $-5.4$' (line 172) stands: figs_almost.json gives −5.3873 at θ = 292.40187 and −5.3784 on the plotting grid, and almost_dip.json of item part6_almost_dip agrees to 1e−8. The background mean and minimum 0.51 and 0.069 at height 300, and 1.59 and 0.34 at height 10^5 (line 127) stand: 0.5143, 0.06904, 1.5876, 0.34205. The kernel heights 2.1 and 6.4 (line 172) stand.
* Not a chapter file and not edited here, since this item may only add files: in the first paragraph of riesz_program/README.md, replace 'ri_riesz.pdf, ag_almost.pdf and ag_cancel.pdf are the notes' fig_riesz.pdf, fig_almost.pdf and fig_cancel.pdf' with 'ri_riesz.pdf is the notes' fig_riesz.pdf, and ag_almost.pdf and ag_cancel.pdf are the notes' fig_almost.pdf and fig_cancel.pdf redrawn by code/figs_almost.py (added for the book, output data/figs_almost.json) with the density written v_alpha/pi'

## For the authors

* fig/ag_almost.pdf and fig/ag_cancel.pdf are now drawn by riesz_program/code/figs_almost.py (added for the book). The notes' figures were drawn from an evaluation of the almost.py sum at step 0.005 that the package did not keep. The script reproduces it exactly, and the old files are kept in misc/repro/reconstruct/part6_ag_figures/.
* Notation in the new panels: the titles write v_α/π and v_0.6/π. The left legend of ag_almost writes the exponents in α, as Theorem ag:thm:almost does; the notes' figure used σ. The dotted marker now sits at the hypothetical height 292.40177 (it was at 292.4000). The Guth and Maynard curve starts at α = 0.7 exactly (it started at 0.70067). Both differences are below a tenth of a point on the page.
* The ag_cancel legends still say 'with a zero at β = 0.7' (and 0.9), as in the notes, although high.py adds the pair β and 1 − β. The caption already says 'a hypothetical pair of zeros', so the legend was kept.
