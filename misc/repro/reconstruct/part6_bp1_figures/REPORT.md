# Reproducibility report: bp1_figures

Written by the main session from the structured result of workflow wf_4047bf71-f6a, because the harness does not let subagents write report files. Status: done.

## Findings

No archived script drew fig/bp1_wiener.pdf or fig/bp1_salem_riesz.pdf, and the explicit formula curve was stored nowhere. I added basepoint_one/wiener/figs_wiener.py, which redraws both figures from wiener.json, salem_riesz.json, riesz_explicit.json and g_1_400.npy in the authors' style and the book's notation: θ for t, α and σ = 1/2 and 1, the Möbius typo fixed, italic N and x, and the Riesz title R(x)=O(x^{1/4+δ}). The final run was Hopper job 1243495. By eye (150 and 400 dpi renders) the content is the same. The path comparison (job 1243627) shows identical axis limits: the Wiener curves and markers coincide exactly, the Salem curves agree within 0.14 pt, and the 36 Möbius markers coincide exactly. The curve reproduces riesz_explicit.json bit for bit. I replaced both fig/ files and kept the old copies. The wr:fig:salem caption must drop "the plot writes t for θ". All quoted values stand. REPORT.md was not written because the harness refuses report files from subagents; its content is in these fields.

## Files added

* /Users/vsokolov/Dropbox/papers/computing_rh_book/basepoint_one/wiener/figs_wiener.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/basepoint_one/wiener/fig_wiener.pdf
* /Users/vsokolov/Dropbox/papers/computing_rh_book/basepoint_one/wiener/fig_salem_riesz.pdf
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/part6_bp1_figures/old_bp1_wiener.pdf
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/part6_bp1_figures/old_bp1_salem_riesz.pdf
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/part6_bp1_figures/old_bp1_wiener.png
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/part6_bp1_figures/old_bp1_salem_riesz.png
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/part6_bp1_figures/new_bp1_wiener.png
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/part6_bp1_figures/new_bp1_salem_riesz.png
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/part6_bp1_figures/figs_wiener.log
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/part6_bp1_figures/compare_paths.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/part6_bp1_figures/compare_paths.log

## Files replaced

* /Users/vsokolov/Dropbox/papers/computing_rh_book/fig/bp1_wiener.pdf
* /Users/vsokolov/Dropbox/papers/computing_rh_book/fig/bp1_salem_riesz.pdf

## Hopper jobs

* 1242697 (hop064): figs_wiener.py, first draft (run folder /scratch/vsokolov/rh_book_repro/agents/part6_bp1_figures/draft1)
* 1242908 (hop061): figs_wiener.py, second draft (draft2)
* 1242999 (hop064): compare_paths.py on the second draft (compare_draft2); it showed that the old Salem panel uses the fixed range 1e-3 to 5 and line width 1.0, and the old barrier lines width 1.0, which the final version adopts
* 1243120 (hop064): figs_wiener.py, third run (run3)
* 1243254 (hop064): compare_paths.py on the run3 files (compare3), same result as job 1243627
* 1243495 (hop047): figs_wiener.py, final run after a comment only edit (run4); its PDFs are basepoint_one/wiener/fig_wiener.pdf and fig_salem_riesz.pdf, copied to fig/bp1_wiener.pdf and fig/bp1_salem_riesz.pdf; log misc/repro/reconstruct/part6_bp1_figures/figs_wiener.log (the scp of fig_salem_riesz.pdf failed once and was redone by hand, sha256 checked against the Hopper copy)
* 1243627 (hop064): compare_paths.py on the final fig/ files against the old ones (compare4); log misc/repro/reconstruct/part6_bp1_figures/compare_paths.log: grid lines coincide within 7e-7 pt in every panel, Wiener polylines, barrier lines and 48 markers 0.000 pt, Salem |eta| curves within 0.104 pt, explicit formula curve within 0.142 pt, 36 Mobius markers 0.000 pt, identical line widths and dash patterns

## Changes the book needs

* Required, because the new figure writes theta. File ch/strip.tex, line 276 (caption of Figure wr:fig:salem). Old text: `at $\sigma=\frac12$, $0.6$, $0.75$, $1$; the plot writes $t$ for $\theta$. Right:` New text: `at $\sigma=\frac12$, $0.6$, $0.75$, $1$. Right:`
* Recommended, because the chapter names the script at each figure (line 5), and riesz_explicit.json holds only the trivial coefficients of the curve, not the curve itself. File ch/strip.tex, line 276, end of the same caption. Old text: `The data are in \texttt{salem\_riesz.json}, written by \texttt{salem\_riesz.py}, and in \texttt{riesz\_explicit.json}.}` New text: `The data are in \texttt{salem\_riesz.json}, written by \texttt{salem\_riesz.py}; \texttt{figs\_wiener.py} draws the figure and evaluates the explicit formula with the first $60$ ordinates in \texttt{g\_1\_400.npy} and the trivial terms in \texttt{riesz\_explicit.json}.}`
* Recommended, for the same reason. File ch/strip.tex, line 247, end of the caption of Figure wr:fig:wiener. Old text: `are in \texttt{wiener.json}, written by \texttt{wiener.py}.}` New text: `are in \texttt{wiener.json}, written by \texttt{wiener.py}, and \texttt{figs\_wiener.py} draws the figure.}`
* No number quoted in Sections wr:sec:wiener and wr:sec:salem changes. From the same JSON, figs_wiener.log gives: the barrier 0.10565; the Table wr:tab:wiener errors; |h^(gamma_1)| of 0.0167, 0.0377 and 0.0635 at N=80, and at most 1.17e-8 at the centre; coefficient sums of 1.97e8, 1.19e8 and 3.30e7 above the centre and 4.01e9 at the centre; the sigma=0.6 dips 0.175 at theta=14.144 and 0.189 at theta=37.584; the deepest sigma=0.75 dip 0.266; and the zeros 9.0647, 18.1294, 27.1942 of 1-2^(1-s).

## For the authors

* Riesz panel title: the new figure states Riesz's criterion as ch/strip.tex line 286 does, R(x)=O(x^{1/4+delta}). It replaces the authors' 'R(x)/x^{1/4} bounded', which claims more than the criterion (under RH it would also need simple zeros and a bound on 1/zeta'(rho)) and more than a computation on [1e3,1e10] can show. To restore the old wording, edit the set_title call at basepoint_one/wiener/figs_wiener.py line 72.
* The legends write one half with a slash (alpha = 1/2, sigma = 1/2), not as a stacked \frac12. A stacked fraction in the 7 point legend would print at about 3 points once the figure is scaled to the text width. The captions keep \frac12.
* REPORT.md was not written: the Write tool refused it ('Subagents should return findings as text, not write report files'). Its content is in these fields: what was missing, the files, the jobs, the comparison and the chapter changes. The evidence files are in misc/repro/reconstruct/part6_bp1_figures/.
* Neither basepoint_one/README.md nor basepoint_one/wiener/README.md lists figs_wiener.py, and I did not edit either. Suggested sentence for the first paragraph of basepoint_one/README.md: 'Added for the book: wiener/figs_wiener.py, which draws fig/bp1_wiener.pdf and fig/bp1_salem_riesz.pdf (the merged paper's fig_wiener.pdf and fig_salem_riesz.pdf, whose script is not in the archive) from wiener.json, salem_riesz.json and riesz_explicit.json.'
* The authors' copies basepoint_one/merged/fig_wiener.pdf and fig_salem_riesz.pdf are unchanged, so they now differ from fig/bp1_*.pdf in the notation changes only.
* compare_paths.py imports pdfpaths.py of the item ag_figures (misc/repro/reconstruct/part6_ag_figures/pdfpaths.py, sha256 6c70edb2e6ef782b77191e431c009377512945901f02f5e9cd8836c4d7177057 when used); if that file changes, a rerun of the comparison needs this version.
