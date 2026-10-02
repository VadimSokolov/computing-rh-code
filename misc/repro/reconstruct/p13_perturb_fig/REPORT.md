# Reproducibility report: perturb_fig

Written by the main session from the structured result of workflow wf_65af78f3-bcf, because the harness does not let subagents write report files. Status: done.

## Findings

No script in code/ draws fig/perturb.pdf. plots3.py section 8, which appA names, underflows (all 300 values are 0), and a rerun of plots3.py would overwrite the figure with an empty panel (job 1244206). The archived figure plots log10|W_F-W|/W on logspace(-8,-2,400) with ylim(-3000,5), for zeros 1496 and 1497 (1977.1739, 1977.2714), with W the sum over the first 1700 zeros. Recomputed in mpmath at 40 digits and drawn with Matplotlib 3.10.8, its PDF content stream is byte identical to the archived one; only the date differs (jobs 1244301, 1244302). W from the explicit formula would move the first vertex by 0.73 in log10. No printed number changes: the maximum is 3.71e-11 at t=7.13e-7 (3.66e-11 at 7.24e-7 with the exact W), and the change is below 1e-1000 from t=5.88e-4. The authors should add code/perturb.py and data/perturb.json, delete plots3.py section 8, and fix the caption and appA. REPORT.md was not written because the harness blocks report files from subagents; the details are in these fields.

## Files added

* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_perturb_fig/perturb.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_perturb_fig/run_final.slurm
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_perturb_fig/perturb.json
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_perturb_fig/perturb_w1700.pdf
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_perturb_fig/perturb_w1700.png
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_perturb_fig/perturb_exact.pdf
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_perturb_fig/perturb_exact.png
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_perturb_fig/compare_fig.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_perturb_fig/compare_fig.json
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_perturb_fig/pdf_paths.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_perturb_fig/plots3_section8.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_perturb_fig/plots3_section8.pdf
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_perturb_fig/plots3_section8-1.png
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_perturb_fig/old_perturb.pdf
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_perturb_fig/old_perturb-1.png
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_perturb_fig/old_perturb.content.txt
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_perturb_fig/old_perturb.paths.json
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_perturb_fig/g_1_400.npy
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_perturb_fig/g_401_1000.npy
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_perturb_fig/g_1001_1700.npy
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_perturb_fig/code_proposal/perturb.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_perturb_fig/code_proposal/run_proposal.slurm
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_perturb_fig/code_proposal/perturb.json
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_perturb_fig/code_proposal/perturb.pdf
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_perturb_fig/logs/run_final_slurm-1244301.log
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_perturb_fig/logs/run_first16_slurm-1243749.log
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_perturb_fig/logs/proposal_slurm-1244011.log
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_perturb_fig/logs/proposal_slurm-1244302.log

## Hopper jobs

* 1243028: pdf_paths.py read the archived fig/perturb.pdf. It is Matplotlib 3.10.8, created 2026-09-25 18:06:33 UTC, 6.5 by 3.5 in. One curve of log10|W_F-W|/W on numpy.logspace(-8,-2,400), ylim(-3000,5), lw 1.6, default colour #1f77b4. rcParams are those of deconvplot.py, lplots.py and kmono.py (plots3.py would draw #1b4f72). Title (size 9) 'Moving the pair near height 1977 off the line:' / 'relative change in the heat trace'. Matplotlib kept 42 vertices; the curve is cut at the page edge after crossing -3000 at t=1.767e-3.
* 1243499: first run of perturb.py with the venv's Matplotlib 3.11.2. Same data; the axes box moves by 1 to 3 pt.
* 1243688: first compare_fig.py. The archived vertices match the curve with W as the sum over 1700 zeros.
* 1243749: perturb.py plus compare_fig.py with Matplotlib 3.10.8 (private pip --no-deps install in /scratch/vsokolov/rh_book_repro/agents/p13_perturb_fig_mpl/mpl3108; the shared venv is untouched). Content stream identical to the archive. Superseded by 1244301 after a comment edit.
* 1244301 (final, run_final.slurm, 16 cores, 192 s): perturb.py plus compare_fig.py. perturb_w1700.pdf has a content stream byte identical to the archive. All 40 grid vertices agree with the 1700 zero curve to 7.95e-6 in log10 (PDF rounding is 1.7e-5). With W from the explicit formula or the two term expansion, the first vertex (t=1e-8) is off by 0.7324. Maximum |W_F-W|/W: 3.7055e-11 at t=7.134e-7 with W1700; 3.6619e-11 at t=7.242e-7 with the explicit formula and with the two term expansion. W1700 against exact W: 1669.34 vs 9014.80 at 1e-8, 1430.14 vs 2337.82 at 1e-7, 710.770 vs 720.132 at 7e-7, 575.918 vs 577.494 at 1e-6, equal to 1e-16 for t>=1e-5. Relative change 3.70e-1699 at t=1e-3; below 1e-1000 from t=5.882e-4. W_F-W changes sign at t=1.27897e-7 (about 1/(2g^2)); the grid shows only a dip to -13.0 there. The explicit formula reproduces W(0.002), W(0.005), W(0.01) of data/arith.json to 15 digits and the ch12 two term errors 4.0e-5, 1.46e-4, 6.53e-4 at t=1e-5, 1e-4, 1e-3.
* 1243969: first test of the proposed code/perturb.py. It failed before running (no /usr/bin/time on the node).
* 1244011: proposed code/perturb.py on one core, 158 s. Content stream identical to the archive.
* 1244302 (final, run_proposal.slurm): proposed code/perturb.py on one core, 162 s, plus compare_fig.py. Content stream identical. The whole PDF differs from fig/perturb.pdf only in /CreationDate and in the startxref offset it shifts (14641 vs 14647).
* 1244206: plots3_section8.py, which is code/plots3.py lines 92 to 103 verbatim. 300 of 300 values of W_F-W are exactly 0 in double precision, so the panel is empty.

## Changes the book needs

* code/perturb.py (new file): copy misc/repro/reconstruct/p13_perturb_fig/code_proposal/perturb.py (52 lines). It loads g_1_400.npy, g_401_1000.npy and g_1001_1700.npy and picks the pair as plots3.py section 8 does. It evaluates W_F-W at 40 digits with W the sum over 1700 zeros, and finds the maximum, the sign change and the crossing of 1e-1000. It evaluates the explicit formula at t=1e-8, 1e-7, 7e-7, 1e-6 and 1e-5, then writes perturb.json and fig/perturb.pdf, in under 3 minutes on one core. With Matplotlib 3.10.8 it reproduces the archived fig/perturb.pdf byte for byte apart from the date (job 1244302).
* code/plots3.py lines 92 to 103 (section 8): delete them; no new text. No later section uses the names they define. OLD: # 8 perturbation invisibility (Theorem D type) gp=g[(g>1000)&(g<2000)]; k=int(np.argmin(np.diff(gp))); A_,B_=gp[k],gp[k+1]; gg=(A_+B_)/2 tt=np.logspace(-3,-0.5,300) W=np.array([np.sum(np.exp(-g**2*x)) for x in tt]) a=gg**2-1/16-1j*gg/2 dW=np.array([2*np.real(np.exp(-a*x))-np.exp(-A_**2*x)-np.exp(-B_**2*x) for x in tt]) fig,ax=plt.subplots(figsize=(6.5,3.5)) ax.loglog(tt,np.abs(dW)/W+1e-300,lw=1.6) ax.set_xlabel('t'); ax.set_ylabel(r'$|W_F-W|/W$'); ax.set_ylim(1e-40,1) ax.set_title('Moving a pair near height %.0f off the line changes the heat trace by'%gg+'\nat most this relative amount',fontsize=9) fig.tight_layout(); fig.savefig(F+'perturb.pdf'); plt.close() print('pair',A_,B_)
* data/perturb.json (new file): copy misc/repro/reconstruct/p13_perturb_fig/code_proposal/perturb.json (20 KB). It holds the pair and its indices, the 400 point curve, max_rel 3.7055e-11 at t_max 7.134e-7, max_rel_Wexp 3.6619e-11 at 7.242e-7, t_zero 1.279e-7, t_below_1e_minus_1000 5.882e-4, and W_check.
* ch/ch13.tex line 30 (caption of fig:ch13:perturb). OLD: \caption{Relative change $|W_F-W|/W$ of the heat trace when the closest pair of zeros between heights $1000$ and $2000$ (near $1977.22$) is moved to distance $\frac14$ from the line, computed from $1700$ zeros. The change is at most $4\times10^{-11}$, near $t=7\times10^{-7}$, and falls below $10^{-1000}$ by $t=10^{-3}$.} NEW: \caption{Relative change $|W_F-W|/W$ of the heat trace, on a $\log_{10}$ scale, when the closest pair of zeros between heights $1000$ and $2000$ (near $1977.22$) is moved to distance $\frac14$ from the line, computed at $40$ digits with $W$ replaced by its sum over the first $1700$ zeros. The change is at most $4\times10^{-11}$, near $t=7\times10^{-7}$, and falls below $10^{-1000}$ by $t=10^{-3}$. The sum over $1700$ zeros is smaller than $W$, by $1.3$ percent at $t=7\times10^{-7}$ and by a factor of $5.4$ at $t=10^{-8}$, so the curve bounds the change from above; with $W$ from the explicit formula \eqref{eq:ch2:Warith} the maximum is $3.66\times10^{-11}$, at $t=7.2\times10^{-7}$.}
* ch/ch13.tex line 26: delete the AUTHORS comment, which the changes above answer; no new text. OLD: % AUTHORS (review U09 B14, F12): no script in code/ produces fig/perturb.pdf (code/plots3.py, section 8, evaluates t in [10^(-3), 10^(-0.5)], where the change underflows in double precision). The maximum 3.66 x 10^(-11) at t = 7.2 x 10^(-7) was recomputed with W from the two term small time expansion; the sum over the first 1700 zeros falls short of W by about 1 percent at t = 7 x 10^(-7) and by a factor of about 5 at t = 10^(-8). Please add the script and say how W was evaluated at small t.
* ch/appA.tex line 45 (row of tab:app:scripts). OLD: \ref{ch:finite} & \raggedright \texttt{plots3.py} & \raggedright Figure~\ref{fig:ch13:perturb}\tabularnewline NEW: \ref{ch:finite} & \raggedright \texttt{perturb.py} (\texttt{perturb.json}) & \raggedright Figure~\ref{fig:ch13:perturb}\tabularnewline
* ch/appA.tex line 21. OLD: Nine scripts load the first $1700$ ordinates of the zeros from NEW: Ten scripts load the first $1700$ ordinates of the zeros from (perturb.py is one more; plots3.py still loads them). Recount if other round 2 scripts that load g_1_400.npy, such as hankel.py or spacing.py, are adopted.
* README.md line 5 (optional). OLD: code/ holds the scripts behind every table and figure (plots*.py, deconvplot.py, lehmerflowplot.py draw the figures); data/ holds their outputs. NEW: code/ holds the scripts behind every table and figure (plots*.py, deconvplot.py, lehmerflowplot.py and perturb.py draw the figures); data/ holds their outputs.
* Unchanged: fig/perturb.pdf (the redraw is identical, so it was not replaced), ch/ch13.tex line 24 (pair near 1977.22, smallest normalised spacing), ch/ch15.tex line 54 (0.089; 1977.17 and 1977.27). A grep of ch/*.tex for 1977, 10^{-11}, 10^{-1000}, 7\times10^{-7}, 3.66 and perturb finds no other sentence that quotes the figure.

## For the authors

* REPORT.md was not written: the harness refused report .md files from this subagent. Its content (what was wrong, the jobs, the comparison, the exact changes) is in this structured output.
* Choose how W is shown. The proposal keeps the archived choice (W is the sum over 1700 zeros, an upper bound on the change) and says so in the caption. The alternative is W from the explicit formula: misc/repro/reconstruct/p13_perturb_fig/perturb_exact.pdf, visually identical, first vertex 0.73 lower in log10, maximum 3.66e-11 at 7.2e-7. That would need the caption to say so.
* Byte identity needs Matplotlib 3.10.8. The venv's 3.11.2 gives the same data, but the axes box moves by 1 to 3 pt.
* Output path convention: the proposed code/perturb.py writes fig/perturb.pdf, like plots.py and plots2.py. plots3.py, deconvplot.py and kmono.py write book/fig/. Pick one.
* If the authors would rather keep everything in plots3.py, the body of the proposal can replace section 8, but ax.semilogx needs color='#1f77b4' because plots3.py's colour cycle starts with #1b4f72. That variant was not run.
* Optional: the vertical range of 3000 hides what the caption quotes. For t below 1e-5 the curve (-10.4 to -12.9) sits on the top axis, and W_F-W changes sign at t=1.28e-7. An inset for t in [1e-8,1e-5] would show both; that would be a new figure.
* The archived figure came from a script that is not in the archive. The file was written 25 s after fig/kernel.pdf and uses the deconvplot.py style of rcParams.
