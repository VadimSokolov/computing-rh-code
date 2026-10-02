# Reproducibility report: kappa_cell

Written by the main session from the structured result of workflow wf_65af78f3-bcf, because the harness does not let subagents write report files. Status: done.

## Findings

5.43665 is right. I computed κ_cell from its definition (Ξ′/Ξ at the cell boundaries, Arb, certified zeros): 5.4366509594507135570609604 ± 2.3×10⁻⁴⁵. The book's Pólya route agrees to 2×10⁻⁴⁸. The other numbers of ch15:22 are also right: 5.415 (spacing.json, reproduced bit for bit), observed slope 5.43664 = κ_cell − 0.527ε², 0.4 percent from the zeros above γ_1700 (0.0215569; the smooth tail beyond T* recovers κ_cell to 1.3×10⁻⁹), 0.3981063, 0.3347, −0.225 and +0.250, and the end cells. The computation belongs in code/plots.py beside XipXi200: a tested six line block that writes kappa_cell and cell_leaks to data/extra.json. Side finding: in tab:tv the cell TV at ε = 1 and 0.5 should read 4.992 and 2.655 (np.interp error in code/analyse.py). The harness refused REPORT.md, so its content is in these fields.

## Files added

* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kappa_cell/kappa_cell.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kappa_cell/kappa_cell.json
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kappa_cell/interp_check.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kappa_cell/interp_check.json
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kappa_cell/plots_extra_patched.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kappa_cell/run_plots_extra.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kappa_cell/extra_patched.json
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kappa_cell/compare_extra.txt
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kappa_cell/summary.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kappa_cell/summary.txt
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kappa_cell/probe_api.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kappa_cell/probe_speed.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kappa_cell/run_hopper.sh
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kappa_cell/logs/main_run1.log
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kappa_cell/logs/main_run2.log
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kappa_cell/logs/interp_run1.log
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kappa_cell/logs/interp_run2.log
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kappa_cell/logs/interp_run3.log
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kappa_cell/logs/plots_extra_run1.log
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kappa_cell/logs/summary_run1.log
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kappa_cell/logs/probe_api_slurm-1243193.log
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/p13_kappa_cell/logs/probe_speed_slurm-1243285.log

## Hopper jobs

* 1243193 probe_api.py: python-flint 0.9 API probe (acb_series zeta, acb.zeta_zeros)
* 1243285 probe_speed.py: speed of acb.zeta_zeros up to the 100000th zero (not needed in the end)
* 1243541 kappa_cell.py first run, superseded by 1243696 (its check against part1.json read the references at 15 digits; every other number identical)
* 1243696 kappa_cell.py: kappa_cell in Arb at exact and float64 cells, the Polya route block, truncated sums with the smooth tail, exact d_cell(eps) at 10 resolutions (kappa_cell.json, logs/main_run2.log)
* 1243703 interp_check.py: failed on dictionary keys (numpy rounding)
* 1243727 interp_check.py: linear interpolation diagnosis of analysis.json tv_cells
* 1243743 interp_check.py with cubic Hermite interpolation (interp_check.json)
* 1243767 run_plots_extra.py with plots_extra_patched.py: the patched plots.py section writing extra.json (extra_patched.json, compare_extra.txt)
* 1243806 summary.py: every quoted number, rounded as printed (summary.txt)

## Changes the book needs

* ch/ch15.tex line 22, first two sentences. OLD: Truncated to the first $1700$ ordinates, the zero sum gives $5.415$ for $\kappa_{\rm cell}$, against the observed slope $d_{\rm cell}(0.005)/0.005=5.437$. Since the remainder is $O(\varepsilon^3)$, the observed slope equals $\kappa_{\rm cell}$ up to $O(\varepsilon^2)$, and the difference of $0.4$ percent is the contribution of the zeros above $\gamma_{1700}$, which the form in $\Xi'/\Xi$ includes. NEW: Truncated to the first $1700$ ordinates, the zero sum gives $5.415$, against $\kappa_{\rm cell}=5.43665$ from the form in $\Xi'/\Xi$ and the observed slope $d_{\rm cell}(0.005)/0.005=5.43664$. Since the remainder is $O(\varepsilon^3)$, the observed slope equals $\kappa_{\rm cell}$ up to $O(\varepsilon^2)$, here $1.3\times10^{-5}$, and the difference of $0.4$ percent between $5.415$ and $5.43665$ is the contribution of the zeros above $\gamma_{1700}$, which the form in $\Xi'/\Xi$ includes. The rest of line 22 stays: 0.3981063, 0.3347, −0.225, +0.250 and the end cells claim were recomputed and are right.
* ch/ch15.tex line 22, optional, once extra.json holds cell_leaks: in the last sentence replace `(Figure~\ref{fig:ch15:cells}), with the two largest` by `(Figure~\ref{fig:ch15:cells}), and in the form in $\Xi'/\Xi$ from $-0.2243\,\varepsilon$ to $+0.2513\,\varepsilon$, with the two largest` (exact range −0.22427 at cell 72 to +0.25128 at cell 79).
* ch/ch15.tex lines 24 and 25: once code/plots.py holds the block, delete the comment `% AUTHORS (review U10 B15): the value 5.43665 of the cell constant kappa_cell from Xi'/Xi at the 79 cell boundaries is in no data file; add it to a script if the text is to quote it.` and the empty line after it.
* ch/ch15.tex line 29, caption of fig:ch15:cells. OLD: Their absolute sum, $5.415$, falls short of the observed slope $5.437$ because the zeros above $\gamma_{1700}$ are omitted. NEW: Their absolute sum, $5.415$, falls short of $\kappa_{\rm cell}=5.43665$ because the zeros above $\gamma_{1700}$ are omitted.
* ch/appA.tex line 47, row of ch:spacing in tab:app:scripts. OLD: \ref{ch:spacing} & \raggedright \texttt{plots3.py} & \raggedright Figure~\ref{fig:ch15:cells}\tabularnewline NEW: \ref{ch:spacing} & \raggedright \texttt{plots.py} (\texttt{extra.json}), \texttt{plots3.py} & \raggedright the cell constant $\kappa_{\rm cell}$; Figure~\ref{fig:ch15:cells}\tabularnewline
* ch/appA.tex tab:app:constants, new row after line 76 ($\Xi'(200)/\Xi(200)$ & $-1.2506879$ & \texttt{extra.json}\\): $\kappa_{\rm cell}$, the cell constant of the $79$ cells below height $200$ & 5.43665 & \texttt{extra.json}\\
* No change needed: ch15:3, ch07:92, ch03:149, ch01:47 (5.44 is κ_cell to three digits; 5.44/79 = 0.069); ch07:92 range [0.99888,1.00126] at ε = 0.005 (exact [0.9988786, 1.0012564]); ch07:79, ch09:97, ch09:108 (0.3981063, 1.2506879). fig/cells.pdf needs no change: the exact cell errors exceed the plotted truncated ones by 3.16×10⁻⁴ per unit of cell length, at most 0.0056.
* Side finding, only together with the code/analyse.py fix in the author flags: ch/ch07.tex line 33, Table tab:tv cell TV at ε = 1, `4.995` becomes `4.992`; line 34, at ε = 0.5, `2.656` becomes `2.655` (exact d_cell 4.9919803 and 2.6551924; archived tv_cells 4.9945013 and 2.6556972). The other six entries (1.351, 0.5431, 0.2718, 0.1087, 0.05437, 0.02718) are right.

## For the authors

* REPORT.md is missing because the harness refuses report .md files from subagents. The report is in these fields, and every quoted number is in misc/repro/reconstruct/p13_kappa_cell/summary.txt (job 1243806) and kappa_cell.json.
* code/plots.py: after line 104 (out=dict(XipXi200=XiprimeOverXi, pred_coeff=-XiprimeOverXi/np.pi)) insert: # first order cell constant of ch07 and ch15: Xi'/Xi at the cell boundaries b_0=0, b_k=(g_k+g_{k+1})/2, b_79=200 bc=[0.0]+[float((Z[k]+Z[k+1])/2) for k in range(78)]+[200.0] Lp=[] for t in bc: xa,xb=xi_pair(S,acb(0,t)); Lp.append(float((-(xb.imag))/xa.real)) leaks=[-(Lp[k+1]-Lp[k])/np.pi for k in range(79)] out['kappa_cell']=float(np.sum(np.abs(leaks))); out['cell_leaks']=leaks Then rerun plots.py to regenerate data/extra.json. I ran this block verbatim (plots_extra_patched.py, job 1243767). XipXi200, pred_coeff and prime are unchanged. kappa_cell = 5.436650959450781. cell_leaks: sum 0.3981063280907891, min −0.224273 (cell 72), max 0.251279 (cell 79), each within 1.1×10⁻¹³ of Arb.
* Side finding (code/analyse.py, outside this item): tv_cells reads M_eps at the cell boundaries by np.interp on grids of step 0.125 (ε = 1) and 0.0625 (ε = 0.5). The error is up to 1.1×10⁻⁴ per boundary, while the grid values themselves are exact to 5×10⁻¹³, so tab:tv prints 4.995 and 2.656 where the exact values are 4.992 and 2.655. Fix: line 3 `from scipy.interpolate import CubicSpline, CubicHermiteSpline`; line 44 ` Mc = CubicHermiteSpline(th, M, F/PI)(cells)`. Tested: it reproduces the exact d_cell to 8×10⁻⁷ at ε = 1 and 5×10⁻⁸ at 0.5, and better at smaller ε. Then rerun analyse.py and the plotting scripts (fig/rates.pdf and fig/bench.pdf change invisibly) and apply the tab:tv change.
* data/spacing.json still has no script in code/ (appA:23). misc/repro/reconstruct/spacing/spacing.py reproduces it bit for bit, and this run confirmed it again: coef and all 79 leaks are identical.
