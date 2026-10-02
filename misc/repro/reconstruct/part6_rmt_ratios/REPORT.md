# Reproducibility report: rmt_ratios

Written by the main session from the structured result of workflow wf_4047bf71-f6a, because the harness does not let subagents write report files. Status: done.

## Findings

Line 237 of ch/almost.tex quotes the exact ratios E|Z_N|^{2k}/N^{k^2}. They are the field ratio of killing2.py, carried to N=10^2..10^4 and to k=3. Only 0.1716 and 1/12 were archived. rmt_ratios.py computes them in five independent ways: the mpmath product from killing2, an exact rational product with its closed form, the beta decomposition moments, Toeplitz determinants from Heine's identity, and Arb balls. All five agree, and the exact routes match for every N up to 10^4. The exact k=3 ratios are 1.184183e-3, 1.50930019e-4, 1.18902795e-4 and 1.16053609e-4. For k=2 they are 0.1716, 0.09019401, 0.084001919 and 0.0834000192. The Barnes constants are 1/12 and 1/8640 = 1.1574074e-4. Every quoted value stands. REPORT.md was not written because the harness blocks report files from subagents, so its content is in these fields.

## Files added

* /Users/vsokolov/Dropbox/papers/computing_rh_book/riesz_program/code/rmt_ratios.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/riesz_program/data/rmt_ratios.json
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/part6_rmt_ratios/slurm-1242657.log

## Hopper jobs

* 1242427: probe of python-flint 0.9.0 (acb.barnes_g, arb.lgamma) and mpmath 1.3.0 on hop064 (Intel), folder agents/part6_rmt_ratios_probe; no result used
* 1242657: rmt_ratios.py with input killing.json on hop064 (partition normal, feature intel), COMPLETED, exit 0, 5 s, folder /scratch/vsokolov/rh_book_repro/agents/part6_rmt_ratios; local script and JSON are byte identical to the Hopper copies (sha256 checked)

## Changes the book needs

* ch/almost.tex line 237: no value changes. 0.1716, 0.0902, 0.0840 and 0.0834 (k=2), 1.18e-3, 1.51e-4, 1.189e-4 and 1.1605e-4 (k=3), 1/12 and 1.1574e-4 are all the exact values correctly rounded. The per value table is in riesz_program/data/rmt_ratios.json (ratios[].quoted, exact_rounded_like_quoted, quoted_ok).
* Recommended provenance edit, values unchanged. The paragraph credits killing2.py, which computes only N=10, 20 and k=1, 2. ch/almost.tex line 237, old: `The exact values show how slowly $\E|Z_N|^{2k}/N^{k^2}$ approaches the Barnes constant:` new: `The exact values, computed by \texttt{rmt\_ratios.py}, show how slowly $\E|Z_N|^{2k}/N^{k^2}$ approaches the Barnes constant:`

## For the authors

* Optional (authors decide): the closed form gives the rate of the slow convergence. ch/almost.tex line 237, old: `against $\mathsf G(4)^2/\mathsf G(7)=1.1574\times10^{-4}$. The Haar estimates` new: `against $\mathsf G(4)^2/\mathsf G(7)=1.1574\times10^{-4}$. For integer $k$ the product telescopes to $\E|Z_N|^{2k}=\prod_{i=0}^{k-1}\binom{N+k+i}{k}\big/\binom{k+i}{k}$, so that $\E|Z_N|^{2k}/N^{k^2}=\frac{\mathsf G(1+k)^2}{\mathsf G(1+2k)}\prod_{i=0}^{k-1}\prod_{m=1}^k\big(1+\frac{i+m}N\big)$ exactly, and the relative excess over the Barnes constant is $k^3/N+O(N^{-2})$, still $0.27\%$ for $k=3$ at $N=10^4$. The Haar estimates`. The output field excess_times_N_over_k3 confirms it: 1.0012 for k=3 and 1.0003 for k=2 at N=10^4.
* The layout paragraph of riesz_program/README.md lists the scripts added for the book (v1min.py, ks_ak.py). code/rmt_ratios.py and its output data/rmt_ratios.json should be added there. It was not edited here, to avoid clashing with concurrent edits.
* misc/repro/reconstruct/part6_rmt_ratios/REPORT.md was not written because the harness refuses report files from subagents. The main session can write it from these fields. The folder holds the job log slurm-1242657.log. The Hopper jobs are 1242427 (probe) and 1242657 (run). Nothing in ch/, main.tex, fig/ or the authors' scripts and JSON was modified, and no figure is involved.
* Outside this item: in ch/almost.tex line 293, the constant c = 0.70 and the predicted variance 1.48 on [500,1000] are in no data file. ks.json has only var 1.5274 and pred_var 0.945, and the two values rest on the hand computation of check M5(g) in misc/reviews/riesz/moments-math.md.
