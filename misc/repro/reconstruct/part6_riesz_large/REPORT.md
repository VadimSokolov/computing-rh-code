# Reproducibility report: riesz_large

Written by the main session from the structured result of workflow wf_4047bf71-f6a, because the harness does not let subagents write report files. Status: done.

## Findings

The quoted value R_50(10^10) = 37.441418 - 26.042497 i came from an agent's run on the laptop (session output tasks/bya04v9k5.output). The laptop is an Apple M2, where numpy longdouble is plain double. On a Hopper Intel node with 80 bit longdouble, riesz(1e10, 200, 50) gives 37.44141600803035 - 26.042494604612816 i (job 1242771, bit identical to the earlier repro job 1226919). That is 3.07e-9 from the explicit formula over 300 zeros; the quoted value is 1.275e-8 away. An Arb recomputation confirms the explicit value in twisted.json to 1.8e-14. A double precision variant run on Intel reproduces the laptop value to 3e-13. The K scan shows the remaining gap is the truncated Mobius tail. Rounding 1/zeta(2+100i) to double adds 2.25e-10. The untwisted control matches check_large.json bit for bit. Line 47 needs 37.441416 - 26.042495 i and 3.1e-9.

## Files added

* /Users/vsokolov/Dropbox/papers/computing_rh_book/riesz_program/code/check_large_twisted.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/riesz_program/data/check_large_twisted.json
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/part6_riesz_large/run.log
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/part6_riesz_large/check_large_twisted.json

## Hopper jobs

* 1242771: check_large_twisted.py on hop062 (Intel x86_64, numpy 2.4.6, longdouble eps 1.08e-19), run folder /scratch/vsokolov/rh_book_repro/agents/part6_riesz_large, 21:16:52 to 21:20:27 on 26 September, COMPLETED
* 1226919 task 1 (earlier Part VI repro run, not mine): riesz_large.py 1e10 200 50 on hop050, gave the same value bit for bit; its log /scratch/vsokolov/rh_book_repro/riesz/logs/riesz_large_tw.out was never archived

## Changes the book needs

* ch/riesz.tex, line 47 (checked 21:21 on 26 September; Section ri:sec:filters, paragraph beginning 'The cost of $R(x)$ or $R_h(x)$'; the old text occurs once). OLD: Its twisted option agrees as well: at $x=10^{10}$ with $h=50$ it gives $R_{50}(x)=37.441418-26.042497\,i$, within $1.3\times10^{-8}$ of \eqref{ri:eq:twexplicit} over $300$ zeros in $|R_h(x)|/x^{1/4}$. NEW: Its twisted option agrees as well: at $x=10^{10}$ with $h=50$ it gives $R_{50}(x)=37.441416-26.042495\,i$, within $3.1\times10^{-9}$ of \eqref{ri:eq:twexplicit} over $300$ zeros in $|R_h(x)|/x^{1/4}$ (\texttt{check\_large\_twisted.py}, output \texttt{check\_large\_twisted.json}).

## For the authors

* REPORT.md was not written. The harness refused the Write of misc/repro/reconstruct/part6_riesz_large/REPORT.md ('Subagents should return findings as text, not write report files'), so the report content is in this output. The data are in riesz_program/data/check_large_twisted.json and the job log is misc/repro/reconstruct/part6_riesz_large/run.log.
* The rest of line 47 stands. 'Extended precision' holds for the archived runs: on Intel the untwisted control is bit identical to check_large.json. The K scan supports the stated order K^(-7/2) of the tail remainder. With the 1/zeta rounding removed, |R_50 - explicit|/x^(1/4) is 1.65e-8, 3.01e-9, 1.8e-10, 3.3e-12 and 8.4e-12 at K = 100, 200, 400, 800 and 1600.
* Optional clarification, ch/riesz.tex line 47. After 'implements this with a segmented sieve and extended precision' add '(numpy's long double, which is the $80$ bit x87 format on Intel and AMD processors but only double precision on Apple silicon)'. On a Mac the script silently returns the old value 37.441418-26.042497i.
* For the authors: riesz_twisted in riesz_large.py rounds 1/zeta(2+2ih) to double, while the untwisted code keeps 6/pi^2 in extended precision as hi+lo. The error this adds to |R_h|/x^(1/4) grows like x^(3/4): 2.25e-10 here, up to about 3e-9 at x=1e10 for an unlucky h, and about 30 times more at 1e12. A port to x=1e16 (work package 1) should split it the same way.
* riesz_program/README.md: its 'Added for the book' sentence could list code/check_large_twisted.py (output data/check_large_twisted.json). I did not edit the README because other agents may be editing it.
