# Reproducibility report: hankel

Written by the main session from the structured result of workflow wf_36eed84b-56b, because the harness does not let subagents write report files. Status: reproduced.

## Findings

No writer for data/hankel.json exists anywhere: not in the book, its git history or other projects under ~/Dropbox/papers. The file holds, for t = 0.01 and 0.05, the ratio c0c2/c1^2 and the normalised even and odd leading Hankel minors for N = 2..5 (Table tab:ch12:hankel). I reconstructed the writer as hankel.py (mpmath at 60 digits, 1700 float64 zeros from g_*.npy, nstr at 4 and 10 digits). Its output is byte identical to the archive (md5 0e2cd78b...): 18 of 18 numbers agree, 0 differences. Arb enclosures (relative radius under 5e-71) show every archived string and all 18 book table entries are correctly rounded. The worst deviation, 4.1e-4, is plain four digit rounding. The result does not change with 80 zeros, with 25 digit arithmetic, or with float64 zeros (3.5e-15). REPORT.md was not written because the harness refuses report files from subagents; its content is below and in hankel_check.json.

## Scripts

* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/hankel/hankel.py
* /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/hankel/hankel_check.py

## Hopper jobs

* 1241775: hankel.py (reconstruction), hop058, intel, 2 s, COMPLETED, in /scratch/vsokolov/rh_book_repro/agents/reconstruct_hankel
* 1241966: hankel_check.py (Arb enclosures, precision and truncation scans, comparisons), hop062, intel, 16 s, COMPLETED, same folder; outputs fetched: hankel.json, hankel_check.json, slurm-1241775.log, slurm-1241966.log

## Changes the book needs

* No number in the book needs correction: all 18 entries of Table tab:ch12:hankel (ch/ch12.tex lines 55 to 58) are the correctly rounded true values (e.g. line 58: 1.79472e-5, 1.22490666e-13, 1.71497467e-26, 3.81297e-42); the ratios are 1.13465478943 and 1.00000811308.
* ch/appA.tex line 23: once a writer is added, remove hankel.json from the list of data files that no script writes.
* ch/appA.tex line 44 (Table tab:app:scripts, row for ch:heat): add hankel.py (hankel.json) to the scripts column and Table tab:ch12:hankel to the right column.
* Optional, ch/ch12.tex line 65 and ch/ch22.tex line 60: the powers 10^-5, 10^-14, 10^-27, 10^-44 (and 10^-44 by size five) are those of the even block; the odd block gives 10^-5, 10^-13, 10^-26, 10^-41 (3.81e-42); write 'the even block'.
* Optional, ch/ch12.tex line 47: the zeros are double precision ordinates and the sixty digits are the arithmetic; 80 zeros and 25 digit arithmetic already give every printed digit (Exercise at line 162: the first row does not change with 100 zeros, not even by 5e-71 relative).

## For the authors

* Add a writer for data/hankel.json to code/ (e.g. copy misc/repro/reconstruct/hankel/hankel.py, whose header marks it as a reconstruction, or the authors' own script if it turns up) and record the command 'python3 hankel.py', run beside g_1_400.npy, g_401_1000.npy and g_1001_1700.npy. These zero files are not in code/ (misc/repro/zeros.py writes them; identical copies are in riesz_program/zeros/). The same holds for arith.py, kmono.py and plots3.py.
* The three digit table cannot be regenerated from the stored four digit strings. The stored values 1.795e-5, 1.225e-13 and 1.715e-26 round half up to 1.80e-5, 1.23e-13 and 1.72e-26, while the book correctly prints 1.79e-5, 1.22e-13 and 1.71e-26. The table was rounded from more digits than the file keeps; consider storing six digits (mp.nstr(x, 6)).
* P17 pinpoints (round 1 U09 item B20): Theorem 7.1, Numerical observation 7.2, Theorem 7.3 and Theorem 7.4 do resolve, with matching subjects, in the September 2026 draft ~/Dropbox/papers/thorin/incoming/2026-09-14/thorin_replacement-from-rh_lamperti_note/paper/thorin_xi.tex (Section 7), but not in arXiv:1708.02653v25. Observation 7.2 there uses 120 zeros, and its 10 printed values (t = 0.01 minors, ratios 1.1347 and 1.0000081) agree with the enclosures, so the claim at ch/ch12.tex line 47 holds for that draft. Cite that version or P26. Only these four numbers were checked.
* REPORT.md is missing: the harness refused report files from this subagent. The orchestrator should write /Users/vsokolov/Dropbox/papers/computing_rh_book/misc/repro/reconstruct/hankel/REPORT.md from this result. The full number by number comparison (Arb midpoints to 20 digits, radii, precision scan, truncation bounds, book table, round 1 values, P17 observation) is in misc/repro/reconstruct/hankel/hankel_check.json.
