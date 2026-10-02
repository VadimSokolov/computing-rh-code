# spot_checks_report.md: archive, rerun and book

## 1. Archived report against the regenerated report

Archived report: 60 lines, 4388 bytes. Regenerated report: 95 lines, 6773 bytes.

Identical lines: 60. Changed blocks: 0. Lines only in the regenerated report: 35. Lines only in the archive: 0.

Numeric tokens compared on the common lines (labels included): 176; beyond relative tolerance 1e-9: 0; largest relative difference: 0.0.

The archive equals the first 57 lines of the regenerated report followed by its last 3 lines (blank line, summary header, 'checks failed: 0'): True.

Lines 59 to 93 of the regenerated report are not in the archive (35 lines; sections: Proposition prop:hmfalse (assumption (H) is false), Table tab:symwin (zeros of xi in the doubled windows), Round 2 corrections).


## 2. Regenerated report against the full rerun (logs/spot_checks.out)

Byte identical: True. Lines 95 and 95. Numeric tokens compared: 319; beyond 1e-9: 0; largest relative difference: 0.0.


## 3. Numbers the book quotes from the report

112 quoted values or statements checked: 112 agree at the printed precision, 0 do not or were not found.

A value agrees when the book value is the rerun value correctly rounded at the last digit the book prints (allowing half a unit of the last digit the report prints; the book prints S(gamma_289^-) and gamma_289 truncated in rz3, marked trunc). Line numbers refer to the chapter files in book_snapshot/ (sha256 in comparison.json).

| id | where | quantity | book | rerun | verdict |
|---|---|---|---|---|---|
| R01 | rz4:497 | c, three ways (Table rz:tab:constants) | 2.6860917096 | 2.68609170961 | agrees |
| R01 | rz4:497 | c, three ways (Table rz:tab:constants) | 2.6860917096 | 2.68609170961 | agrees |
| R01 | rz4:497 | c, three ways (Table rz:tab:constants) | 2.6860917096 | 2.68609170961 | agrees |
| R02 | rz4:498 | tau_0, theta(tau_0), max M (Table rz:tab:constants) | 6.2898360 | 6.28983598884 | agrees |
| R02 | rz4:498 | tau_0, theta(tau_0), max M (Table rz:tab:constants) | -3.5309728 | -3.53097282902 | agrees |
| R02 | rz4:498 | tau_0, theta(tau_0), max M (Table rz:tab:constants) | 1.1239436 | 1.12394355932 | agrees |
| R03 | rz4:499 | M(gamma_1^2 -) (Table rz:tab:constants) | 0.5502528 | 0.550252829469 | agrees |
| R04 | rz4:500 | tau_pi (Table rz:tab:constants) | 3.4362182 | 3.43621822609 | agrees |
| R05 | rz4:501 | x_* (Table rz:tab:constants) | 0.0032173906 | 0.00321739055955 | agrees |
| R06 | rz4:502 | F(0.003), F(0.0035), F(0.004) (Table rz:tab:constants) | 0.0366850 | 0.0366850179113 | agrees |
| R06 | rz4:502 | F(0.003), F(0.0035), F(0.004) (Table rz:tab:constants) | -0.0442285 | -0.0442284923318 | agrees |
| R06 | rz4:502 | F(0.003), F(0.0035), F(0.004) (Table rz:tab:constants) | -0.1145624 | -0.11456244439 | agrees |
| R07 | rz4:503 | P(0.004) (Table rz:tab:constants) | 1.99e-13 | 1.98858765094e-13 | agrees |
| R08 | rz4:504 | min S over the first 300, at gamma_289 (Table rz:tab:constants) | -1.1454808 | -1.14548076977 | agrees |
| R09 | rz4:505 | max S over the first 300 and its ordinate gamma_213 (Table rz:tab:constants) | 1.0975638 | 1.09756376206 | agrees |
| R09 | rz4:505 | max S over the first 300 and its ordinate gamma_213 (Table rz:tab:constants) | 415.4552 | 415.455215 | agrees |
| R10 | rz4:464 | w(0.1), Table rz:tab:num (the script computes the zeros side only) | 0.8625096 | 0.86250960302 | agrees |
| R11 | rz4:465 | w(1), Table rz:tab:num (zeros side only) | 0.4914032 | 0.491403169211 | agrees |
| R12 | rz4:466 | W_0(0.1), Table rz:tab:num (the script computes one side, the integral of theta') | -0.86250960 | -0.862509600916 | agrees |
| R13 | rz4:477 | Prop. rz:prop:zetaright row of Table rz:tab:num, left and right side | -0.587879 | -0.587879014707 | agrees |
| R13 | rz4:477 | Prop. rz:prop:zetaright row of Table rz:tab:num, left and right side | -0.587340 | -0.587339889106 | agrees |
| R14 | rz4:450 | weight of the omitted terms at n=3000 (the book's sign: sqrt(s) minus (1-n^-sqrt(s))/log n) | 1.607 | 1.6071503 (-zr_wt) | agrees |
| R15 | rz4:450 | gap between the two sides of the zetaright row | 5.39e-4 | 0.0005391256 (-zr_err) | agrees |
| R16 | rz4:489 | w at x = 0.001, 0.01, 0.1, 1, 10 | 1.00025 | 1.00025003125 | agrees |
| R16 | rz4:489 | w at x = 0.001, 0.01, 0.1, 1, 10 | 1.00249 | 1.00249472811 | agrees |
| R16 | rz4:489 | w at x = 0.001, 0.01, 0.1, 1, 10 | 0.86251 | 0.86250960302 | agrees |
| R16 | rz4:489 | w at x = 0.001, 0.01, 0.1, 1, 10 | 0.49140 | 0.491403169211 | agrees |
| R16 | rz4:489 | w at x = 0.001, 0.01, 0.1, 1, 10 | 0.21524 | 0.215237315084 | agrees |
| R17 | rz4:489 | c/(2 sqrt pi), limit of x^{3/2} f(x) | 0.75773 | 0.757732481509 | agrees |
| R18 | rz4:88 | kappa(rho_1) to kappa(rho_4), real and imaginary parts (xi form) | -4.762 | -4.762466 | agrees |
| R18 | rz4:88 | kappa(rho_1) to kappa(rho_4), real and imaginary parts (xi form) | -1.944 | -1.94353 | agrees |
| R18 | rz4:88 | kappa(rho_1) to kappa(rho_4), real and imaginary parts (xi form) | -8.673 | -8.672921 | agrees |
| R18 | rz4:88 | kappa(rho_1) to kappa(rho_4), real and imaginary parts (xi form) | -0.627 | -0.6272364 | agrees |
| R18 | rz4:88 | kappa(rho_1) to kappa(rho_4), real and imaginary parts (xi form) | -10.149 | -10.14905 | agrees |
| R18 | rz4:88 | kappa(rho_1) to kappa(rho_4), real and imaginary parts (xi form) | -3.872 | -3.871717 | agrees |
| R18 | rz4:88 | kappa(rho_1) to kappa(rho_4), real and imaginary parts (xi form) | -15.165 | -15.16508 | agrees |
| R18 | rz4:88 | kappa(rho_1) to kappa(rho_4), real and imaginary parts (xi form) | 3.705 | 3.704673 | agrees |
| R18b | rz4:88 | closed form of Theorem rz:thm:kappaexact equals the xi form for rho_1 to rho_4 | closed form = xi form | k1re==c1re and k1im==c1im and k2re==c2re and k2im==c2im and k3re==c3re and k3im==c3im and k4re==c4re and k4im==c4im | agrees |
| R19 | rz4:112 | minus Re kappa(rho_j), caption of Table rz:tab:traj | 4.762 | 4.762466 (-k1re) | agrees |
| R19 | rz4:112 | minus Re kappa(rho_j), caption of Table rz:tab:traj | 8.673 | 8.672921 (-k2re) | agrees |
| R19 | rz4:112 | minus Re kappa(rho_j), caption of Table rz:tab:traj | 10.149 | 10.14905 (-k3re) | agrees |
| R19 | rz4:112 | minus Re kappa(rho_j), caption of Table rz:tab:traj | 15.165 | 15.16508 (-k4re) | agrees |
| R20 | rz4:88 | Re kappa < 0 for the first 100 zeros | no flip below index 101 | min(flip_ks) > 100 | agrees |
| R21 | rz4:86 | first zero with Re kappa > 0 | 213 | min(flip_ks) == 213 | agrees |
| R22 | rz4:130 | number of zeros with Re kappa > 0 among the first 460 | 3 | 3 | agrees |
| R23 | rz4:132 | ordinates of the three zeros with Re kappa > 0 | 415.4552 | 415.455215 | agrees |
| R23 | rz4:132 | ordinates of the three zeros with Re kappa > 0 | 527.9036 | 527.9036416 | agrees |
| R23 | rz4:132 | ordinates of the three zeros with Re kappa > 0 | 650.6687 | 650.6686839 | agrees |
| R24 | rz4:134 | kappa at gamma_213, gamma_289, gamma_379 | 233.96 | 233.9579 | agrees |
| R24 | rz4:134 | kappa at gamma_213, gamma_289, gamma_379 | -2298.91 | -2298.912 | agrees |
| R24 | rz4:134 | kappa at gamma_213, gamma_289, gamma_379 | 604.16 | 604.1578 | agrees |
| R24 | rz4:134 | kappa at gamma_213, gamma_289, gamma_379 | 3074.08 | 3074.083 | agrees |
| R24 | rz4:134 | kappa at gamma_213, gamma_289, gamma_379 | 1164.60 | 1164.599 | agrees |
| R24 | rz4:134 | kappa at gamma_213, gamma_289, gamma_379 | -4296.91 | -4296.907 | agrees |
| R25 | rz4:143 | Re kappa(rho_213), proof of Corollary rz:cor:Tfalse | 233.96 | 233.9579 | agrees |
| R26 | rz4:136 | S just below and above gamma_213 and gamma_289 (max S minus 1, max S, min S, min S plus 1) | 0.098 | 0.09756376206 (maxS-1) | agrees |
| R26 | rz4:136 | S just below and above gamma_213 and gamma_289 (max S minus 1, max S, min S, min S plus 1) | 1.098 | 1.09756376206 | agrees |
| R26 | rz4:136 | S just below and above gamma_213 and gamma_289 (max S minus 1, max S, min S, min S plus 1) | -1.145 | -1.14548076977 | agrees |
| R26 | rz4:136 | S just below and above gamma_213 and gamma_289 (max S minus 1, max S, min S, min S plus 1) | -0.145 | -0.14548076977 (minS+1) | agrees |
| R27 | rz4:136 | midpoints of the jumps at gamma_213 and gamma_289 | 0.598 | 0.59756376206 (maxS-mpf(1)/2) | agrees |
| R27 | rz4:136 | midpoints of the jumps at gamma_213 and gamma_289 | -0.645 | -0.64548 | agrees |
| R28 | rz4:136 | number of the first 460 zeros with S(gamma^-) < -1 | 7 | 7 | agrees |
| R29 | rz4:136 | indices of those zeros | 127, 196, 233, 289, 368, 380, 401 | low_ks == [127, 196, 233, 289, 368, 380, 401] | agrees |
| R30 | rz4:136 | their midpoints S(gamma^-)+1/2 | -0.506 | -0.5063 | agrees |
| R30 | rz4:136 | their midpoints S(gamma^-)+1/2 | -0.505 | -0.50518 | agrees |
| R30 | rz4:136 | their midpoints S(gamma^-)+1/2 | -0.537 | -0.53724 | agrees |
| R30 | rz4:136 | their midpoints S(gamma^-)+1/2 | -0.556 | -0.55649 | agrees |
| R30 | rz4:136 | their midpoints S(gamma^-)+1/2 | -0.538 | -0.53825 | agrees |
| R30 | rz4:136 | their midpoints S(gamma^-)+1/2 | -0.514 | -0.51375 | agrees |
| R31 | rz4:136 | largest midpoint modulus with Re kappa < 0 (the script scans negative midpoints only) | 0.556 | 0.556493 (-mid_neg_max) | agrees |
| R31b | rz4:136 | it is attained at gamma_368 | 368 | abs(mids[368] - mid_neg_max) < mpf('1e-5') | agrees |
| R32 | rz4:136 | midpoint at gamma_213 (the script does not print the one at gamma_379) | 0.598 | 0.59756376206 (maxS-mpf(1)/2) | agrees |
| R33 | rz4:199 | zeros of xi with 1 < gamma < 80 and 80 < gamma < 140 (Table rz:tab:symwin) | 21 | 21 | agrees |
| R33 | rz4:199 | zeros of xi with 1 < gamma < 80 and 80 < gamma < 140 (Table rz:tab:symwin) | 27 | 27 | agrees |
| R34 | rz4:201 | zeros of xi with 140 < gamma < 280 and 390 < gamma < 430 (Table rz:tab:symwin) | 78 | 78 | agrees |
| R34 | rz4:201 | zeros of xi with 140 < gamma < 280 and 390 < gamma < 430 (Table rz:tab:symwin) | 26 | 26 | agrees |
| R35 | rz4:210 | first and last index in 390 < gamma < 430, and the count | 196 | 196 | agrees |
| R35 | rz4:210 | first and last index in 390 < gamma < 430, and the count | 221 | 221 | agrees |
| R35 | rz4:210 | first and last index in 390 < gamma < 430, and the count | 26 | 26 | agrees |
| R36 | rz4:210 | gamma_195 and gamma_222 | 388.85 | 388.84613 | agrees |
| R36 | rz4:210 | gamma_195 and gamma_222 | 430.33 | 430.32875 | agrees |
| R40 | rz3:72 | S(gamma_289^-) and gamma_289, printed truncated | -1.14548 | -1.14548076977 | agrees |
| R40 | rz3:72 | S(gamma_289^-) and gamma_289, printed truncated | 527.9036 | 527.9036416 | agrees |
| R41 | rz3:81 | S(gamma_127^-) = midpoint minus 1/2 | -1.0063 | -1.0063 (mid127-mpf(1)/2) | agrees |
| R42 | rz3:81 | most negative S(gamma^-) among the first 460 | at gamma_289 | min(mids, key=lambda k: mids[k]) == 289 | agrees |
| R43 | rz3:89 | f(0+) | 0.36509 | 0.365088627202 | agrees |
| R44 | rz3:91 | c/(2 sqrt pi) | 0.75773 | 0.757732481509 | agrees |
| R45 | rz3:126 | c/(2 sqrt pi) | 0.75773 | 0.757732481509 | agrees |
| R46 | rz3:101 | f'(0+) | 0.09127 | 0.0912721568006 | agrees |
| R47 | rz3:126 | f(0+) | 0.365089 | 0.365088627202 | agrees |
| R48 | rz3:161 | maximum of f_H and where (Prop. rz:prop:hmfalse) | 1.11125 | 1.111254 | agrees |
| R48 | rz3:161 | maximum of f_H and where (Prop. rz:prop:hmfalse) | 0.09765 | 0.0976518 | agrees |
| R49 | rz3:200 | negative b_n for 2 <= n <= 60 | 42 | 42 | agrees |
| R49 | rz3:200 | negative b_n for 2 <= n <= 60 | 42 | 42 | agrees |
| R49 | rz3:200 | negative b_n for 2 <= n <= 60 | 42 | 42 | agrees |
| R49 | rz3:200 | negative b_n for 2 <= n <= 60 | 42 | 42 | agrees |
| R50 | rz3:305 | min of Re zeta'/zeta over (0,60] (Prop. rz:prop:onlyhalf discussion) | -1.049 | -1.04918 | agrees |
| R50 | rz3:305 | min of Re zeta'/zeta over (0,60] (Prop. rz:prop:onlyhalf discussion) | -0.995 | -0.9949883 | agrees |
| R50 | rz3:305 | min of Re zeta'/zeta over (0,60] (Prop. rz:prop:onlyhalf discussion) | -0.867 | -0.8666873 | agrees |
| R50 | rz3:305 | min of Re zeta'/zeta over (0,60] (Prop. rz:prop:onlyhalf discussion) | -0.769 | -0.7688592 | agrees |
| R50 | rz3:305 | min of Re zeta'/zeta over (0,60] (Prop. rz:prop:onlyhalf discussion) | -0.721 | -0.721326 | agrees |
| R51 | rz3:35 | x_* (caption of Figure rz:fig:w) | 0.00322 | 0.00321739055955 | agrees |
| R52 | rz3:7 | c | 2.68609 | 2.68609170961 | agrees |
| R53 | rz3:30 | minus W_0(0.004) = minus F(0.004)/(4 sqrt(0.004 pi)) by Lemma rz:lem:W0 | exceeds 0.25 | -F004/(4*sqrt(pi*mpf('0.004'))) > mpf('0.25') | agrees |
| R60 | rz2:7 | b | 0.860973 | 0.860972775375 | agrees |
| R61 | rz2:264 | c | 2.68609 | 2.68609170961 | agrees |
| R70 | rz1:228 | tau_0 | 6.28984 | 6.28983598884 | agrees |
| R71 | rz1:232 | c | 2.68609 | 2.68609170961 | agrees |
| R72 | rz1:234 | max M | 1.12394 | 1.12394355932 | agrees |
| R73 | rz1:234 | M(gamma_1^2 -) | 0.55025 | 0.550252829469 | agrees |
| R80 | ch12:153 | c | 2.6860917 | 2.68609170961 | agrees |
| R90 | appA:75 | c, constants table of Appendix A (cites spot_checks_report.md) | 2.6860917096 | 2.68609170961 | agrees |

## 4. Reference values printed by spot_checks.py that are no longer in the book

The script prints 55 reference strings; 4 of them (b_2, b_3, b_4, b_6) are its own evaluations of formulas in the text, and the other 51 were looked up in the chapters. 47 occur in the current text; 4 do not:

- zr_trunc: rerun -0.587339889106, script prints "paper: -0.587379"; "0.587379" does not occur in rz1 to rz4, appA or ch12
- zr_err: rerun -0.0005391256, script prints "row's printed gap: -5.00e-4"; "5.00\times10^{-4}" does not occur in rz1 to rz4, appA or ch12
- zr_wt: rerun -1.6071503, script prints "round 1 assumed 1"; wording of a reply to a referee; no number to look for
- grid_worst: rerun -0.729142, script prints "paper: below -0.729"; "0.729" does not occur in rz1 to rz4, appA or ch12

## 5. Monte Carlo lines

- h=0.1: Monte Carlo 0.891195, closed form 0.891613, difference -0.000418, standard error 0.000695, z = -0.601; the script's tolerance 0.004 is 5.75 standard errors.
- h=1: Monte Carlo 0.224115, closed form 0.221429, difference 0.00269, standard error 0.000928, z = 2.89; the script's tolerance 0.004 is 4.31 standard errors.
- h=2: Monte Carlo 0.038555, closed form 0.0378576, difference 0.000697, standard error 0.000427, z = 1.63; the script's tolerance 0.004 is 9.37 standard errors.
