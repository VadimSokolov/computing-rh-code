Layout in this repository (September 2026). The authors sent this package after the merged basepoint one paper, and then two revisions of notes/riesz_thorin.tex; the first version, integrated first, is kept as misc/reviews/riesz/riesz_thorin_v1.tex. The notes are integrated as Chapters ch:riesz (notes/riesz_thorin.tex without its section on the moments of zeta, and work package 1 of notes/strategy.tex) and ch:almost (notes/almost_ggc.tex, work packages 3 and 4, and the section of riesz_thorin.tex on the moments of zeta, the Barnes G function and zero density) of Part VI; the rest of the strategy note goes to Chapter ch:open of Part VII. The figures are in fig/: ri_riesz.pdf is the notes' fig_riesz.pdf; ag_almost.pdf and ag_cancel.pdf are the notes' fig_almost.pdf and fig_cancel.pdf redrawn by code/figs_almost.py (added for the book, output data/figs_almost.json) with the density written v_alpha/pi; and code/figs_strategy.py splits the three panels of fig_strategy.pdf into ri_filters.pdf, ag_margins.pdf and st_strategy.pdf. The zero lists in zeros/ are the same files as in basepoint_one/xi. Added for the book: code/v1min.py, which certifies that the Thorin density v_1 of T_1 attains its minimum only at the origin (output data/v1min.json, Proposition bp:prop:v1min); code/explicit_bounds.py (output data/explicit_bounds.json), which certifies in ball arithmetic the gaps of the first 1800 zeros, the window count above height 2000 and the other numbers of the proof of Theorem ag:thm:explicit, and the constants of the proof of Theorem sm:thm:strip (3); misc/repro/ks_ak.py, which recomputes the arithmetic factors that code/ks.py hard-codes; code/almost_dip.py (output data/almost_dip.json), which finds every negative local minimum of the density of almost.py at the basepoints 0.6 and 0.7, among them the dip to -5.387 at alpha = 0.7 near the height 292.40 that Chapter ch:almost quotes; code/check_large_twisted.py (output data/check_large_twisted.json), which runs the twisted option of riesz_large.py at x = 1e10 and h = 50 and compares it with the explicit formula over 300 zeros (Section ri:sec:filters); code/ggc_offset.py (output data/ggc_offset.json), which splits the errors of Table ri:tab:ggccoef, computed by ggc_coef.py, into the error of the Weyl tail and the quadrature error; code/rmt_ratios.py (input data/killing.json, output data/rmt_ratios.json), which computes the exact ratios E|Z_N|^{2k}/N^{k^2} for Haar unitary matrices that Chapter ch:almost compares with the Barnes constants; and code/figs_almost.py (output ag_almost.pdf and ag_cancel.pdf, copied to fig/, and data/figs_almost.json), which draws Figures ag:fig:almost and ag:fig:cancel from almost.json, more.json, high.json and the zero lists. Every script of the authors' package was rerun on Hopper (misc/repro/compare_riesz_report.md), and the scripts added for the book ran there too (for the last five, see misc/repro/reconstruct/part6_*/REPORT.md).

The authors' README follows.

# The Riesz condition and a computational strategy for the Thorin path to RH

Nick Polson and Vadim Sokolov, September 2026. Package for Vadim.

## Notes (notes/)
riesz_thorin.pdf   Riesz's criterion and the Thorin representation of 1/zeta(2k) (13 pages)
                   R(x) = -sum (-x)^k/((k-1)! zeta(2k)) = x sum mu(n) n^-2 e^{-x/n^2}; RH iff R(x) = O(x^{1/4+eps});
                   basepoint version: no zeros with Re s > alpha iff R(x) = O(x^{alpha/2+eps});
                   Thorin form: 1/zeta(2k) = factor * E exp(-(2k-1/2)^2 T), T the centre clock (no hypothesis);
                   GGC form: every 1/zeta(2k) is a value of the Laplace transform of the basepoint one clock X_1
                   (a GGC with no hypothesis) at (2k-1)^2, built from the primes by the 2 sin^2 lemma, and factors as
                   a stable 1/2 value times a GGC value; Riesz's criterion = X_alpha is a GGC for every alpha > 1/2;
                   hitting times: each coefficient is a killing probability P(X_1 < E_{(2k-1)^2}); the stable factor is
                   a Brownian hitting probability at the Wald defect level; half Cauchy duality through the sech clock;
                   moments of zeta, Barnes G (Nikeghbali and Yor) and zero density exponents;
                   Keating Snaith moments: arithmetic factors a_k, moments to T = 1000, and the identity
                   log|xi(1/2+it)| = log|xi(1+it)| - pi int_{1/2}^1 rho_a(t) da (moments = exponential moments of Thorin density);
                   why it cannot supply the bound (positive only where Re 2k >= 1).
strategy.pdf       A computational strategy for the Thorin path to RH (6 pages): five work packages.
almost_ggc.pdf     Is there a path to RH? Almost generalised gamma convolutions (9 pages): context for the strategy.

## Key computational findings
1. The classical Riesz function is blind above the first few zeros: zero rho enters with weight
   |Gamma(1 - rho/2) / 2 zeta'(rho)|, which decays like exp(-pi gamma / 4) (fitted slope 0.786).
2. Twisted Riesz functions R_T(x) = x sum mu(n) n^-2 (x/n^2)^{iT} e^{-x/n^2} are still equivalent to RH
   (Mellin transform Gamma(1-s+iT)/zeta(2s)) but see the zeros near height 2T. At T = 50 the weight peaks at
   gamma = 101.3 (0.124 against 2.9e-35 untwisted); the Mobius series matches the explicit formula to 1.2e-8.
   They are computed from mu(n) only, so they check zeros independently of any zeta evaluation.
3. riesz_large.py (segmented Mobius sieve, extended precision, exact tail) matches the explicit formula to about
   2e-9 in R(x)/x^{1/4} for x = 1e8, 1e10, 1e11, 1e12 (42 s at 1e12 in Python). The twisted option matches too.
4. Positivity margins of the Thorin density at midpoints between unfolded zeros follow the circular unitary
   ensemble (medians agree to 1 percent; zeta slightly more rigid in the lower tail).

## Code (code/) and how to run
Run ./setup.sh once (links zeros/ and data/ into code/), then from code/:
  python3 riesz_large.py 1e12            Riesz function at x (default K = 200, N = K sqrt(x))
  python3 riesz_large.py 1e10 200 50     twisted Riesz function R_T with T = 50
  python3 check_large.py                 validation against the explicit formula (60 zeros + trivial terms)
  python3 twisted.py                     twisted Riesz function T = 50: weights, series, explicit formula
  python3 sens.py                        Riesz weights; margin statistics, zeta against CUE
  python3 salem_riesz.py                 Salem kernel |eta(sigma+it)|, Riesz function from the Mobius series
  python3 riesz_thorin.py                Riesz coefficients and R(x) as averages over the centre clock
  python3 ggc_coef.py                    Riesz coefficients in GGC form (Thorin integral, stable factor, Levy form)
  python3 killing.py; python3 killing2.py   killing probabilities, Riesz series from them, half Cauchy duality,
                                         Keating Snaith moments (Haar unitary and beta decomposition)
  python3 ks.py                          a_k, moments of zeta on the line, Selberg variance, Thorin form of log|xi| (about 2 min)
  python3 almost.py, more2.py, high.py   almost GGC experiments (exponents, negative mass, cancellation lemma)
  python3 figs_riesz.py, figs_strategy.py   figures
Requirements: Python 3 with numpy, scipy, matplotlib, mpmath, sympy. The Mobius sieve to N needs about N bytes.

## Suggested next steps (strategy.pdf, work packages)
WP1  Port riesz_large.py to C or numba; run R(x) and R_T(x) for x up to 1e14 then 1e16, T on a grid up to 500;
     compare with the explicit formula over computed zeros; beyond 1e16 compute sum_{n<=N} mu(n)/n^2 by
     Mertens function style algorithms and use quadruple precision.
WP2  Certified Thorin densities rho_alpha at alpha - 1/2 = 1e-1, 1e-2, 1e-3 on windows at heights 1e6 to 1e10
     (Arb via python-flint) and the zero density identity for window certification.
WP3  Margin statistics from Odlyzko's tables near the 1e12th, 1e21st and 1e22nd zeros; CUE extreme value fit;
     convert to the width of the dangerous strip via the cancellation lemma.
WP4  Explicit almost GGC defect tables with explicit zero density estimates.
WP5  Test candidate lemmas numerically before pursuing them (Conrey and Li lesson).
