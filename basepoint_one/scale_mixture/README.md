# The reciprocal Xi function on the real axis: normal scale mixture, Thorin form, extended GGCs in the strip

Nick Polson and Vadim Sokolov, September 2026.

scale_mixture_xi.tex / .pdf   the note (8 pages); compile with pdflatex twice (needs fig_strip.pdf)

Main points
  xi(a)/xi(a+|u|) is the characteristic function of Y_a = sqrt(2 X_a) Z for every a >= 1/2, unconditionally.
  The sine squared lemma turns its Levy form in |u| into the Thorin form exp(-int log(1+u^2/t^2) rho_a(t) dt).
  Y_a is an extended GGC iff rho_a >= 0: unconditionally for a >= 1, iff a zero free half plane in the strip,
  iff RH at the centre (then Y = sum L_k / gamma_k, Laplace variables). Unconditionally an extended GGC up to
  a remainder of size log H / H, H = 3e12. Cauchy component of scale xi'/xi(a); density tail b_a/(pi y^2).
  Poisson trace K(y) = int y (4 pi t^3)^(-1/2) exp(-y^2/4t) W(t) dt; RH iff K is completely monotone.

Scripts (run from this folder; zeta ordinates g_*.npy included)
  comp.py     Thorin densities, Thorin form check, Poisson traces, Monte Carlo of Y_a   -> comp.json
  improve.py  Poisson trace as subordinated heat trace, positivity scan of rho_1          -> improve.json
  strip.py    strip checks: Thorin form and product over zeros at a = 0.6, 0.75, remainder bound -> strip.json
  dens2.py    densities of Y_a by Fourier inversion (updates strip.json)
  figs2.py    figure fig_strip.pdf
Requirements: Python 3 with numpy, scipy, matplotlib, mpmath.
