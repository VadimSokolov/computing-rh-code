# Sources for Part VI: the basepoint one and the critical strip

Layout in this repository (September 2026). The authors' merged paper, which combines the five notes below, is in `merged/` (basepoint_one_merged.tex and .pdf, with its figures and its own README). The scripts and outputs of the notes are in `xi/` (the basepoint one note), `eta/` (its Dedekind eta analogue), `scale_mixture/` (normal scale mixtures and extended GGCs in the strip), `wiener/` (completeness after Wiener, Salem and Riesz), `process_map/` (the map of stochastic processes and the boundaries across basepoints), and `../sato_bondesson/` (the Sato and Bondesson note). The merged package's code folders xi_basepoint_one and eta_basepoint_one correspond to `xi/` and `eta/` here and are byte identical. The book's figures are copied to fig/bp1_*.pdf, fig/bp1eta_*.pdf and fig/sb_*.pdf.

Added for the book (not part of the authors' package): `wiener/figs_wiener.py`, which draws fig/bp1_wiener.pdf and fig/bp1_salem_riesz.pdf (the merged paper's fig_wiener.pdf and fig_salem_riesz.pdf, whose script is not in the archive) from wiener.json, salem_riesz.json, riesz_explicit.json and g_1_400.npy (run in `wiener/`: python3 figs_wiener.py, which writes fig_wiener.pdf and fig_salem_riesz.pdf there); and `wiener/riesz_explicit.py`, a reconstruction of the computation that `wiener/README.md` calls inline in the paper's description, which writes riesz_explicit.json from salem_riesz.json and g_1_400.npy, byte identical to the archived file (run in `wiener/`: python3 salem_riesz.py; python3 riesz_explicit.py). Also added: `xi/heat_primes.py`, which reads no file and writes heat_primes.json: the thresholds of Table bp:tab:heat in floating point and the ball arithmetic certificate of Proposition bp:prop:heatcert (needs python-flint; run in `xi/`: python3 heat_primes.py [workers], about 90 seconds with the default 32 worker processes).

The authors' README of the first two notes follows.

# The basepoint one: xi(1)/xi(1+sqrt s) and its Dedekind eta analogue

Nick Polson and Vadim Sokolov, September 2026.

xi_basepoint_one/    The Thorin Measure of xi(1)/xi(1+sqrt s) (18 pages, tex and pdf)
    Thorin density as Cauchy kernels at the zeros (positive unconditionally), infinite mean
    and stable one half tail 0.01303 x^(-1/2); Levy form of xi(1+s) with drift sum 1/rho and
    signed Levy measure (gamma density, prime atoms, pole) converted by the sine squared lemma;
    extension to xi(1-eps+s) equivalent to a zero free half plane; explicit error formula over
    zeros; subordination X1 = T~ + H_{T~/sqrt2}; density formula; Kent's theorem; exact Brownian
    boundary over Roberts' three zones (flat level -0.19277, tangent approximation fails at
    t = 0.091); Wald identity E b(X1) = sum 1/rho / sqrt 2.

    compute.py      Thorin density, sine squared check with primes to 10^7, no go growth
    mc.py           Monte Carlo of the subordination identity
    improve.py      Talbot density, subordination density, other basepoints
    levy.py         Levy form of xi(1+s) and its conversion
    below.py        compensated formula below the basepoint one
    errformula.py   explicit error formula over the zeros and the hybrid correction
    rs1.py          exact Brownian boundary (Volterra), rs1c.py tangent boundary,
                    rs1d.py Monte Carlo check of the boundary
    figs*.py        figures;  *.json outputs;  g_*.npy first 1700 zeta ordinates

eta_basepoint_one/   The Dedekind Eta Analogue for L(s, chi_12) (8 pages, tex and pdf)
    zeros12.py      297 zeros of L(s, chi_12) below height 340 (saved as g12.npy)
    base.py         constants, power sums, Thorin density
    arith.py        Levy form, sine squared conversion, below one, error formula
    hit.py, hit2.py subordination Monte Carlo, exact boundary, tangent boundary, checks
    figs.py         figures;  exploratory/ an earlier exploratory script

Requirements: Python 3 with numpy, scipy, matplotlib, mpmath, sympy.
Run each script from its own folder; the zero files are placed there.
Prime sieves to 10^7 need about 1 GB of memory; Monte Carlo seeds are fixed.
Compile each note with pdflatex twice.
