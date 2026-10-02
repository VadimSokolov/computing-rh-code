A map of stochastic processes for the Riemann Xi function (8 pages): process_map_xi.tex / .pdf

comp.py      Biane Pitman Yor Bessel identity, first Wiener gamma run, subordinated BM at alpha = 1
comp2.py     Riemann OU (short run), Wiener gamma with exact cell averages, subordinated BM across alpha
ou.py        long exact simulation of the Riemann Ornstein Uhlenbeck process (6.8 million jumps)
bnd_alpha.py exact Brownian boundary at a basepoint: python3 bnd_alpha.py 0.75  (writes bnd_0.75.json)
figs.py      sample path figure fig_paths.pdf
Data: g_*.npy (1700 zeta ordinates), alpha1.json and rs1.json (basepoint one note), volterra.json (centre boundary).
Requirements: numpy, scipy, matplotlib, mpmath.
