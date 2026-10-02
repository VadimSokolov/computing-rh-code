# peak_check.py (round 2 item kernel_fig, not part of the authors' package). A side check of the numbers that ch/ch04.tex
# quotes about the curves of fig/kernel.pdf ("for delta <= 0.02 [the maximum of |Phi(x+i eta)|] is attained at the points
# x = +-(1/2) log 2 ... at delta = 0.02, |Phi(i eta)| = 0.012 against a maximum of 3938"). The figure itself samples x on
# the grid linspace(-3, 3, 601) of code/plots3.py, whose step 0.01 misses the peak; here the maximum is located by a fine
# search, with the evaluator code/tilted.py at the precision of plots3.py (ctx.prec = 100, series cut at 2^-120), and
# the maximum on the grid of step 0.002 over |x| <= 6 that the caption of Table tab:ch4:Meta names is recomputed.
import numpy as np
from flint import acb, arb, ctx
import tilted
ctx.prec = 100
for d in [0.6854, 0.3927, 0.1, 0.02, 0.01]:
    eta = arb.pi() / 4 - arb(d)
    f = lambda x: abs(complex(tilted.Phi_c(acb(float(x), eta), 120)))
    xs = np.linspace(0, 3, 3001)
    v = np.array([f(x) for x in xs])
    k = int(v.argmax()); lo, hi = xs[max(k - 1, 0)], xs[min(k + 1, len(xs) - 1)]
    for _ in range(60):  # golden section on the bracket of the grid maximum
        a = hi - (hi - lo) * 0.6180339887498949; b = lo + (hi - lo) * 0.6180339887498949
        if f(a) > f(b): hi = b
        else: lo = a
    xm = 0.5 * (lo + hi)
    grid = np.linspace(-3, 3, 601); vg = np.array([f(x) for x in grid])
    # the grid of Table tab:ch4:Meta: step 0.002 over |x| <= 6 (|Phi| is even in x on these lines, so x >= 0 suffices)
    g2 = np.arange(0, 3001) * 0.002; v2 = np.array([f(x) for x in g2])
    print('d=%-6s eta=%.4f  |Phi(i eta)|=%.6g  max=%.6g at x=%.6f  ((1/2)log 2=%.6f)  grid max of the figure %.6g at x=%+.2f  max on the grid of step 0.002: %.6g at x=%.3f'
          % (d, np.pi / 4 - d, f(0.0), f(xm), xm, 0.5 * np.log(2), vg.max(), grid[vg.argmax()], v2.max(), g2[v2.argmax()]))
