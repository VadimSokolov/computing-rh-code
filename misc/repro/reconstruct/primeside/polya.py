"""RECONSTRUCTION, NOT THE AUTHORS' CODE.

Compatibility shim for the module `polya` that code/primeside.py imports
(`from polya import xi_pair`, line 3) and that is missing from the archive.
It exists only so that the archived code/primeside.py can run unchanged and
regenerate data/results_prime.json. Written for misc/repro/reconstruct/primeside/.

Interface, as primeside.py uses it (lines 28 and 29):
    F, F1 = xi_pair(sigma, b)      sigma = alpha - 1/2 (real), b = numpy array of ordinates
    pol = (F1[0] / F[0]).real      = Re (xi'/xi)(alpha + i b[0])
so F[k] = xi(1/2 + sigma + i b[k]) and F1[k] = xi'(1/2 + sigma + i b[k]), returned here as
numpy complex128 arrays of the same shape as b.

The values come from the archived evaluator code/polya_flint.py (trapezoid rule for the
Polya integral, ball arithmetic in python-flint) with setup(420, 200, 3): 420 bits,
h = 1/200, |u| <= 3. These are the parameters of Section sec:ch4:trapezoid, used by part1.py
and part2.py, and the ones code/plots.py uses when it recomputes the same Polya column of
results_prime.json for data/extra.json. The node values are computed once, on first use.
"""
import numpy as np
from flint import acb, ctx
from polya_flint import setup, xi_pair as _xi_pair_flint

PREC_BITS, H_INV, UMAX = 420, 200, 3
_S = None


def _nodes():
    global _S
    if _S is None:
        _S = setup(PREC_BITS, H_INV, UMAX)
    return _S


def xi_pair(sigma, b):
    S = _nodes()
    ctx.prec = PREC_BITS
    b = np.atleast_1d(np.asarray(b, dtype=float))
    F = np.empty(b.shape, dtype=complex)
    F1 = np.empty(b.shape, dtype=complex)
    for k, bk in np.ndenumerate(b):
        A, B = _xi_pair_flint(S, acb(float(sigma), float(bk)))
        F[k] = complex(A)
        F1[k] = complex(B)
    return F, F1
