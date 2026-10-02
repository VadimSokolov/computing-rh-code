"""Common definitions for the torus functions of Theorem rz:thm:torus (ch/rz4.tex).

beta_k = k^2 prod_{j<=n, j!=k} (1 - k^2/j^2)^{-2},  alpha_k = -2 beta_k sum_{j<=n, j!=k} k^2/(j^2 - k^2),
P_n(x, w) = sum_k beta_k k^{-2x} exp(-i sum_p v_p(k) w_p),  P_n'(x, w) = -sum_k 2 log k beta_k k^{-2x} exp(...),
r(x) = sum_k |alpha_k| k^{-2x},  r'(x) = sum_k 2 log k |alpha_k| k^{-2x}.
"""
from fractions import Fraction as Fr
import math
import numpy as np


def exact_coeffs(n):
    """Exact rationals (alpha_k, beta_k), k = 1..n."""
    out = []
    for k in range(1, n + 1):
        prod = Fr(1)
        s = Fr(0)
        for j in range(1, n + 1):
            if j != k:
                prod *= 1 - Fr(k * k, j * j)
                s += Fr(k * k, j * j - k * k)
        beta = Fr(k * k) / (prod * prod)
        out.append((-2 * beta * s, beta))
    return out


def primes_upto(n):
    return [p for p in range(2, n + 1) if all(p % q for q in range(2, int(p ** 0.5) + 1))]


def vp(k, p):
    e = 0
    while k % p == 0:
        k //= p
        e += 1
    return e


def exponent_matrix(n):
    ps = primes_upto(n)
    V = np.array([[vp(k, p) for p in ps] for k in range(1, n + 1)], dtype=float)
    return ps, V


class Torus:
    def __init__(self, n):
        self.n = n
        ab = exact_coeffs(n)
        self.alpha = np.array([float(a) for a, b in ab])
        self.beta = np.array([float(b) for a, b in ab])
        self.k = np.arange(1, n + 1, dtype=float)
        self.logk = np.log(self.k)
        self.ps, self.V = exponent_matrix(n)
        self.K = len(self.ps)

    def c(self, x):
        x = np.asarray(x, dtype=float)
        return self.beta * np.exp(-2.0 * np.multiply.outer(x, self.logk))

    def r(self, x):
        return float(np.sum(np.abs(self.alpha) * np.exp(-2 * x * self.logk)))

    def rp(self, x):
        return float(np.sum(2 * self.logk * np.abs(self.alpha) * np.exp(-2 * x * self.logk)))

    def PP(self, x, w):
        """P and P' at arrays x (N,), w (N, K)."""
        c = self.c(x)
        e = np.exp(-1j * (w @ self.V.T))
        P = np.sum(c * e, axis=-1)
        P1 = np.sum((-2 * self.logk) * c * e, axis=-1)
        return P, P1

    def PP_grad(self, x, w):
        """P, P' and their gradients in (x, w) at a single point."""
        c = self.beta * np.exp(-2.0 * x * self.logk)
        e = np.exp(-1j * (self.V @ w))
        ce = c * e
        P = ce.sum()
        P1 = (-2 * self.logk * ce).sum()
        dP = np.empty(1 + self.K, dtype=complex)
        dP1 = np.empty(1 + self.K, dtype=complex)
        dP[0] = P1
        dP1[0] = (4 * self.logk ** 2 * ce).sum()
        dP[1:] = (-1j) * (self.V.T @ ce)
        dP1[1:] = (-1j) * (self.V.T @ (-2 * self.logk * ce))
        return P, P1, dP, dP1
