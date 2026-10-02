#!/usr/bin/env python3
"""Interval enclosures E1 to E6 used by Theorems rz:thm:id, rz:thm:floor, rz:thm:nmono and
rz:thm:bernsteincone, with every tail bounded and every printed inequality decided by a
directed comparison of an interval endpoint with the book's threshold.

Writes enclosures_report.md beside itself; with the argument --json it also writes
enclosures_certificates.json (every check, and every enclosure to 40 digits rounded outward).

Revised in September 2026. The changes from the first version of this program:

  1. E1 and E6: each cell width, and each product of a width with an interval supremum or
     infimum, is an interval. The first version formed (hi - lo) * g(cell).b as the product
     of two mp.mpf numbers, rounded to nearest, and added it to the interval sum as an exact
     point.
  2. E3 and E5: the omitted terms j >= J of the positive series are bounded by a geometric
     series and added. The first version summed j < 400 (E3) and j < 600 (E5) and dropped
     the rest, and its comment that the E3 remainder is below e^{-810} was false: the first
     omitted term alone is about 12.2 e^{-810}.
  3. Every inequality is decided on exact rationals: the interval endpoint, converted exactly,
     against the book's threshold, read exactly from its decimal string. Upper bounds are
     printed rounded up and lower bounds rounded down. The first version printed mp.nstr,
     which rounds to nearest, and compared with 0 or 1 rather than with the book's numbers.
  4. E4 compares with the constants the proof of Theorem rz:thm:nmono uses, C_0 = 0.54993
     and r_1 = C_0 sqrt(2 pi) sqrt(2) e^{log 2 - 1.6015} < 0.786, and reports that the relaxed
     exponent 0.908 in place of 1.6015 - log 2, which the text does not use, would give more
     than 0.786 at n = 1. The first version compared a sharper quantity (the computed C_0 and
     the exponent (log 2)^2/0.3) with 1.
  5. Added: the other numerical steps of the same proofs (Trudgian's gap, the monotonicity
     conditions, the constants 2.6 and 0.27 and the large y bounds of Theorem rz:thm:nmono,
     and the uniform bound for the zeros above T_0 at large y), and two sided Darboux
     enclosures of F on a finer grid, for the remark after Theorem rz:thm:id and for Table
     rz:tab:constants.

Each line of the report names, in brackets, the statement of the book that uses it.

The enclosures, with x1 = 0.004:
  E1  F(x1) < -0.111, where F(x) = 2 theta'(0) + int_0^inf g(u, x) du,
        g(u, x) = e^{-u/4} (1 - e^{-u^2/16x}) / (1 - e^{-u}),
        2 theta'(0) = psi(1/4) - log pi = -gamma - pi/2 - 3 log 2 - log pi.
      Upper bound: g(u) <= u(1+u)/(16x) on [0, e0], e0 = mpf('1e-4') (just below 1e-4);
      upper Riemann sum of interval suprema on 12000 cells of [e0, 200], with steps 2e-4,
      1e-3, 5e-3 and 0.05 on [e0, 0.5], [0.5, 3], [3, 20] and [20, 200];
      g(u) <= e^{-u/4}/(1 - e^{-200}) on [200, inf).
  E2  P(x1) < 2e-13, P(x) = (4 pi x)^{-1/2} sum_{n>=2} Lambda(n) n^{-1/2} e^{-(log n)^2/4x}:
      n = 2 in interval arithmetic; for n >= 3, (log n)^2/4x1 >= k log n with k = log 3/(4 x1),
      so the term is at most n^{-1/2-k} log n <= n^{1-k}, and sum_{n>=3} n^{1-k} <= 2^{2-k}/(k-2).
  E3  E(x1) < 0.0011 and E(x1) < |F(x1)|/(4 sqrt(pi x1)), where
        E(x) = sum_{j>=0} a_j(x),  a_j(x) = 2 log(51+j) e^{-x((50+j)^2 - 1/4)}.
      j < J = 400 in interval arithmetic. For j >= J, log(52+j)/log(51+j) <= 1 + 1/((51+j) log(51+j)),
      so a_{j+1}/a_j <= q_J(x) = (1 + 1/((51+J) log(51+J))) e^{-x(101+2J)} < 1 and
      sum_{j>=J} a_j <= a_J/(1 - q_J).
  E4  C_0 = 0.54993 bounds e^{(log 2)^2/6y} sum_m Lambda(m) m^{-1/2} e^{-(log m)^2/6y} for
      y <= 0.05 (each term increases with y; the bound is taken at y = 0.05): m <= 40 in
      interval arithmetic, m > 40 as in E2; r_1 = 0.54993 sqrt(2 pi) sqrt(2) e^{log 2 - 1.6015} < 0.786.
  E5  E(0.003) < 0.0162 and e^{-0.003 * 14.14^2} > 0.548: as E3 with J = 600.
  E6  F(0.003) > 0.0329: lower Riemann sum of interval infima on the cells of E1; the parts
      [0, e0] and [200, inf), where g >= 0, are dropped.
"""
import json
import sys
from decimal import Decimal, localcontext, ROUND_CEILING, ROUND_FLOOR
from fractions import Fraction
from pathlib import Path

import mpmath
from mpmath import iv, mp, mpf

iv.dps = 30
mp.dps = 30

HERE = Path(__file__).resolve().parent
lines = []
checks = []
values = {}


def out(s=""):
    print(s, flush=True)
    lines.append(s)


# exact rational arithmetic on interval endpoints ---------------------------------------------
def _rational(t):
    sign, man, exp, bc = t
    if not man:
        if exp == 0:
            return Fraction(0)
        raise ValueError("interval endpoint is not finite")
    v = Fraction(int(man)) * Fraction(2) ** int(exp)
    return -v if sign else v


def lo(X):
    return _rational(iv.mpf(X)._mpi_[0])


def hi(X):
    return _rational(iv.mpf(X)._mpi_[1])


def dec(r, up, digits=12):
    """The rational r as a decimal string, rounded toward +inf (up) or toward -inf."""
    with localcontext() as c:
        c.prec = digits
        c.rounding = ROUND_CEILING if up else ROUND_FLOOR
        return str(Decimal(r.numerator) / Decimal(r.denominator))


def enc(X, digits=12):
    return "[%s, %s]" % (dec(lo(X), False, digits), dec(hi(X), True, digits))


def keep(key, X, what):
    values[key] = {"what": what, "lower": dec(lo(X), False, 40), "upper": dec(hi(X), True, 40)}


def check(key, label, X, rel, threshold, place, digits=12):
    """Certify 'label rel threshold' for the quantity enclosed by X: for '<' the upper endpoint
    of X must be below the threshold, for '>' the lower endpoint above it (exact rationals)."""
    T = Fraction(threshold)
    if rel == "<":
        ok = hi(X) < T
        shown = "%s <= %s" % (label, dec(hi(X), True, digits))
    else:
        ok = lo(X) > T
        shown = "%s >= %s" % (label, dec(lo(X), False, digits))
    status = "CERTIFIED" if ok else "FAILS"
    out("%-3s %s, so %s %s %s: %s  [%s]" % (key, shown, label, rel, threshold, status, place))
    checks.append({"key": key, "claim": "%s %s %s" % (label, rel, threshold), "place": place,
                   "bound": shown, "certified": ok})
    return ok


def check_between(key, labelA, A, labelB, B, place, digits=12):
    """Certify A < B for two enclosed quantities: upper(A) < lower(B)."""
    ok = hi(A) < lo(B)
    shown = "%s <= %s < %s <= %s" % (labelA, dec(hi(A), True, digits), dec(lo(B), False, digits), labelB)
    status = "CERTIFIED" if ok else "FAILS"
    out("%-3s %s: %s  [%s]" % (key, shown, status, place))
    checks.append({"key": key, "claim": "%s < %s" % (labelA, labelB), "place": place,
                   "bound": shown, "certified": ok})
    return ok


def check_exact(key, text, ok, place):
    status = "CERTIFIED" if ok else "FAILS"
    out("%-3s %s (exact rational arithmetic): %s  [%s]" % (key, text, status, place))
    checks.append({"key": key, "claim": text, "place": place, "bound": "exact", "certified": ok})
    return ok


# the integrand of F and Riemann sums ---------------------------------------------------------
def g(u, x):
    return iv.exp(-u / 4) * (1 - iv.exp(-u * u / (16 * x))) / (1 - iv.exp(-u))


def grid(start, segments):
    """Edges a = start; for (b, h): while a < b: append a; a += h, in mp.mpf arithmetic, as in
    the first version of this program, so that the cells of E1 and E6 are unchanged."""
    edges = []
    a = start
    for b, h in segments:
        while a < b:
            edges.append(a)
            a = a + h
    edges.append(segments[-1][0])
    return edges


def darboux(xs, edges, factored):
    """Lower and upper Riemann sums of g(., x) over the cells [e_i, e_{i+1}] for each x in xs,
    with the interval infimum and supremum of g on each cell; widths and products are
    intervals. Every edge is an exactly representable binary number."""
    low = [iv.mpf(0) for _ in xs]
    up = [iv.mpf(0) for _ in xs]
    for a, b in zip(edges[:-1], edges[1:]):
        cell = iv.mpf([a, b])
        width = iv.mpf(b) - iv.mpf(a)
        if factored:
            common = iv.exp(-cell / 4) / (1 - iv.exp(-cell))
            sq = cell * cell
        for i, x in enumerate(xs):
            if factored:
                val = common * (1 - iv.exp(-sq / (16 * x)))
            else:
                val = g(cell, x)
            low[i] += width * val.a
            up[i] += width * val.b
    return low, up


def head(e0, x):
    """int_0^{e0} u(1+u)/(16x) du, a bound for int_0^{e0} g(u, x) du."""
    return (e0 * e0 / 2 + e0 ** 3 / 3) / (16 * x)


TAIL200 = 4 * iv.exp(-iv.mpf(50)) / (1 - iv.exp(-iv.mpf(200)))   # int_200^inf e^{-u/4}/(1-e^{-200}) du

E1_GRID = [(mpf('0.5'), mpf('2e-4')), (mpf(3), mpf('1e-3')), (mpf(20), mpf('5e-3')), (mpf(200), mpf('0.05'))]
EDGES = grid(mpf('1e-4'), E1_GRID)


# the off line tail E(x) ----------------------------------------------------------------------
def a_term(j, x):
    t = iv.mpf(50 + j)
    return 2 * iv.log(t + 1) * iv.exp(-x * (t * t - iv.mpf(1) / 4))


def offline_tail(x, J):
    partial = iv.mpf(0)
    for j in range(J):
        partial += a_term(j, x)
    s = iv.mpf(51 + J)
    q = (1 + 1 / (s * iv.log(s))) * iv.exp(-x * (101 + 2 * J))
    if not hi(q) < 1:
        raise AssertionError("geometric ratio is not below 1")
    rem = a_term(J, x) / (1 - q)
    return partial, rem, q


def prime_power_base(m):
    """p if m = p^k with p prime and k >= 1, else None."""
    p = next(d for d in range(2, m + 1) if m % d == 0)
    while m % p == 0:
        m //= p
    return p if m == 1 else None


# ============================================================================================
out("mpmath %s, backend %s, interval precision %d decimal digits (prec %d bits)"
    % (mpmath.__version__, mpmath.libmp.BACKEND, iv.dps, iv.prec))
out("Every CERTIFIED line compares an interval endpoint, converted exactly to a rational, with the")
out("book's threshold as an exact rational. Upper bounds are printed rounded up, lower bounds down.")
out()

x1 = iv.mpf(1) / 250          # encloses 0.004
x6 = iv.mpf(3) / 1000         # encloses 0.003
two_thp0 = -iv.euler - iv.pi / 2 - 3 * iv.log(2) - iv.log(iv.pi)
keep("two_theta_prime_0", two_thp0, "2 theta'(0) = psi(1/4) - log pi")
c_const = -two_thp0 / 2
keep("c", c_const, "c = -theta'(0)")

# E1 and E6 on the grid of E1, with directed arithmetic ------------------------------------
e0 = iv.mpf(EDGES[0])
out("E1  grid: %d cells on [%s, 200] (first edge mpf('1e-4') = %s)"
    % (len(EDGES) - 1, mp.nstr(EDGES[0], 5), repr(EDGES[0])))
low_o, up_o = darboux([x1, x6], EDGES, factored=False)
int_up_1 = head(e0, x1) + up_o[0] + TAIL200
F1_up = two_thp0 + int_up_1
F1_low = two_thp0 + low_o[0]
F6_low = two_thp0 + low_o[1]
F6_up = two_thp0 + head(e0, x6) + up_o[1] + TAIL200
keep("E1_int_g_upper", int_up_1, "upper bound for int_0^inf g(u, 0.004) du, grid of E1")
keep("E1_F_0.004", iv.mpf([F1_low.a, F1_up.b]), "F(0.004), lower and upper Riemann sums, grid of E1")
keep("E6_F_0.003", iv.mpf([F6_low.a, F6_up.b]), "F(0.003), lower and upper Riemann sums, grid of E1")
out("E1  2 theta'(0) in %s" % enc(two_thp0, 25))
out("E1  int_0^inf g(u, 0.004) du <= %s" % dec(hi(int_up_1), True))
out("E1  F(0.004) in %s (lower and upper Riemann sums on the same cells)" % enc(iv.mpf([F1_low.a, F1_up.b])))
check("E1", "F(0.004)", F1_up, "<", "-0.111", "Theorem rz:thm:id")
absF_book = iv.mpf(111) / 1000 / (4 * iv.sqrt(iv.pi * x1))
check("E1", "0.111/(4 sqrt(pi x1))", absF_book, ">", "0.247",
      "Theorem rz:thm:id; the rounded F(x1) < -0.111 alone gives 0.247")
out()

# E2 --------------------------------------------------------------------------------------------
pref = 1 / iv.sqrt(4 * iv.pi * x1)
t2 = iv.log(2) / iv.sqrt(2) * iv.exp(-iv.log(2) ** 2 / (4 * x1))
k2 = iv.log(3) / (4 * x1)
assert lo(k2) > 2
tail2 = iv.mpf(2) ** (2 - k2) / (k2 - 2)
P_up = pref * (t2 + tail2)
P_low = pref * t2
keep("E2_P_0.004", iv.mpf([P_low.a, P_up.b]), "P(0.004): the term n = 2 alone below, with the n >= 3 bound above")
out("E2  n >= 3 contribute at most (4 pi x1)^{-1/2} * %s" % dec(hi(tail2), True, 6))
check("E2", "P(0.004)", P_up, "<", "2e-13", "Theorem rz:thm:id, the remark after it and Figure rz:fig:w")
check("E2", "(log 2)^2/(4 x1)", iv.log(2) ** 2 / (4 * x1), ">", "0.5",
      "Theorem rz:thm:id: each term of P increases in x on (0, x1]")
check("E2", "P(0.004)", P_low, ">", "1.985e-13", "Table rz:tab:constants, which prints 1.99e-13")
check("E2", "P(0.004)", P_up, "<", "1.995e-13", "Table rz:tab:constants")
out()

# E3 --------------------------------------------------------------------------------------------
part3, rem3, q3 = offline_tail(x1, 400)
E3 = part3 + rem3
keep("E3_partial_j_lt_400", part3, "sum_{j<400} a_j(0.004), the terms summed in interval arithmetic")
keep("E3_remainder_j_ge_400", rem3, "geometric bound for sum_{j>=400} a_j(0.004)")
keep("E3_q_400", q3, "ratio bound q_400(0.004)")
keep("E3_E_0.004", E3, "E(0.004) with the remainder")
a400 = a_term(400, x1)
e810 = iv.exp(-iv.mpf(810))
keep("E3_first_omitted_term", a400, "a_400(0.004) = 2 log 451 e^{-0.004(450^2 - 1/4)}")
out("E3  q_400(0.004) <= %s; remainder sum_{j>=400} a_j(0.004) <= %s"
    % (dec(hi(q3), True, 6), dec(hi(rem3), True, 6)))
check("E3", "a_400(0.004)/e^{-810}", a400 / e810, ">", "12",
      "Theorem rz:thm:id: the first omitted term (change 2 of the docstring)")
check("E3", "E(0.004)", E3, "<", "0.0011", "Theorem rz:thm:id")
rhs3 = -F1_up.b / (4 * iv.sqrt(iv.pi * x1))
check("E3", "|F(x1)|/(4 sqrt(pi x1))", rhs3, ">", "0.247", "Theorem rz:thm:id")
check_between("E3", "E(x1)", E3, "|F(x1)|/(4 sqrt(pi x1))", rhs3, "Theorem rz:thm:id, the inequality the proof uses")
check_exact("E3", "(50+j)^2 - 1/4 >= 2499.75 > 125 = 1/(2 x1) for every j >= 0, so sqrt(x) E(x) decreases on [x1, inf)",
            Fraction("2499.75") > 1 / (2 * Fraction(4, 1000)), "Theorem rz:thm:id")
# Trudgian: N(t+1) - N(t) <= main + errors <= 2 log(t+1) for t >= 50
t50 = iv.mpf(50)
M50 = iv.log((t50 + 1) / (2 * iv.pi)) / (2 * iv.pi)
R50 = 2 * (iv.mpf(112) / 1000 * iv.log(t50 + 1) + iv.mpf(278) / 1000 * iv.log(iv.log(t50 + 1))
           + iv.mpf(251) / 100 + iv.mpf(2) / 10 / t50)
keep("Trudgian_bound_at_50", M50 + R50, "(1/2pi) log(51/2pi) + 2(0.112 log 51 + 0.278 log log 51 + 2.51 + 0.2/50)")
keep("two_log_51", 2 * iv.log(t50 + 1), "2 log 51")
out("E3  at t = 50: main term increment plus errors in %s, 2 log 51 in %s"
    % (enc(M50 + R50, 8), enc(2 * iv.log(t50 + 1), 8)))
check("E3", "2 log 51 - (bound for N(51) - N(50))", 2 * iv.log(t50 + 1) - M50 - R50, ">", "0.86",
      "Theorem rz:thm:id ('about 7.0 at t=50 against 2 log 51 = 7.86')")
check("E3", "2 - 1/(2 pi) - 0.224 - 0.556/log 51", 2 - 1 / (2 * iv.pi) - iv.mpf(224) / 1000
      - iv.mpf(556) / 1000 / iv.log(t50 + 1), ">", "0",
      "Theorem rz:thm:id ('the gap widens'): the derivative of the gap is this over t+1 plus 0.4/t^2")
check("E3", "the bound at t = 50", M50 + R50, "<", "7.05", "Theorem rz:thm:id ('about 7.0')")
check("E3", "the bound at t = 50", M50 + R50, ">", "6.95", "Theorem rz:thm:id ('about 7.0')")
check("E3", "2 log 51", 2 * iv.log(t50 + 1), ">", "7.855", "Theorem rz:thm:id ('7.86')")
check("E3", "2 log 51", 2 * iv.log(t50 + 1), "<", "7.865", "Theorem rz:thm:id ('7.86')")
out()

# E4 and the constants of Theorem rz:thm:nmono ----------------------------------------------------
y4 = iv.mpf(1) / 20
S4 = iv.mpf(0)
for m in range(2, 41):
    p = prime_power_base(m)
    if p is not None:
        S4 += iv.log(p) / iv.sqrt(m) * iv.exp(-(iv.log(m) ** 2 - iv.log(2) ** 2) / (6 * y4))
k4 = iv.log(41) / (6 * y4)
assert lo(k4) > 2
tail4 = iv.exp(iv.log(2) ** 2 / (6 * y4)) * iv.mpf(40) ** (2 - k4) / (k4 - 2)
C0enc = S4 + tail4
keep("E4_C0_sum_at_0.05", C0enc, "e^{(log 2)^2/0.3} sum_m Lambda(m) m^{-1/2} e^{-(log m)^2/0.3}")
out("E4  terms m > 40 contribute at most %s" % dec(hi(tail4), True, 6))
check("E4", "C0 sum at y = 0.05", C0enc, "<", "0.54993", "Theorem rz:thm:nmono ('C_0=0.54993 bounds ...')")
C0 = iv.mpf(54993) / 100000
sqrt2pi = iv.sqrt(2 * iv.pi)
r1 = C0 * sqrt2pi * iv.sqrt(2) * iv.exp(iv.log(2) - iv.mpf(16015) / 10000)
keep("r1_book", r1, "r_1 = 0.54993 sqrt(2 pi) sqrt(2) e^{log 2 - 1.6015}, as the book defines it")
check("E4", "(log 2)^2/(6 * 0.05)", iv.log(2) ** 2 / (6 * y4), ">", "1.6015",
      "Theorem rz:thm:nmono ('(log2)^2/6y >= 1.6015 n' for y <= y_n)")
check("E4", "r_1 = C_0 sqrt(2 pi) sqrt 2 e^{log 2 - 1.6015}", r1, "<", "0.786", "Theorem rz:thm:nmono ('r_n <= r_1 < 0.786')")
check("E4", "log 2 - 1.6015", iv.log(2) - iv.mpf(16015) / 10000, "<", "-0.908",
      "Theorem rz:thm:nmono: a relaxed exponent, not used in the text")
r1rel = C0 * sqrt2pi * iv.sqrt(2) * iv.exp(-iv.mpf(908) / 1000)
keep("r1_relaxed", r1rel, "C_0 sqrt(2 pi) sqrt(2) e^{-0.908}, the relaxed exponent at n = 1")
check("E4", "C_0 sqrt(2 pi) sqrt 2 e^{-0.908}", r1rel, ">", "0.786",
      "Theorem rz:thm:nmono: the relaxed exponent 0.908 would not give r_1 < 0.786")
astar = iv.log(C0 * sqrt2pi * iv.sqrt(2) / (iv.mpf(786) / 1000))
keep("relaxed_exponent_needed", astar, "the least a with C_0 sqrt(2 pi) sqrt(2) e^{-a} <= 0.786")
check("E4", "least usable relaxed exponent a*", astar, "<", "0.90835", "Theorem rz:thm:nmono: the relaxed exponent 0.90835 would give r_1 < 0.786")
check("E4", "1.6015 - log 2", iv.mpf(16015) / 10000 - iv.log(2), ">", "0.90835",
      "Theorem rz:thm:nmono: e^{(log 2 - 1.6015) n} <= e^{-0.90835 n}")
check("E4", "r_{n+1}/r_n <= sqrt(3/2) e^{log 2 - 1.6015}", iv.sqrt(iv.mpf(3) / 2) * iv.exp(iv.log(2) - iv.mpf(16015) / 10000),
      "<", "0.5", "Theorem rz:thm:nmono ('which decreases in n')")
R1prog = C0enc * sqrt2pi * iv.sqrt(2) * iv.exp(iv.log(2) - iv.log(2) ** 2 / (iv.mpf(3) / 10))
keep("R1_sharper", R1prog, "r_1 with the enclosed C_0 sum and the exponent log 2 - (log 2)^2/0.3")
out("E4  with the enclosed C_0 sum and the exponent log 2 - (log 2)^2/0.3 in place of 0.54993 and log 2 - 1.6015, r_1 lies in %s"
    % enc(R1prog, 10))
G32 = iv.sqrt(iv.pi) / 2
first_const = iv.mpf(214) / 1000 * G32 * iv.mpf(20) ** (iv.mpf(3) / 2) / (2 * iv.pi)
check("E4", "0.214 Gamma(3/2) 20^{3/2}/(2 pi)", first_const, ">", "2.6", "Theorem rz:thm:nmono ('> 2.6')")
first_rel = (1 - r1rel) * G32 * iv.mpf(20) ** (iv.mpf(3) / 2) / (2 * iv.pi)
check("E4", "(1 - C_0 sqrt(2 pi) sqrt 2 e^{-0.908}) Gamma(3/2) 20^{3/2}/(2 pi)", first_rel, ">", "2.6",
      "Theorem rz:thm:nmono: 2.6 would hold even with the relaxed exponent 0.908")
negs = iv.mpf(1) / (8 * 3 * iv.pi) + iv.exp(iv.mpf(5) / 100 / 4) / 4
check("E4", "2^{-3}/(3 pi) + 4^{-1} e^{0.05/4} (n = 1, the largest case)", negs, "<", "0.27", "Theorem rz:thm:nmono ('< 0.27')")
T0 = iv.mpf(3 * 10 ** 12)
T0q = Fraction(3 * 10 ** 12)
NN = 10 ** 10
check_exact("E4", "y_n T_0^2 >= (0.05/10^10)(3 10^12)^2 = 4.5e13 for n <= 10^10",
            Fraction(5, 100) / NN * T0q ** 2 >= Fraction("4.5e13"), "Theorem rz:thm:nmono, large y")
check("E4", "2 * 10^10 * log(T_0 + 1)", 2 * NN * iv.log(T0 + 1), "<", "5.8e11", "Theorem rz:thm:nmono, large y")
check_exact("E4", "n/y <= 20 n^2 <= 2e21 < T_0^2 for n <= 10^10", 20 * Fraction(NN) ** 2 < T0q ** 2, "Theorem rz:thm:nmono, large y")
# large y: the zeros above T_0 contribute at most 4 log(T_0+1) (T_0+1)^{2n} e^{-A y}, A = T_0^2 - 1/4
A = T0 * T0 - iv.mpf(1) / 4
# y_n A > 2 n log(T_0+1) + 6 with y_n = 1/(20 n) means A > 40 log(T_0+1) n^2 + 120 n; n* is the positive root
lT1 = iv.log(T0 + 1)
nstar = (iv.sqrt(120 ** 2 + 160 * lT1 * A) - 120) / (80 * lT1)
keep("n_star", nstar, "largest n with y_n (T_0^2 - 1/4) > 2 n log(T_0+1) + 6")
check("E4", "n* with y_n (T_0^2 - 1/4) = 2 n log(T_0+1) + 6", nstar, ">", "8.8e10",
      "Theorem rz:thm:nmono ('This sufficient condition holds while n < 8.8e10')")
first_term = iv.log(8 * iv.pi * lT1 / (c_const * G32))
check("E4", "log(8 pi log(T_0+1)/(c Gamma(3/2)))", first_term, "<", "6",
      "Theorem rz:thm:nmono, large y ('its first term is below 6')")
qmax = iv.exp(2 * NN / T0 + 1 / (T0 * iv.log(T0)) - (2 * T0 + 1) / (20 * NN))
keep("largey_qmax", qmax, "ratio bound b_{j+1}/b_j for n <= 10^10, y >= 1/(20 n)")
check("E4", "b_{j+1}/b_j for n <= 10^10, y >= y_n", qmax, "<", "1e-13", "Theorem rz:thm:nmono, large y: consecutive blocks of zeros above T_0")
check_exact("E4", "20 N (N + 1/2) < T_0^2 - 1/4 for N = 10^10 (tail to floor ratio decreases in y >= y_n)",
            20 * Fraction(NN) * (Fraction(NN) + Fraction(1, 2)) < T0q ** 2 - Fraction(1, 4), "Theorem rz:thm:nmono, large y")
Lrat = (iv.log(8 * iv.pi * iv.log(T0 + 1) / (c_const * G32)) + 2 * NN * iv.log(T0 + 1) - A / (20 * NN))
keep("largey_log_ratio", Lrat, "upper bound for log(tail/floor term), n <= 10^10, y >= y_n")
check("E4", "log(tail / (c Gamma(n+1/2) y^{-n-1/2}/(2 pi)))", Lrat, "<", "-4.4e13", "Theorem rz:thm:nmono, large y: the tail is below the positive term")
n88 = 88 * 10 ** 9
q88 = iv.exp(2 * iv.mpf(n88) / T0 + 1 / (T0 * iv.log(T0)) - (2 * T0 + 1) / (20 * iv.mpf(n88)))
check("E4", "b_{j+1}/b_j at n = 8.8e10", q88, "<", "0.5", "Theorem rz:thm:nmono ('while n < 8.8e10'), large y")
L88 = iv.log(8 * iv.pi * iv.log(T0 + 1) / (c_const * G32)) + 2 * n88 * iv.log(T0 + 1) - A / (20 * iv.mpf(n88))
check("E4", "log ratio bound at n = 8.8e10", L88, "<", "0", "Theorem rz:thm:nmono ('while n < 8.8e10'), large y")
check_exact("E4", "20 n^2 < T_0^2 at n = 8.8e10", 20 * Fraction(n88) ** 2 < T0q ** 2, "Theorem rz:thm:nmono ('while n < 8.8e10')")
out()

# E5 --------------------------------------------------------------------------------------------
y5 = iv.mpf(3) / 1000
part5, rem5, q5 = offline_tail(y5, 600)
E5 = part5 + rem5
first = iv.exp(-y5 * (iv.mpf(1414) / 100) ** 2)
keep("E5_partial_j_lt_600", part5, "sum_{j<600} a_j(0.003), the terms summed in interval arithmetic")
keep("E5_remainder_j_ge_600", rem5, "geometric bound for sum_{j>=600} a_j(0.003)")
keep("E5_E_0.003", E5, "E(0.003) with the remainder")
keep("E5_first_zero_term", first, "e^{-0.003 * 14.14^2}")
out("E5  q_600(0.003) <= %s; remainder sum_{j>=600} a_j(0.003) <= %s" % (dec(hi(q5), True, 6), dec(hi(rem5), True, 6)))
check("E5", "E(0.003)", E5, "<", "0.0162", "Theorem rz:thm:bernsteincone")
check("E5", "e^{-0.003 * 14.14^2}", first, ">", "0.548", "Theorem rz:thm:bernsteincone")
check_between("E5", "E(0.003)", E5, "e^{-0.003 * 14.14^2}", first, "Theorem rz:thm:bernsteincone, the inequality the proof uses")
check_exact("E5", "(50+j)^2 - 1/4 - 14.14^2 >= 2299.8104 > 0 for j >= 0, so E(y) e^{y 14.14^2} decreases in y",
            Fraction("2499.75") - Fraction("14.14") ** 2 > 0, "Theorem rz:thm:bernsteincone")
out()

# E6 --------------------------------------------------------------------------------------------
out("E6  F(0.003) in %s (lower and upper Riemann sums on the grid of E1)" % enc(iv.mpf([F6_low.a, F6_up.b])))
check("E6", "F(0.003)", F6_low, ">", "0.0329", "Theorem rz:thm:bernsteincone")
out()

# a finer grid: the remark after Theorem rz:thm:id and Table rz:tab:constants -------------------
FINE = grid(mpf(2) ** -15, [(mpf(1) / 2, mpf(2) ** -15), (mpf(3), mpf(2) ** -12), (mpf(20), mpf(2) ** -11),
                            (mpf(200), mpf(2) ** -7)])
xsF = [iv.mpf(1) / 250, iv.mpf(35) / 10000, iv.mpf(3) / 1000]
namesF = ["0.004", "0.0035", "0.003"]
tableF = ["-0.1145624", "-0.0442285", "0.0366850"]
out("F   finer grid: %d cells, steps 2^-15, 2^-12, 2^-11, 2^-7 on [2^-15, 1/2], [1/2, 3], [3, 20], [20, 200]"
    % (len(FINE) - 1))
lowF, upF = darboux(xsF, FINE, factored=True)
e0F = iv.mpf(FINE[0])
encF = {}
for i, (x, nm, tv) in enumerate(zip(xsF, namesF, tableF)):
    Fl = two_thp0 + lowF[i]
    Fu = two_thp0 + head(e0F, x) + upF[i] + TAIL200
    encF[nm] = (Fl, Fu)
    keep("F_fine_" + nm, iv.mpf([Fl.a, Fu.b]), "F(%s), finer grid" % nm)
    out("F   F(%s) in %s" % (nm, enc(iv.mpf([Fl.a, Fu.b]), 10)))
    tq = Fraction(tv)
    check_exact("F", "Table value F(%s) = %s lies in the enclosure" % (nm, tv), lo(Fl) <= tq <= hi(Fu), "Table rz:tab:constants")
Fu4 = encF["0.004"][1]
check("F", "-W_0(x1) = |F(x1)|/(4 sqrt(pi x1))", -Fu4.b / (4 * iv.sqrt(iv.pi * x1)), ">", "0.25",
      "the remark after Theorem rz:thm:id ('-W_0(x_1) already exceeds 0.25')")
out("F   with E1's grid alone, -W_0(x1) >= %s only, which does not reach 0.25"
    % dec(lo(rhs3), False, 8))
out()

failed = [c for c in checks if not c["certified"]]
out("%d checks, %d certified, %d failed" % (len(checks), len(checks) - len(failed), len(failed)))

report = HERE / "enclosures_report.md"
with open(report, "w") as fh:
    fh.write("# Interval enclosures E1 to E6\n\n")
    fh.write("Output of `enclosures.py` (mpmath.iv, outward rounding; "
             "tails bounded; directed comparisons on exact rationals).\n\n```\n")
    fh.write("\n".join(lines) + "\n```\n")
if "--json" in sys.argv[1:]:
    with open(HERE / "enclosures_certificates.json", "w") as fh:
        json.dump({"mpmath": mpmath.__version__, "iv_dps": iv.dps, "checks": checks, "values": values},
                  fh, indent=1)
if failed:
    raise SystemExit("some checks failed")
