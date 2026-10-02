# Interval enclosures E1 to E6

Output of `enclosures.py` (mpmath.iv, outward rounding; tails bounded; directed comparisons on exact rationals).

```
mpmath 1.3.0, backend python, interval precision 30 decimal digits (prec 103 bits)
Every CERTIFIED line compares an interval endpoint, converted exactly to a rational, with the
book's threshold as an exact rational. Upper bounds are printed rounded up, lower bounds down.

E1  grid: 12000 cells on [0.0001, 200] (first edge mpf('1e-4') = mpf('0.0000999999999999999999999999999999982'))
E1  2 theta'(0) in [-5.372183419225665582232958, -5.372183419225665582232957]
E1  int_0^inf g(u, 0.004) du <= 5.26113527223
E1  F(0.004) in [-0.118065274462, -0.111048147003] (lower and upper Riemann sums on the same cells)
E1  F(0.004) <= -0.111048147003, so F(0.004) < -0.111: CERTIFIED  [Theorem rz:thm:id]
E1  0.111/(4 sqrt(pi x1)) >= 0.247547221116, so 0.111/(4 sqrt(pi x1)) > 0.247: CERTIFIED  [Theorem rz:thm:id; the rounded F(x1) < -0.111 alone gives 0.247]

E2  n >= 3 contribute at most (4 pi x1)^{-1/2} * 1.28372E-22
E2  P(0.004) <= 1.98858765667E-13, so P(0.004) < 2e-13: CERTIFIED  [Theorem rz:thm:id, the remark after it and Figure rz:fig:w]
E2  (log 2)^2/(4 x1) >= 30.0283133698, so (log 2)^2/(4 x1) > 0.5: CERTIFIED  [Theorem rz:thm:id: each term of P increases in x on (0, x1]]
E2  P(0.004) >= 1.98858765094E-13, so P(0.004) > 1.985e-13: CERTIFIED  [Table rz:tab:constants, which prints 1.99e-13]
E2  P(0.004) <= 1.98858765667E-13, so P(0.004) < 1.995e-13: CERTIFIED  [Table rz:tab:constants]

E3  q_400(0.004) <= 0.0272246; remainder sum_{j>=400} a_j(0.004) <= 2.09444E-351
E3  a_400(0.004)/e^{-810} >= 12.2351637271, so a_400(0.004)/e^{-810} > 12: CERTIFIED  [Theorem rz:thm:id: the first omitted term (change 2 of the docstring)]
E3  E(0.004) <= 0.00105399390764, so E(0.004) < 0.0011: CERTIFIED  [Theorem rz:thm:id]
E3  |F(x1)|/(4 sqrt(pi x1)) >= 0.247654596403, so |F(x1)|/(4 sqrt(pi x1)) > 0.247: CERTIFIED  [Theorem rz:thm:id]
E3  E(x1) <= 0.00105399390764 < 0.247654596403 <= |F(x1)|/(4 sqrt(pi x1)): CERTIFIED  [Theorem rz:thm:id, the inequality the proof uses]
E3  (50+j)^2 - 1/4 >= 2499.75 > 125 = 1/(2 x1) for every j >= 0, so sqrt(x) E(x) decreases on [x1, inf) (exact rational arithmetic): CERTIFIED  [Theorem rz:thm:id]
E3  at t = 50: main term increment plus errors in [7.0032129, 7.0032130], 2 log 51 in [7.8636512, 7.8636513]
E3  2 log 51 - (bound for N(51) - N(50)) >= 0.860438315092, so 2 log 51 - (bound for N(51) - N(50)) > 0.86: CERTIFIED  [Theorem rz:thm:id ('about 7.0 at t=50 against 2 log 51 = 7.86')]
E3  2 - 1/(2 pi) - 0.224 - 0.556/log 51 >= 1.47543492026, so 2 - 1/(2 pi) - 0.224 - 0.556/log 51 > 0: CERTIFIED  [Theorem rz:thm:id ('the gap widens'): the derivative of the gap is this over t+1 plus 0.4/t^2]
E3  the bound at t = 50 <= 7.00321295036, so the bound at t = 50 < 7.05: CERTIFIED  [Theorem rz:thm:id ('about 7.0')]
E3  the bound at t = 50 >= 7.00321295035, so the bound at t = 50 > 6.95: CERTIFIED  [Theorem rz:thm:id ('about 7.0')]
E3  2 log 51 >= 7.86365126544, so 2 log 51 > 7.855: CERTIFIED  [Theorem rz:thm:id ('7.86')]
E3  2 log 51 <= 7.86365126545, so 2 log 51 < 7.865: CERTIFIED  [Theorem rz:thm:id ('7.86')]

E4  terms m > 40 contribute at most 1.12796E-17
E4  C0 sum at y = 0.05 <= 0.549924871856, so C0 sum at y = 0.05 < 0.54993: CERTIFIED  [Theorem rz:thm:nmono ('C_0=0.54993 bounds ...')]
E4  (log 2)^2/(6 * 0.05) >= 1.60151004639, so (log 2)^2/(6 * 0.05) > 1.6015: CERTIFIED  [Theorem rz:thm:nmono ('(log2)^2/6y >= 1.6015 n' for y <= y_n)]
E4  r_1 = C_0 sqrt(2 pi) sqrt 2 e^{log 2 - 1.6015} <= 0.785994898103, so r_1 = C_0 sqrt(2 pi) sqrt 2 e^{log 2 - 1.6015} < 0.786: CERTIFIED  [Theorem rz:thm:nmono ('r_n <= r_1 < 0.786')]
E4  log 2 - 1.6015 <= -0.908352819440, so log 2 - 1.6015 < -0.908: CERTIFIED  [Theorem rz:thm:nmono: a relaxed exponent, not used in the text]
E4  C_0 sqrt(2 pi) sqrt 2 e^{-0.908} >= 0.786272261308, so C_0 sqrt(2 pi) sqrt 2 e^{-0.908} > 0.786: CERTIFIED  [Theorem rz:thm:nmono: the relaxed exponent 0.908 would not give r_1 < 0.786]
E4  least usable relaxed exponent a* <= 0.908346328455, so least usable relaxed exponent a* < 0.90835: CERTIFIED  [Theorem rz:thm:nmono: the relaxed exponent 0.90835 would give r_1 < 0.786]
E4  1.6015 - log 2 >= 0.908352819440, so 1.6015 - log 2 > 0.90835: CERTIFIED  [Theorem rz:thm:nmono: e^{(log 2 - 1.6015) n} <= e^{-0.90835 n}]
E4  r_{n+1}/r_n <= sqrt(3/2) e^{log 2 - 1.6015} <= 0.493802190840, so r_{n+1}/r_n <= sqrt(3/2) e^{log 2 - 1.6015} < 0.5: CERTIFIED  [Theorem rz:thm:nmono ('which decreases in n')]
E4  with the enclosed C_0 sum and the exponent log 2 - (log 2)^2/0.3 in place of 0.54993 and log 2 - 1.6015, r_1 lies in [0.7859796723, 0.7859796724]
E4  0.214 Gamma(3/2) 20^{3/2}/(2 pi) >= 2.69975179856, so 0.214 Gamma(3/2) 20^{3/2}/(2 pi) > 2.6: CERTIFIED  [Theorem rz:thm:nmono ('> 2.6')]
E4  (1 - C_0 sqrt(2 pi) sqrt 2 e^{-0.908}) Gamma(3/2) 20^{3/2}/(2 pi) >= 2.69631704174, so (1 - C_0 sqrt(2 pi) sqrt 2 e^{-0.908}) Gamma(3/2) 20^{3/2}/(2 pi) > 2.6: CERTIFIED  [Theorem rz:thm:nmono: 2.6 would hold even with the relaxed exponent 0.908]
E4  2^{-3}/(3 pi) + 4^{-1} e^{0.05/4} (n = 1, the largest case) <= 0.266407524810, so 2^{-3}/(3 pi) + 4^{-1} e^{0.05/4} (n = 1, the largest case) < 0.27: CERTIFIED  [Theorem rz:thm:nmono ('< 0.27')]
E4  y_n T_0^2 >= (0.05/10^10)(3 10^12)^2 = 4.5e13 for n <= 10^10 (exact rational arithmetic): CERTIFIED  [Theorem rz:thm:nmono, large y]
E4  2 * 10^10 * log(T_0 + 1) <= 574592668092, so 2 * 10^10 * log(T_0 + 1) < 5.8e11: CERTIFIED  [Theorem rz:thm:nmono, large y]
E4  n/y <= 20 n^2 <= 2e21 < T_0^2 for n <= 10^10 (exact rational arithmetic): CERTIFIED  [Theorem rz:thm:nmono, large y]
E4  n* with y_n (T_0^2 - 1/4) = 2 n log(T_0+1) + 6 >= 88496524862.1, so n* with y_n (T_0^2 - 1/4) = 2 n log(T_0+1) + 6 > 8.8e10: CERTIFIED  [Theorem rz:thm:nmono ('This sufficient condition holds while n < 8.8e10')]
E4  log(8 pi log(T_0+1)/(c Gamma(3/2))) <= 5.71479553698, so log(8 pi log(T_0+1)/(c Gamma(3/2))) < 6: CERTIFIED  [Theorem rz:thm:nmono, large y ('its first term is below 6')]
E4  b_{j+1}/b_j for n <= 10^10, y >= y_n <= 9.42021553198E-14, so b_{j+1}/b_j for n <= 10^10, y >= y_n < 1e-13: CERTIFIED  [Theorem rz:thm:nmono, large y: consecutive blocks of zeros above T_0]
E4  20 N (N + 1/2) < T_0^2 - 1/4 for N = 10^10 (tail to floor ratio decreases in y >= y_n) (exact rational arithmetic): CERTIFIED  [Theorem rz:thm:nmono, large y]
E4  log(tail / (c Gamma(n+1/2) y^{-n-1/2}/(2 pi))) <= -4.44254073319E+13, so log(tail / (c Gamma(n+1/2) y^{-n-1/2}/(2 pi))) < -4.4e13: CERTIFIED  [Theorem rz:thm:nmono, large y: the tail is below the positive term]
E4  b_{j+1}/b_j at n = 8.8e10 <= 0.0350694729863, so b_{j+1}/b_j at n = 8.8e10 < 0.5: CERTIFIED  [Theorem rz:thm:nmono ('while n < 8.8e10'), large y]
E4  log ratio bound at n = 8.8e10 <= -57220884421.5, so log ratio bound at n = 8.8e10 < 0: CERTIFIED  [Theorem rz:thm:nmono ('while n < 8.8e10'), large y]
E4  20 n^2 < T_0^2 at n = 8.8e10 (exact rational arithmetic): CERTIFIED  [Theorem rz:thm:nmono ('while n < 8.8e10')]

E5  q_600(0.003) <= 0.0201861; remainder sum_{j>=600} a_j(0.003) <= 4.50227E-550
E5  E(0.003) <= 0.0161664660484, so E(0.003) < 0.0162: CERTIFIED  [Theorem rz:thm:bernsteincone]
E5  e^{-0.003 * 14.14^2} >= 0.548911089772, so e^{-0.003 * 14.14^2} > 0.548: CERTIFIED  [Theorem rz:thm:bernsteincone]
E5  E(0.003) <= 0.0161664660484 < 0.548911089772 <= e^{-0.003 * 14.14^2}: CERTIFIED  [Theorem rz:thm:bernsteincone, the inequality the proof uses]
E5  (50+j)^2 - 1/4 - 14.14^2 >= 2299.8104 > 0 for j >= 0, so E(y) e^{y 14.14^2} decreases in y (exact rational arithmetic): CERTIFIED  [Theorem rz:thm:bernsteincone]

E6  F(0.003) in [0.0329774347638, 0.0404070606964] (lower and upper Riemann sums on the grid of E1)
E6  F(0.003) >= 0.0329774347638, so F(0.003) > 0.0329: CERTIFIED  [Theorem rz:thm:bernsteincone]

F   finer grid: 84479 cells, steps 2^-15, 2^-12, 2^-11, 2^-7 on [2^-15, 1/2], [1/2, 3], [3, 20], [20, 200]
F   F(0.004) in [-0.1151125453, -0.1140120396]
F   Table value F(0.004) = -0.1145624 lies in the enclosure (exact rational arithmetic): CERTIFIED  [Table rz:tab:constants]
F   F(0.0035) in [-0.04479168702, -0.04366495970]
F   Table value F(0.0035) = -0.0442285 lies in the enclosure (exact rational arithmetic): CERTIFIED  [Table rz:tab:constants]
F   F(0.003) in [0.03610472483, 0.03726569364]
F   Table value F(0.003) = 0.0366850 lies in the enclosure (exact rational arithmetic): CERTIFIED  [Table rz:tab:constants]
F   -W_0(x1) = |F(x1)|/(4 sqrt(pi x1)) >= 0.254264536878, so -W_0(x1) = |F(x1)|/(4 sqrt(pi x1)) > 0.25: CERTIFIED  [the remark after Theorem rz:thm:id ('-W_0(x_1) already exceeds 0.25')]
F   with E1's grid alone, -W_0(x1) >= 0.24765459 only, which does not reach 0.25

48 checks, 48 certified, 0 failed
```
