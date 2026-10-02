"""RECONSTRUCTION AID, NOT THE AUTHORS' CODE.

Checks the 'gamma' column of results_prime.json (code/primeside.py computes the archimedean
integral int_0^60 (1 - cos bx) e^{-alpha x}/(e^{2x} - 1) dx with scipy.integrate.quad in
double precision) against the closed form
    int_0^oo (1 - cos bx) e^{-alpha x}/(e^{2x} - 1) dx = (1/2) Re[psi(1 + (alpha + i b)/2) - psi(1 + alpha/2)],
obtained by expanding 1/(e^{2x} - 1) = sum_{k>=1} e^{-2kx}, evaluated in mpmath at 40 digits.
It then recomputes the bracket c + gamma + pole + primes with the exact archimedean piece and
reports bracket - exact, the error of the prime side in Table tab:prime (ch02).
Usage: python archimedean_check.py BOOK_results_prime.json RERUN_results_prime.json
"""
import json, sys
import mpmath as mp

mp.mp.dps = 40
book = json.load(open(sys.argv[1]))
rerun = json.load(open(sys.argv[2]))
out = []
print(f"{'alpha':>5s} {'b':>10s} {'gamma closed form':>22s} {'scipy - closed (book)':>22s} {'scipy - closed (rerun)':>22s} "
      f"{'book bracket - exact':>21s} {'bracket_cf - exact':>19s}")
for rb, rr in zip(book, rerun):
    a, b = mp.mpf(rb['alpha']), mp.mpf(rb['b'])
    g = (mp.re(mp.digamma(1 + (a + 1j * b) / 2)) - mp.digamma(1 + a / 2)) / 2
    # tail of the integral beyond x = 60, which primeside.py cuts off: below e^{-(alpha+2)60}
    br_cf = mp.mpf(rb['c']) + g + mp.mpf(rb['pole']) + mp.mpf(rb['primes'])
    row = dict(alpha=rb['alpha'], b=rb['b'], gamma_closed=float(g),
               scipy_minus_closed_book=float(mp.mpf(rb['gamma']) - g),
               scipy_minus_closed_rerun=float(mp.mpf(rr['gamma']) - g),
               book_bracket_minus_exact=rb['bracket'] - rb['exact'],
               bracket_closed_minus_exact=float(br_cf - mp.mpf(rb['exact'])))
    out.append(row)
    print(f"{rb['alpha']:5.1f} {rb['b']:10.6f} {mp.nstr(g, 18):>22s} {row['scipy_minus_closed_book']:22.3e} "
          f"{row['scipy_minus_closed_rerun']:22.3e} {row['book_bracket_minus_exact']:21.3e} {row['bracket_closed_minus_exact']:19.3e}")
json.dump(out, open('archimedean_check.json', 'w'), indent=1)
