# New check script (item arith_exact): recomputes the ordinates of zeros_1700.txt (Arb, zeros_hp.py) with
# mpmath.zetazero at 75 digits, the routine behind the archived float64 lists (misc/repro/zeros.py) and the one
# Appendix app:code names for the ordinates of the zeros. One chunk of 50 zeros per array task.
# Usage: python3 zeros_mp_check.py TASK   (TASK = 0..33; writes zeros_mp_TASK.json)
import sys, json, mpmath as mp
task = int(sys.argv[1]); lo = 1 + 50 * task; hi = min(1700, lo + 49)
S = open('zeros_1700.txt').read().split()
mp.mp.dps = 75
out = []
for n in range(lo, hi + 1):
    z = mp.zetazero(n)
    a = mp.mpf(S[n - 1])
    out.append(dict(n=n, mpmath=mp.nstr(z.imag, 72), rel_diff=mp.nstr(abs(z.imag - a) / a, 3), re_minus_half=mp.nstr(z.real - mp.mpf(1) / 2, 3)))
worst = max(out, key=lambda r: mp.mpf(r['rel_diff']))
json.dump(dict(task=task, n_lo=lo, n_hi=hi, worst=worst, zeros=out), open('zeros_mp_%02d.json' % task, 'w'), indent=1)
print(task, lo, hi, 'max relative difference', worst['rel_diff'], 'at n =', worst['n'])
