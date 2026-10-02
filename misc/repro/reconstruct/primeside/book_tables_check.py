"""RECONSTRUCTION AID, NOT THE AUTHORS' CODE.

Checks the numbers printed in Table tab:prime (ch/ch02.tex, lines 160 to 165) and Table
tab:split (ch/ch06.tex, lines 40 to 45, and the text at lines 28 and 32) against a
results_prime.json and an extra.json, by rounding each computed value to the printed digits.
The printed values are copied here by hand from the chapters as they stood on 2026-09-26.
Usage: python book_tables_check.py results_prime.json extra.json
"""
import json, sys

prime = {(r['alpha'], r['b']): r for r in json.load(open(sys.argv[1]))}
extra = {(r['alpha'], r['b']): r for r in json.load(open(sys.argv[2]))['prime']}
G1 = 14.134725
# Table tab:prime: alpha, theta, sine squared bracket, Polya route, error of prime side
tab_prime = [(3.0, G1, '0.561447238789', '0.561447238789', '1.1e-13'),
             (3.0, 50.0, '1.120967014740', '1.120967015289', '5.5e-10'),
             (2.0, G1, '0.766398694727', '0.766399169607', '4.7e-07'),
             (2.0, 50.0, '1.224475714058', '1.224476216908', '5.0e-07'),
             (1.5, G1, '1.065762504228', '1.067140080666', '1.4e-03'),
             (1.5, 50.0, '1.385288214458', '1.386706213218', '1.4e-03')]
# Table tab:split: alpha, theta, drift, archimedean, pole (signed), primes, v, kappa
tab_split = [(3.0, 2.0, '0.1144', '0.0533', '-0.2500', '0.1988', '0.1165', '3.85'),
             (3.0, 50.0, '0.1144', '1.2594', '-0.4992', '0.2464', '1.1210', '0.67'),
             (2.0, 2.0, '0.0691', '0.0859', '-0.8000', '0.7154', '0.0704', '21.5'),
             (2.0, 50.0, '0.0691', '1.3989', '-0.9996', '0.7561', '1.2245', '1.43'),
             (1.5, 2.0, '0.0461', '0.1137', '-1.8824', '1.7683', '0.0470', '77.6'),
             (1.5, 50.0, '0.0461', '1.4863', '-1.9998', '1.8527', '1.3867', '2.78')]


def fmt(x, printed):
    if 'e' in printed:
        return f'{x:.1e}'
    d = len(printed.split('.')[1]) if '.' in printed else 0
    return f'{x:.{d}f}'


n = bad = 0
for src_name, src in [('results_prime.json', prime), ('extra.json', extra)]:
    print(f'== against {src_name}')
    for a, t, br, po, er in tab_prime:
        r = src[(a, t)]
        for lab, val, pr in [('bracket', r['bracket'], br), ('polya', r['polya'], po), ('error', abs(r['bracket'] - r['polya']), er)]:
            n += 1; ok = fmt(val, pr) == pr; bad += not ok
            if not ok: print(f'   tab:prime ({a}, {t}) {lab}: printed {pr}, computed {val!r} -> {fmt(val, pr)}')
    for a, t, c, g, p, q, v, k in tab_split:
        r = src[(a, t)]
        kap = (abs(r['pole']) + abs(r['primes'])) / abs(r['polya'])
        for lab, val, pr in [('drift', r['c'], c), ('archimedean', r['gamma'], g), ('pole', r['pole'], p),
                             ('primes', r['primes'], q), ('v', r['polya'], v), ('kappa', kap, k)]:
            n += 1; ok = fmt(val, pr) == pr; bad += not ok
            if not ok: print(f'   tab:split ({a}, {t}) {lab}: printed {pr}, computed {val!r} -> {fmt(val, pr)}')
    r = src[(1.5, 2.0)]
    kap = {a: (abs(src[(a, 2.0)]['pole']) + abs(src[(a, 2.0)]['primes'])) / abs(src[(a, 2.0)]['polya']) for a in (3.0, 2.0, 1.5)}
    for lab, val, pr in [('ch06:32 kappa_3(2)', kap[3.0], '3.9'), ('ch06:32 kappa_2(2)', kap[2.0], '21.5'),
                         ('ch06:32 kappa_1.5(2)', kap[1.5], '77.6'), ('ch06:28 tail at (1.5, 2)', r['polya'] - r['bracket'], '1.3e-03')]:
        n += 1; ok = fmt(val, pr) == pr; bad += not ok
        if not ok: print(f'   {lab}: printed {pr}, computed {val!r} -> {fmt(val, pr)}')
print(f'{n} printed numbers checked (each table against both files), {bad} disagree after rounding')
