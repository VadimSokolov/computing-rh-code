# Collects the outputs of v1_grid_mp.py and zeros_mp.py and compares them with riesz_program/data/v1min.json.
import json, glob, math
g = [json.load(open(f)) for f in sorted(glob.glob('grid_*.json'))]
z = [json.load(open(f)) for f in sorted(glob.glob('zeros_*.json'))]
ref = json.load(open('v1min.json'))
out = {'grid_tasks': len(g), 'grid_points': sum(r['points'] for r in g), 'b1': g[0]['b1'] if g else None}
if g:
    m = min(g, key=lambda r: r['min_excess'])
    out['min_excess'] = m['min_excess']; out['min_excess_at'] = m['min_excess_at']
    out['points_with_v1_le_b1'] = sum(r['points_with_v1_le_b1'] for r in g)
    m2 = min((r for r in g if r['min_on_13p5_2000'] is not None), key=lambda r: r['min_on_13p5_2000'])
    out['min_on_13p5_2000'] = m2['min_on_13p5_2000']; out['min_on_13p5_2000_at'] = m2['min_on_13p5_2000_at']
    out['reference_profile_min'] = ref['profile_min_on_13.5_T0']
if z:
    zs = sorted(t for r in z for t in r['mpmath'])
    gaps = [(zs[i + 1] - zs[i], zs[i], zs[i + 1]) for i in range(len(zs) - 1)]
    upto = [x for x in gaps if x[1] <= 2000 + 7]
    mg = max(upto)
    far = max(zs[0] - 13.5, max(x[0] for x in upto) / 2)
    out.update(zeros=len(zs), largest_ordinate=zs[-1], max_gap_below_2000=mg[0], max_gap_between=[mg[1], mg[2]],
               max_distance_to_nearest_zero=far, x_star=ref['x_star'], range2_ok=far < ref['x_star'],
               max_abs_diff_mpmath_arb=max(r['max_abs_diff'] for r in z),
               reference=dict(max_gap=ref['max_gap'], max_gap_between=ref['max_gap_between'], max_distance_to_nearest_zero=ref['max_distance_to_nearest_zero']))
json.dump(out, open('v1_check.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
