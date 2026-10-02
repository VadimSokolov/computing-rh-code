# Independent check of Proposition bp:prop:v1min with mpmath (not Arb): evaluates
# v_1(theta) = Re (xi'/xi)(1 + i theta) on theta = STEP*k, 1 <= k <= KMAX, in chunk TASK of NTASK,
# and records the smallest excess v_1 - b_1 overall and the smallest value on [13.5, 2000].
import json, sys
import mpmath as mp
mp.mp.dps = 25
TASK, NTASK = int(sys.argv[1]), int(sys.argv[2])
STEP, KMAX = mp.mpf('0.005'), 420000                  # grid up to theta = 2100
b1 = 1 + mp.euler / 2 - mp.log(4 * mp.pi) / 2
def v1(th):
    s = mp.mpc(1, th)
    return (1 / (1 + th * th) - mp.log(mp.pi) / 2 + mp.re(mp.digamma(s / 2)) / 2
            + mp.re(mp.zeta(s, derivative=1) / mp.zeta(s)))
lo = 1 + TASK * KMAX // NTASK; hi = (TASK + 1) * KMAX // NTASK
best, best2, nbad = (mp.inf, None), (mp.inf, None), 0
for k in range(lo, hi + 1):
    th = STEP * k
    e = v1(th) - b1
    if e <= 0:
        nbad += 1
    if e < best[0]:
        best = (e, th)
    if 13.5 <= th <= 2000 and e + b1 < best2[0]:
        best2 = (e + b1, th)
out = dict(task=TASK, theta_lo=float(STEP * lo), theta_hi=float(STEP * hi), points=hi - lo + 1, b1=float(b1),
           min_excess=float(best[0]), min_excess_at=float(best[1]), points_with_v1_le_b1=nbad,
           min_on_13p5_2000=None if best2[1] is None else float(best2[0]), min_on_13p5_2000_at=None if best2[1] is None else float(best2[1]))
json.dump(out, open('grid_%03d.json' % TASK, 'w'))
print(json.dumps(out))
