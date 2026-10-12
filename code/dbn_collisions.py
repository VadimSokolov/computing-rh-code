# dbn_collisions.py: the collisions of the real zeros of H_t under the backward de Bruijn Newman flow among the first 200 gaps,
# for Section sec:ch21:collisions and the section on Biane's theorem that follows it (ch/ch21.tex). Added for the book in October
# 2026, from the first author's note on the semicircle law and the Thorin framework (misc/reviews/nick-wigner/, version 7,
# Sections 7 and 8), recomputed from H_t itself (the note integrates the zero dynamics of the first 800 or 1200 zeros, a
# truncation that makes its collision times about 10^{-3} late).
# H_t(z) = (1/4) int_0^oo e^{t v^2/4} Phi(v) cos(z v/2) dv (Section sec:ch21:flow with v = 2u, Phi the kernel of eq:ch2:Phi), by
# the trapezoid rule with step 1/200 on [0, 3.6] at 200 digits. The integrand is analytic in |Im v| < pi/4, so for z <= 900
# and -2.05 <= t <= 0 the quadrature error is far below the size of H_t, which is about e^{-pi z/8}.
# The real zeros that start at z_j(0) = 2 gamma_j, j <= 215, are followed from t = 0 back to t = -2.05: an Euler predictor with
# dz/dt = H_zz/H_z and Newton's method at each step, with steps that keep every zero well inside its gaps. When two neighbours
# come within 0.05 of each other, their collision is located as the double zero H_t = dH_t/dz = 0 by Newton's method in (z, t)
# (H_t = -H_zz holds exactly for the trapezoid sum), recorded and removed. The top five zeros only follow and are dropped if
# they fail; collisions are reported for the pairs (k, k+1) with k <= 200 and counted for k <= 120 and k <= 200. The Newton
# steps for the zeros run in parallel.
# For each colliding pair the script also records, from H_0(z) = Xi(z/2)/8: the gap g = gamma_{k+1} - gamma_k, the isolated
# pair time -g^2/2, Laguerre's quantity Lambda = -(log H_0)'' (equal to Lambda_Xi(z/2)/4) at the critical point c of H_0 in the
# gap and its minimum m(0) over the gap, attained at u, the share P(0) = (u - z_k)^{-2} + (u - z_{k+1})^{-2} of the two zeros
# and R(0) = m(0) - P(0), the times tau_c = -1/Lambda(c) and tau_B = -1/m(0), and the frozen R predictor
# -(6 R(0))^{-1} log(m(0)/(m(0) - 6 R(0))) when m(0) > 6 R(0).
# Reference ordinates from Arb (acb.zeta_zeros). Run from any folder: python3 dbn_collisions.py; writes dbn_collisions.json
# there. Uses SLURM_CPUS_PER_TASK processes (default 1).
import json
import multiprocessing as mpr
import os
import time

import mpmath as mp
from flint import acb, ctx

t0 = time.time()
out = {}


def rec(k, v):
    out[k] = v
    print(k, v, round(time.time() - t0, 1), flush=True)
    json.dump(out, open('dbn_collisions.json', 'w'), indent=1)


def s(x, n=12):
    return mp.nstr(x, n)


KT = 215          # zeros tracked
KF = 5            # the top KF zeros only follow
KREP = 200        # collisions reported for the pairs (k, k+1), k <= KREP
DPS = 200
mp.mp.dps = DPS
ctx.prec = 600
GAM = [mp.mpf(z.imag.mid().str(160, radius=False)) for z in acb.zeta_zeros(1, KT)]
ctx.prec = 53
HD = mp.mpf(1) / 200


def phi_nodes(h, V, dps):
    # (v, Phi(v)) at v = j h, j = 0, ..., V/h, with Phi(v) = 2 sum_n (2 pi^2 n^4 e^{9v/2} - 3 pi n^2 e^{5v/2}) exp(-pi n^2 e^{2v})
    nodes = []
    with mp.workdps(dps):
        for j in range(int(mp.nint(V / h)) + 1):
            v = j * h
            q = mp.exp(-mp.pi * mp.exp(2 * v))
            a9 = mp.exp(9 * v / 2)
            a5 = mp.exp(5 * v / 2)
            tot = mp.mpf(0)
            qn = q
            step = q ** 3
            n = 1
            while True:   # qn = q^(n^2)
                term = (2 * mp.pi ** 2 * n ** 4 * a9 - 3 * mp.pi * n ** 2 * a5) * qn
                tot += term
                if n > 3 and abs(term) < abs(tot) * mp.eps * mp.mpf(10) ** -20:
                    break
                qn *= step
                step *= q * q
                n += 1
            nodes.append((v, 2 * tot))
    return nodes


NODES = phi_nodes(HD, mp.mpf('3.6'), DPS + 20)
U = [v / 2 for v, _ in NODES]
U2 = [u * u for u in U]
U3 = [u * u * u for u in U]
_cache = {}


def weights(t):
    t = mp.mpf(t)
    if _cache.get('t') != t:
        _cache['t'] = t
        _cache['W'] = [HD * ph * mp.exp(t * v * v / 4) / 4 * (mp.mpf(1) / 2 if v == 0 else 1) for v, ph in NODES]
    return _cache['W']


def parts(z, W):
    # H, H_z, H_zz, H_zzz at real z for the weights W of one t; cos and sin of j z h/2 by rotation
    x = z * HD / 2
    c1, s1 = mp.cos(x), mp.sin(x)
    c, sn = mp.mpf(1), mp.mpf(0)
    H = Hz = Hzz = Hzzz = mp.mpf(0)
    for j in range(len(W)):
        wc = W[j] * c
        ws = W[j] * sn
        H += wc
        Hz -= U[j] * ws
        Hzz -= U2[j] * wc
        Hzzz += U3[j] * ws
        c, sn = c * c1 - sn * s1, sn * c1 + c * s1
    return H, Hz, Hzz, Hzzz


def newton_job(args):
    z, t = args
    W = weights(t)
    P = None
    tol = mp.mpf(10) ** -30
    for _ in range(40):
        P = parts(z, W)
        dz = P[0] / P[1]
        z -= dz
        if abs(dz) < tol * (1 + abs(z)):
            return z, True, P
    return z, False, P


def double_zero(z, t, maxit=40):
    for _ in range(maxit):
        H, Hz, Hzz, Hzzz = parts(z, weights(t))
        a, b, c, d = Hz, -Hzz, Hzz, -Hzzz
        det = a * d - b * c
        dz = (d * H - b * Hz) / det
        dt = (-c * H + a * Hz) / det
        z -= dz
        t -= dt
        if abs(dz) < mp.mpf(10) ** -25 and abs(dt) < mp.mpf(10) ** -25:
            return z, t, parts(z, weights(t)), True
    return z, t, None, False


def lam(z):
    H, Hz, Hzz, _ = parts(z, weights(0))
    return (Hz * Hz - H * Hzz) / (H * H)


def gap_job(k):
    # Laguerre's quantity on the gap (z_k, z_{k+1}) of H_0, k from 1; z_j = 2 gamma_j
    a, b = 2 * GAM[k - 1], 2 * GAM[k]
    g = GAM[k] - GAM[k - 1]
    # critical point: H_z changes sign in the gap; bisection then Newton
    lo, hi = a + (b - a) / 1000, b - (b - a) / 1000
    flo = parts(lo, weights(0))[1]
    for _ in range(60):
        mid = (lo + hi) / 2
        fm = parts(mid, weights(0))[1]
        if (fm > 0) == (flo > 0):
            lo, flo = mid, fm
        else:
            hi = mid
    c = (lo + hi) / 2
    for _ in range(8):
        P = parts(c, weights(0))
        c -= P[1] / P[2]
    # minimum of Lambda over the gap (Lambda is convex there): golden section
    gr = (mp.sqrt(5) - 1) / 2
    lo, hi = a + (b - a) / 1000, b - (b - a) / 1000
    x1, x2 = hi - gr * (hi - lo), lo + gr * (hi - lo)
    f1, f2 = lam(x1), lam(x2)
    for _ in range(110):
        if f1 < f2:
            hi, x2, f2 = x2, x1, f1
            x1 = hi - gr * (hi - lo)
            f1 = lam(x1)
        else:
            lo, x1, f1 = x1, x2, f2
            x2 = lo + gr * (hi - lo)
            f2 = lam(x2)
    u = (lo + hi) / 2
    m0 = lam(u)
    P0 = 1 / (u - a) ** 2 + 1 / (u - b) ** 2
    R0 = m0 - P0
    pred = -mp.log(m0 / (m0 - 6 * R0)) / (6 * R0) if m0 > 6 * R0 else None
    return {'k': k, 'gap': g, 'isolated': -g * g / 2, 'c': c, 'Lambda(c)': lam(c), 'u': u, 'm0': m0, 'P0': P0, 'R0': R0,
            'tau_c': -1 / lam(c), 'tau_B': -1 / m0, 'pred': pred}


if __name__ == '__main__':
    nproc = int(os.environ.get('SLURM_CPUS_PER_TASK', '1'))
    pool = mpr.get_context('fork').Pool(nproc)
    t = mp.mpf(0)
    res = pool.map(newton_job, [(2 * g, t) for g in GAM])
    Z = [r[0] for r in res]
    P = [r[2] for r in res]
    ids = list(range(1, KT + 1))
    rec('processes', nproc)
    rec('largest |z_j(0) - 2 gamma_j|, j <= %d' % KT, s(max(abs(Z[k] - 2 * GAM[k]) for k in range(KT)), 3))
    TEND = mp.mpf('-2.05')
    GTH = mp.mpf('0.05')
    cols, lost = [], []
    steps = 0
    top_ctrl = KT
    while t > TEND and len(Z) > KF + 2:
        steps += 1
        V = [p_[2] / p_[1] for p_ in P]
        ctrl = len(Z) - KF
        top_ctrl = min(top_ctrl, ids[ctrl])
        close = [(Z[k + 1] - Z[k], k) for k in range(ctrl) if Z[k + 1] - Z[k] < GTH and V[k + 1] - V[k] > 0]
        tn = None
        if close:
            g, k = min(close)
            zc, tc, Pc, ok = double_zero((Z[k] + Z[k + 1]) / 2, t - g * g / 8)
            if ok and Z[k] - g < zc < Z[k + 1] + g and t - mp.mpf('0.01') < tc < t:
                cols.append({'ids': [ids[k], ids[k + 1]], 'zc': zc, 'tc': tc, 'res': [Pc[0], Pc[1]]})
                print('collision', ids[k], ids[k + 1], s(zc, 10), s(tc, 10), round(time.time() - t0, 1), flush=True)
                del Z[k:k + 2]
                del ids[k:k + 2]
                del P[k:k + 2]
                V = [p_[2] / p_[1] for p_ in P]
                tn = tc - mp.mpf('1e-9')
                ctrl = len(Z) - KF
            else:
                tn = t - g * g / 32
        if tn is None:
            dt = mp.mpf('0.05')
            for k in range(ctrl):
                r = V[k + 1] - V[k]
                if r > 0:
                    dt = min(dt, mp.mpf('0.3') * (Z[k + 1] - Z[k]) / r)
            for k in range(ctrl):
                gl = Z[k] - Z[k - 1] if k > 0 else Z[0]
                if V[k] != 0:
                    dt = min(dt, mp.mpf('0.2') * min(gl, Z[k + 1] - Z[k]) / abs(V[k]))
            tn = t - dt
        for attempt in range(40):
            res = pool.map(newton_job, [(z + V[k] * (tn - t), tn) for k, z in enumerate(Z)])
            bad = None
            for k, r in enumerate(res):
                if not r[1]:
                    bad = k
                    break
            if bad is None:
                for k in range(len(res) - 1):
                    if not (res[k + 1][0] - res[k][0] > mp.mpf('1e-12')) or (res[k][2][1] > 0) == (res[k + 1][2][1] > 0):
                        bad = k + 1
                        break
            if bad is None:
                break
            if bad >= ctrl:
                lost.append({'ids': ids[bad:], 't': s(tn, 8)})
                Z, ids, P, V = Z[:bad], ids[:bad], P[:bad], V[:bad]
                ctrl = len(Z) - KF
                continue
            tn = t + (tn - t) / 2
        else:
            rec('tracking failed at t', s(t, 10))
            break
        t = tn
        Z = [r[0] for r in res]
        P = [r[2] for r in res]
    rep = [c_ for c_ in cols if c_['ids'][0] <= KREP and c_['ids'][1] == c_['ids'][0] + 1]
    other = [c_ for c_ in cols if c_ not in rep]
    gaps = dict((d['k'], d) for d in pool.map(gap_job, [c_['ids'][0] for c_ in rep]))
    rows = []
    for c_ in sorted(rep, key=lambda c_: -c_['tc']):
        k = c_['ids'][0]
        d = gaps[k]
        rows.append({'pair': c_['ids'], 'z': s(c_['zc'], 12), 'theta = z/2': s(c_['zc'] / 2, 10), 't_coll': s(c_['tc'], 12),
                     'H, H_z there': [s(c_['res'][0], 3), s(c_['res'][1], 3)], 'gap gamma_{k+1} - gamma_k': s(d['gap'], 8),
                     'isolated -gap^2/2': s(d['isolated'], 8), 'critical point c': s(d['c'], 12), 'Lambda(c)': s(d['Lambda(c)'], 10),
                     'tau_c = -1/Lambda(c)': s(d['tau_c'], 8), 'minimiser u': s(d['u'], 12), 'm(0)': s(d['m0'], 10), 'P(0)': s(d['P0'], 10),
                     'R(0)': s(d['R0'], 10), 'tau_B = -1/m(0)': s(d['tau_B'], 8), '6 R(0)/m(0)': s(6 * d['R0'] / d['m0'], 6),
                     'frozen R predictor': s(d['pred'], 8) if d['pred'] is not None else None,
                     't_coll/tau_B': s(c_['tc'] / d['tau_B'], 6), 't_coll/tau_c': s(c_['tc'] / d['tau_c'], 6), 't_coll/isolated': s(c_['tc'] / d['isolated'], 6),
                     'predictor error (relative)': s(d['pred'] / c_['tc'] - 1, 4) if d['pred'] is not None else None})
    rec('collisions among the first %d gaps on [-2.05, 0], latest first' % KREP, rows)
    rec('summary', {'steps': steps, 'final t': s(t, 8), 'tracked zeros left': len(Z), 'dropped at the top': lost,
                    'collisions among the first 120 gaps': sum(1 for c_ in rep if c_['ids'][0] <= 120),
                    'of them at t > -1.6 (first 120 gaps)': sum(1 for c_ in rep if c_['ids'][0] <= 120 and c_['tc'] > mp.mpf('-1.6')),
                    'collisions among the first %d gaps' % KREP: len(rep),
                    'of them at t > -1.6 (first %d gaps)' % KREP: sum(1 for c_ in rep if c_['tc'] > mp.mpf('-1.6')),
                    'pairs (k, k+1) controlled throughout for k <': top_ctrl,
                    'other collisions (pairs of nonadjacent indices, or k > %d)' % KREP: [[c_['ids'], s(c_['zc'], 9), s(c_['tc'], 9)] for c_ in other]})
    out['seconds'] = time.time() - t0
    json.dump(out, open('dbn_collisions.json', 'w'), indent=1)
    print('done', round(time.time() - t0, 1))
