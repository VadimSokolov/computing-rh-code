# The de Bruijn Newman flow for Ramanujan's Delta (Section 8.5.5). H_lam(t) = int_R e^{lam u^2} K(u) cos(t u) du with
# K(u) = e^{6u} Delta(i e^u), even in u, so H_0(t) = xi(1/2 + it, Delta) in the normalisation of Section 8.5.2.
# Quadrature: Gauss Legendre on [0, 4] (K is below 1e-130 beyond 4), 8 panels of 96 nodes, kernel values computed once,
# checked against 192 nodes per panel and against Hecke's formula. Forward (lam = 0.05, ..., 0.5) the zeros below 52 are
# followed by Newton's method. Backward, the critical point of H_lam between each adjacent pair of zeros is followed in lam;
# where the critical value changes sign two zeros have met, located by bisection in lam. Each such event is then tested
# directly: a genuine collision loses two sign changes of H_lam near t_c between lam_c + 0.005 and lam_c - 0.005, and the
# colliding pair is identified by counting the sign changes below t_c - 1.2 (real zeros keep their order). Finally the
# sign changes on (0.5, 52] are counted at eleven values of lam, and the highest zero is followed back until it passes 52.
# Writes dbn_delta.json.
import json, time
import mpmath as mp

mp.mp.dps = 40
t0 = time.time()
NT = 40
sig = [0]*(NT+1)
for d in range(1, NT+1):
    for m in range(d, NT+1, d): sig[m] += d
tau = [0]*(NT+1); tau[1] = 1
for n in range(2, NT+1):
    tau[n] = -24*sum(sig[m]*tau[n-m] for m in range(1, n))//(n-1)
assert tau[2] == -24 and tau[3] == 252 and tau[11] == 534612
def Kdelta(v):
    y = mp.exp(v)
    return mp.exp(6*v)*mp.fsum(tau[n]*mp.exp(-2*mp.pi*n*y) for n in range(1, NT+1))
gl = mp.calculus.quadrature.GaussLegendre(mp.mp)
def make_nodes(degree, panels=8, L=4):
    base = gl.calc_nodes(degree, mp.mp.prec)   # 3*2^(degree-1) nodes on [-1, 1]
    nodes = []
    h = mp.mpf(L)/panels
    for k in range(panels):
        a = k*h; b = a + h
        for x, w in base:
            v = (a+b)/2 + h/2*x
            nodes.append((v, v*v, h*w*Kdelta(v)))   # h/2 for the panel, times 2 for the even extension to R
    return nodes
NODES = make_nodes(6)
def weights(lam, nodes=None):
    nodes = nodes or NODES
    return [(v, c*mp.exp(lam*v2)) for v, v2, c in nodes]
def H(W, t, der=0):
    s0 = s1 = s2 = mp.mpf(0)
    for v, w in W:
        c, s = mp.cos_sin(t*v)
        s0 += w*c
        if der:
            s1 -= w*v*s; s2 -= w*v*v*c
    return (s0, s1, s2) if der else s0
def nsign(W, lo, hi, step):
    lo, hi, step = mp.mpf(lo), mp.mpf(hi), mp.mpf(step)
    n = int(mp.ceil((hi-lo)/step)); cnt = 0; f0 = H(W, lo)
    for k in range(1, n+1):
        f1 = H(W, lo + k*step)
        if f0*f1 < 0: cnt += 1
        f0 = f1
    return cnt

out = {}
# quadrature checks at lam = 0
W0 = weights(mp.mpf(0))
W0fine = weights(mp.mpf(0), make_nodes(7))
def hecke(t):
    s = mp.mpc(6, t)
    return mp.re(mp.fsum(tau[n]*((2*mp.pi*n)**(-s)*mp.gammainc(s, 2*mp.pi*n) + (2*mp.pi*n)**(s-12)*mp.gammainc(12-s, 2*mp.pi*n)) for n in range(1, NT+1)))
chk = {}
for t in [0, 10, 30, 50]:
    a = H(W0, mp.mpf(t)); b = H(W0fine, mp.mpf(t)); h = hecke(mp.mpf(t))
    chk[str(t)] = {'H_0': mp.nstr(a, 20), 'rel_diff_192_nodes': mp.nstr(abs(a-b)/abs(h), 3), 'rel_diff_hecke': mp.nstr(abs(a-h)/abs(h), 3)}
out['quadrature_check_lam0'] = chk
print(chk, round(time.time()-t0, 1), flush=True)
del W0fine

TMAX = 52
def newton_zero(W, t):
    for _ in range(50):
        f, f1, _f2 = H(W, t, der=1)
        dt = f/f1
        if abs(dt) > mp.mpf('0.3'): dt = mp.sign(dt)*mp.mpf('0.3')
        t -= dt
        if abs(dt) < mp.mpf(10)**(-25): break
    return t
def zeros_at(W, step=mp.mpf('0.05'), tmax=TMAX):
    zs = []; t = mp.mpf('0.5'); f0 = H(W, t)
    while t < tmax:
        t1 = t + step; f1 = H(W, t1)
        if f0*f1 < 0:
            zs.append(mp.findroot(lambda x: H(W, x), (t, t1), solver='anderson'))
        t, f0 = t1, f1
    return zs
Z0 = zeros_at(W0)
out['zeros_lam0'] = [mp.nstr(z, 12) for z in Z0]
print('zeros at lam = 0:', out['zeros_lam0'], round(time.time()-t0, 1), flush=True)

# forward
fwd = {}
Zc = list(Z0)
for i in range(1, 11):
    lam = mp.mpf(i)/20; W = weights(lam)
    Zc = [newton_zero(W, z) for z in Zc]
    fwd[mp.nstr(lam, 3)] = {'zeros': [mp.nstr(z, 10) for z in Zc], 'sign_changes_up_to_33': len(zeros_at(W, tmax=33)),
                            'followed_zeros_up_to_33': len([z for z in Zc if z <= 33])}
    print('lam', mp.nstr(lam, 3), fwd[mp.nstr(lam, 3)]['sign_changes_up_to_33'], round(time.time()-t0, 1), flush=True)
out['forward'] = fwd
mono = []
lams = sorted(fwd, key=float)
for k in range(len(Z0)):
    seq = [float(Z0[k])] + [float(fwd[l]['zeros'][k]) for l in lams]
    dif = [b-a for a, b in zip(seq[:-1], seq[1:])]
    mono.append('up' if all(d > 0 for d in dif) else ('down' if all(d < 0 for d in dif) else 'not monotone'))
out['forward_direction_per_zero'] = mono

# backward: candidate collisions from the critical points
def crit(W, t):
    for _ in range(60):
        f, f1, f2 = H(W, t, der=1)
        dt = f1/f2
        if abs(dt) > mp.mpf('0.3'): dt = mp.sign(dt)*mp.mpf('0.3')
        t -= dt
        if abs(dt) < mp.mpf(10)**(-25): break
    return t, H(W, t)
state = []
for a, b in zip(Z0[:-1], Z0[1:]):
    tc, val = crit(W0, (a+b)/2)
    state.append({'tc': tc, 'sign': mp.sign(val), 'event': None})
lam = mp.mpf(0); step = mp.mpf('0.01')
while lam > -2 and any(s['event'] is None for s in state):
    lam1 = lam - step; W = weights(lam1)
    for s in state:
        if s['event'] is not None: continue
        tc, val = crit(W, s['tc'])
        if mp.sign(val) != s['sign']:
            lo, hi, thi = lam1, lam, s['tc']
            for _ in range(25):
                mid = (lo+hi)/2
                tm, vm = crit(weights(mid), thi)
                if mp.sign(vm) == s['sign']: hi, thi = mid, tm
                else: lo = mid
            s['event'] = ((lo+hi)/2, thi)
        else:
            s['tc'] = tc
    lam = lam1
events = []
for s in state:
    if s['event'] is None: continue
    lc, tc = s['event']
    if not any(abs(lc - e[0]) < mp.mpf('1e-4') and abs(tc - e[1]) < mp.mpf('0.05') for e in events):
        events.append((lc, tc))
events.sort(key=lambda e: -e[0])
print('candidate events', len(events), round(time.time()-t0, 1), flush=True)

# test each event, in the order in which it happens as lam decreases
alive = list(range(len(Z0)))     # indices into Z0 of the zeros not yet collided
coll, rejected = [], []
eps = mp.mpf('0.005'); win = mp.mpf('1.2')
for lc, tc in events:
    Wa, Wb = weights(lc + eps), weights(lc - eps)
    above = nsign(Wa, tc - win, tc + win, '0.004'); below = nsign(Wb, tc - win, tc + win, '0.004')
    rec = {'lam_c': mp.nstr(lc, 8), 't_c': mp.nstr(tc, 8), 'local_sign_changes_at_lam_c_plus_0.005': above, 'at_lam_c_minus_0.005': below}
    if above == below + 2 and tc + win < TMAX:
        k = nsign(Wa, '0.5', tc - win, '0.01')
        i, j = alive[k], alive[k+1]
        rec['pair_at_lam0'] = [mp.nstr(Z0[i], 10), mp.nstr(Z0[j], 10)]
        alive = alive[:k] + alive[k+2:]
        coll.append(rec)
    else:
        rejected.append(rec)
    print(rec, round(time.time()-t0, 1), flush=True)
out['collisions'] = coll
out['rejected_events'] = rejected
# global counts of sign changes on (0.5, 33] and (0.5, 52]
cnts = {}
for lam in ['0', '-0.1', '-0.2', '-0.3', '-0.4', '-0.45', '-0.5', '-0.75', '-1', '-1.25', '-1.5']:
    W = weights(mp.mpf(lam))
    c33 = nsign(W, '0.5', '33', '0.01'); c52 = c33 + nsign(W, '33', TMAX, '0.01')
    cnts[lam] = {'up_to_33': c33, 'up_to_52': c52}
    print('lam', lam, cnts[lam], round(time.time()-t0, 1), flush=True)
out['sign_change_counts'] = cnts
# the highest zero below TMAX followed backward by Newton's method, with a sign change check, until it passes TMAX
top = []
z = Z0[-1]
for i in range(1, 101):
    lam = -mp.mpf(i)/100; W = weights(lam)
    z = newton_zero(W, z)
    ok = H(W, z - mp.mpf('0.002'))*H(W, z + mp.mpf('0.002')) < 0
    top.append({'lam': mp.nstr(lam, 3), 'zero': mp.nstr(z, 10), 'sign_change': bool(ok)})
    if z > TMAX or not ok: break
out['highest_zero_backward'] = top
print('highest zero', top[-3:], round(time.time()-t0, 1), flush=True)
out['seconds'] = time.time()-t0
json.dump(out, open('dbn_delta.json', 'w'), indent=1)
print('done', round(time.time()-t0, 1))
