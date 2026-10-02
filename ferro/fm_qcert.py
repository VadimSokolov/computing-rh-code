# Chapter 13, Table 13.6 and Problem 13.4: certified positivity of Q_{K,v} on the whole line, with certified inputs.
#   Q_{K,v}(y) = sum_{j<=K} e_j v^-j He_2j(y),  e_j the elementary symmetric functions of gamma_1^-2, ..., gamma_K^-2,
# and f_{K,v}(x) = g_v(x) Q_{K,v}(x/sqrt v), f_K = f_{K, 2 beta_K}, 2 beta_K = Var X - 2 s_K (Proposition fm:prop:hermite).
# Inputs. The ordinates gamma_1..gamma_N are computed by Arb (acb.zeta_zeros), which isolates the zeros, counts them by
# Turing's method and refines them in ball arithmetic, and Var X = (log xi)''(1/2) comes from the Taylor series of
# log xi at 1/2 (no hypothesis); both at PIN bits, written to fm_zeros_hp.txt by the mode 'zeros'.
# Method. Q is built in the monomial basis in ball arithmetic at a working precision W < PIN (the Hermite sum cancels
# catastrophically, so W is thousands of bits for large K). On a dyadic interval [c - r, c + r] the Taylor coefficients
# T_m = Q^(m)(c)/m! are computed exactly in ball arithmetic (Taylor shift by one polynomial product), and
#   Q(c + h) >= T_0 - sum_{m>=1} |T_m| r^m  for |h| <= r,
# a polynomial identity with no remainder. Intervals where this lower bound is not positive are bisected. The intervals
# cover [0, Y] with Y >= sqrt(8K+2); for |y| >= sqrt(8K+2), Q > 0 is proved in the book, and Q is even, so a complete
# cover proves Q_{K,v} > 0 on the real line, that is f_{K,v} > 0.
# For v = t * 2 beta_K with t a dyadic rational: t = 1 is f_K itself; for the threshold v*(K) of the proposition the
# mode 'table' covers [0, Y] at t_hi (so v*(K) <= t_hi 2 beta_K, positivity being monotone in v) and finds a point with
# Q_{K, t_lo 2 beta_K}(y) < 0 certified (so v*(K) > t_lo 2 beta_K), starting from the bracket of the grid search in
# fm_hermite.json and bisecting with complete covers to width 2^-18 where that bracket is wrong. At t = 1 it also refines
# the cover until the certified lower bound of min Q is within 2e-4 of Q(0).
# Usage: python3 fm_qcert.py zeros N PIN       writes fm_zeros_hp.txt
#        python3 fm_qcert.py table             reads fm_zeros_hp.txt, fm_hermite.json; writes fm_qcert_table.json
#        python3 fm_qcert.py all N             reads fm_zeros_hp.txt; f_K > 0 for every K = 1..N; writes fm_qcert_all.json
#        python3 fm_qcert.py minq 1-200,300    reads fm_zeros_hp.txt; min of Q_{K, 2 beta_K} on the line is Q(0), attained
#                                              only at y = 0, for each listed K; writes fm_qcert_minq.json
import json, math, os, sys, time
from decimal import Decimal
from fractions import Fraction
from multiprocessing import Pool
from flint import arb, acb, arb_poly, acb_series, ctx


def dec(x, n, up):
    """An n significant digit decimal string at most every point of the ball x (up=False) or at least every point
    (up=True): the exact binary endpoint, rounded outward in exact rational arithmetic."""
    m, k = (int(v) for v in (x.upper() if up else x.lower()).man_exp())
    q = Fraction(m)*Fraction(2)**k
    if q == 0:
        return '0'
    d = len(str(abs(q.numerator))) - len(str(q.denominator))
    while Fraction(10)**d > abs(q):
        d -= 1
    while Fraction(10)**(d + 1) <= abs(q):
        d += 1
    j = q/Fraction(10)**(d - n + 1)
    return str(Decimal(math.ceil(j) if up else math.floor(j)).scaleb(d - n + 1))


CPUS = int(os.environ.get('SLURM_CPUS_PER_TASK', '8'))
S0 = 5                          # initial intervals of width 2^-S0
SMAX = 26                       # finest width 2^-SMAX
GOAL = Fraction(2, 10000)       # at t = 1, refine until min lower bound >= Q(0) - GOAL


def zeros_chunk(args):
    start, num, pin = args
    ctx.prec = pin
    Z = acb.zeta_zeros(start, num)
    half = arb(1)/2
    assert all(z.real == half for z in Z)
    return [(start + i, z.imag.str(int(pin*0.30103) + 10, radius=True)) for i, z in enumerate(Z)]


def var_x(pin):
    ctx.prec = pin
    x = acb_series([acb(arb(1)/2), acb(1)], prec=3)
    lx = x.log() + (1 - x).log() - x*(arb.pi().log()/2) + (x/2).lgamma() + (-x.zeta()).log()
    return 2*lx.coeffs()[2].real


def mode_zeros(N, pin):
    t0 = time.time()
    V = var_x(pin)
    step = 10
    tasks = [(s, min(step, N - s + 1), pin) for s in range(1, N + 1, step)]
    with Pool(CPUS) as pool:
        res = [r for chunk in pool.imap(zeros_chunk, tasks) for r in chunk]
    assert [k for k, _ in res] == list(range(1, N + 1))
    ctx.prec = pin
    G = [arb(s) for _, s in res]
    assert all(G[i + 1] > G[i] for i in range(N - 1))
    maxrad = max(float(g.rad()) for g in G)
    with open('fm_zeros_hp.txt', 'w') as f:
        f.write('# prec %d; Var X = (log xi)\'\'(1/2), then gamma_1 .. gamma_%d (Arb, acb.zeta_zeros)\n' % (pin, N))
        f.write(V.str(int(pin*0.30103) + 10, radius=True) + '\n')
        for _, s in res:
            f.write(s + '\n')
    print('zeros: N = %d at %d bits, max radius %.3e, Var X = %s, %.0f s'
          % (N, pin, maxrad, V.str(30), time.time() - t0), flush=True)


def read_inputs(W):
    """Var X and the ordinates as balls at precision W, from fm_zeros_hp.txt."""
    ctx.prec = W
    L = [l.strip() for l in open('fm_zeros_hp.txt') if l.strip() and not l.startswith('#')]
    return arb(L[0]), [arb(s) for s in L[1:]]


def default_prec(K, t):
    """Working precision, from the size of the Hermite terms (as in fm_hermite.py) with a margin for the
    monomial basis and the Taylor shift."""
    if K == 0:
        return 128
    g = [float(x) for x in GF[:K]]
    S = sum(x**-2 for x in g)
    v = float(t)*(VARF - 2*S)
    return int(1.12*3.33*(30 + K*max(0.0, math.log10(32*S/v)))) + 256


class Q_data:
    """Q_{K, t 2beta_K} in the monomial basis, at working precision W."""

    def __init__(self, K, t, W):
        self.K, self.t, self.W = K, t, W
        ctx.prec = W
        V, G = read_inputs(W)
        w = [1/(g*g) for g in G[:K]]
        e = [arb(1)] + [arb(0)]*K
        for wi in w:
            for j in range(K, 0, -1):
                e[j] += wi*e[j - 1]
        S = sum(w, arb(0))
        two_beta = V - 2*S
        v = two_beta*arb(t.numerator)/arb(t.denominator)
        c = [e[j]/v**j for j in range(K + 1)]
        y = arb_poly([0, 1])
        H0, H1 = arb_poly([1]), y
        Q = arb_poly([c[0]])
        for n in range(1, 2*K):
            H0, H1 = H1, y*H1 - n*H0
            if n & 1:
                Q += H1*c[(n + 1) >> 1]
        self.N = N = 2*K
        b = Q.coeffs() + [arb(0)]*(N + 1 - len(Q.coeffs()))
        fact = [1]
        for n in range(1, N + 1):
            fact.append(fact[-1]*n)
        self.U = arb_poly([b[N - i]*fact[N - i] for i in range(N + 1)])    # U[i] = b_{N-i} (N-i)!
        self.b, self.fact = b, fact
        self.two_beta, self.S, self.v = two_beta, S, v
        self.rho = {}

    def value(self, yq):
        """Q at the rational point yq, a ball (Horner)."""
        ctx.prec = self.W
        Y = arb(yq.numerator)/arb(yq.denominator)
        s = arb(0)
        for n in range(self.N, -1, -1):
            s = s*Y + self.b[n]
        return s

    def lower(self, i, s):
        """Lower bound of Q on [i 2^-s, (i+1) 2^-s], and the value at the centre."""
        ctx.prec = self.W
        N = self.N
        c = arb(2*i + 1)/arb(2)**(s + 1)
        if N == 0:
            return self.b[0], self.b[0]
        wk = [arb(1)]
        for k in range(1, N + 1):
            wk.append(wk[-1]*c/k)
        P = (self.U*arb_poly(wk)).coeffs()
        P = P + [arb(0)]*(2*N + 1 - len(P))
        if s not in self.rho:                      # r^m / m!, r = 2^-(s+1)
            r = arb(1)/arb(2)**(s + 1)
            rho = [arb(1)]
            for m in range(1, N + 1):
                rho.append(rho[-1]*r/m)
            self.rho[s] = rho
        rho = self.rho[s]
        T0 = P[N]
        tail = arb(0)
        for m in range(1, N + 1):
            tail += abs(P[N - m])*rho[m]
        return T0 - tail, T0


_cache = {}


def get_data(K, t, W):
    key = (K, t, W)
    if key not in _cache:
        _cache.clear()
        _cache[key] = Q_data(K, t, W)
    return _cache[key]


def cover(task):
    """Cover the initial intervals i0 <= i < i1 (width 2^-S0). Returns the leaves' statistics."""
    K, t, W, i0, i1, goal = task
    D = get_data(K, t, W)
    t1 = time.time()
    stack = [(i, S0) for i in range(i1 - 1, i0 - 1, -1)]
    leaves = 0; minlb = None; at = None; neg = None; undecided = []
    while stack:
        i, s = stack.pop()
        lb, T0 = D.lower(i, s)
        good = lb > 0 and (goal is None or lb > arb(goal))
        if good:
            leaves += 1
            lo = float(lb.lower())
            if minlb is None or lo < minlb:
                minlb, at = lo, (i, s)
            continue
        if T0 < 0:
            neg = (i, s, T0.str(10))
            break
        if s >= SMAX:
            if lb > 0:                          # positive but the goal is out of reach: accept, record
                leaves += 1
                lo = float(lb.lower())
                if minlb is None or lo < minlb:
                    minlb, at = lo, (i, s)
                continue
            undecided.append((i, s, lb.str(6), T0.str(6)))
            if len(undecided) > 5:
                break
            continue
        stack += [(2*i + 1, s + 1), (2*i, s + 1)]
    return dict(K=K, t=str(t), i0=i0, i1=i1, leaves=leaves, minlb=minlb, at=at, neg=neg,
                undecided=undecided, secs=time.time() - t1)


def certify(jobs, pool):
    """jobs: list of (K, t, chunk). Returns per (K, t) summaries."""
    tasks, meta = [], {}
    for K, t, chunk in jobs:
        W = default_prec(K, t)
        n0 = math.isqrt((8*K + 2)*4**S0 - 1) + 1          # n0 2^-S0 >= sqrt(8K+2), exactly
        goal = None
        if t == 1:
            q0 = Q_data(K, t, W).value(Fraction(0))
            goal = float(q0.mid().str(20, radius=False)) - float(GOAL)
            meta[(K, t)] = dict(W=W, Y=n0/2**S0, Q0=q0.str(12))
        else:
            meta[(K, t)] = dict(W=W, Y=n0/2**S0, Q0=None)
        for i0 in range(0, n0, chunk):
            tasks.append((K, t, W, i0, min(n0, i0 + chunk), goal))
    out = {}
    for r in pool.imap_unordered(cover, tasks, chunksize=1):
        key = (r['K'], Fraction(r['t']))
        o = out.setdefault(key, dict(K=r['K'], t=r['t'], leaves=0, minlb=None, at=None, neg=None,
                                     undecided=[], cpu_secs=0.0, chunks=0))
        o['leaves'] += r['leaves']; o['cpu_secs'] += r['secs']; o['chunks'] += 1
        if r['minlb'] is not None and (o['minlb'] is None or r['minlb'] < o['minlb']):
            o['minlb'], o['at'] = r['minlb'], r['at']
        if r['neg'] and not o['neg']:
            o['neg'] = r['neg']
        o['undecided'] += r['undecided']
    for key, o in out.items():
        m = meta[key]
        o['W'] = m['W']; o['Y'] = m['Y']; o['Q0'] = m['Q0']
        o['certified_positive'] = (o['neg'] is None and not o['undecided'])
        if o['at']:
            i, s = o['at']
            o['at'] = [i/2**s, (i + 1)/2**s]
        print(json.dumps(o), flush=True)
    return out


def cover_min(task):
    """Q(y) > Q(0) for 0 < y on the initial intervals i0 <= i < i1. Q is even, so on [0, 2^-s]
      Q(y) - Q(0) = y^2 (b_2 + sum_{m>=2} b_2m y^(2m-2)) >= y^2 (b_2 - sum_{m>=2} |b_2m| 4^-s(m-1)),
    and elsewhere the Taylor lower bound must exceed Q(0)."""
    K, t, W, i0, i1 = task
    D = get_data(K, t, W)
    ctx.prec = W
    q0 = D.b[0]
    t1 = time.time()
    stack = [(i, S0) for i in range(i1 - 1, i0 - 1, -1)]
    leaves = 0; gap = None; at = None; below = None; undecided = []; s_quad = None
    while stack:
        i, s = stack.pop()
        if i == 0:
            y2 = arb(1)/arb(4)**s
            acc = arb(0); p = arb(1)
            for m in range(2, K + 1):
                p *= y2
                acc += abs(D.b[2*m])*p
            if D.b[2] - acc > 0:
                leaves += 1; s_quad = s
                continue
        else:
            lb, T0 = D.lower(i, s)
            if lb > q0:
                leaves += 1
                g = float((lb - q0).lower())
                if gap is None or g < gap:
                    gap, at = g, (i, s)
                continue
            if T0 < q0:
                below = (i, s, (T0 - q0).str(6))
                break
        if s >= SMAX:
            undecided.append((i, s))
            if len(undecided) > 5:
                break
            continue
        stack += [(2*i + 1, s + 1), (2*i, s + 1)]
    return dict(K=K, leaves=leaves, gap=gap, at=at, below=below, undecided=undecided, s_quad=s_quad,
                b2=D.b[2].str(8) if s_quad is not None else None, secs=time.time() - t1)


def load_float_inputs():
    global GF, VARF
    L = [l.strip() for l in open('fm_zeros_hp.txt') if l.strip() and not l.startswith('#')]
    VARF = float(arb(L[0]).mid().str(20, radius=False))
    GF = [float(arb(x).mid().str(20, radius=False)) for x in L[1:]]


if __name__ == '__main__':
    t0 = time.time()
    mode = sys.argv[1]
    if mode == 'zeros':
        mode_zeros(int(sys.argv[2]), int(sys.argv[3]))
        sys.exit(0)
    ctx.prec = 128
    load_float_inputs()
    if mode == 'table':
        H = json.load(open('fm_hermite.json'))['K']
        rows = [r for r in H if r['K'] >= 2 and 'ratio' in r]
        table = []
        V, G = read_inputs(256)
        with Pool(CPUS) as pool:
            def status(K, t):
                return certify([(K, t, 8)], pool)[(K, t)]
            for r in rows:
                K = r['K']
                lo, hi = (Fraction(x) for x in r['ratio'])          # dyadic rationals, exact as floats
                assert all(q.denominator & (q.denominator - 1) == 0 for q in (lo, hi)), (K, lo, hi)
                a = status(K, Fraction(1))
                # bracket v*/2beta_K: a certified negative value at lo, a complete positive cover at hi
                calls = []
                sl = status(K, lo); calls.append((str(lo), 'neg' if sl['neg'] else 'pos' if sl['certified_positive'] else 'und'))
                step = hi - lo
                while not sl['neg']:                    # the threshold is below lo: move down
                    assert sl['certified_positive'], sl
                    hi, lo = lo, lo - step; step *= 2
                    sl = status(K, lo); calls.append((str(lo), 'neg' if sl['neg'] else 'pos'))
                sh = status(K, hi); calls.append((str(hi), 'pos' if sh['certified_positive'] else 'neg' if sh['neg'] else 'und'))
                step = hi - lo
                while not sh['certified_positive']:     # the threshold is above hi: move up
                    assert sh['neg'], sh
                    lo, sl = hi, sh
                    hi = hi + step; step *= 2
                    sh = status(K, hi); calls.append((str(hi), 'pos' if sh['certified_positive'] else 'neg' if sh['neg'] else 'und'))
                while hi - lo > Fraction(1, 2**18):
                    mid = (lo + hi)/2
                    sm = status(K, mid); calls.append((str(mid), 'pos' if sm['certified_positive'] else 'neg' if sm['neg'] else 'und'))
                    if sm['certified_positive']:
                        hi, sh = mid, sm
                    else:
                        assert sm['neg'], sm
                        lo, sl = mid, sm
                ctx.prec = 256
                two_beta = V - 2*sum((1/(g*g) for g in G[:K]), arb(0))
                i, sd, val = sl['neg']
                table.append(dict(K=K, two_beta=two_beta.str(12), f_K_positive=a['certified_positive'],
                                  minQ_lower=a['minlb'], Q0=a['Q0'], t_lo=str(lo), t_hi=str(hi),
                                  ratio=[float(lo), float(hi)], grid_ratio=r['ratio'],
                                  negative_at_t_lo=dict(y=(2*i + 1)/2**(sd + 1), Q=val), positive_at_t_hi=sh['certified_positive'],
                                  vstar_lower=dec(two_beta*arb(lo.numerator)/lo.denominator, 10, False),
                                  vstar_upper=dec(two_beta*arb(hi.numerator)/hi.denominator, 10, True),
                                  calls=calls, cover_t1=a, cover_t_hi=sh))
                print('K = %d: f_K > 0 %s, min Q >= %.6f (Q(0) = %s); v*/2beta in (%s, %s] = (%.7f, %.7f], grid bracket %s; '
                      'v* in (%s, %s]' % (K, a['certified_positive'], a['minlb'] or float('nan'), a['Q0'], lo, hi, float(lo),
                                          float(hi), r['ratio'], table[-1]['vstar_lower'], table[-1]['vstar_upper']), flush=True)
        ok = all(x['f_K_positive'] and x['positive_at_t_hi'] for x in table)
        json.dump(dict(ok=ok, S0=S0, SMAX=SMAX, rows=table, secs=time.time() - t0), open('fm_qcert_table.json', 'w'), indent=1)
        print('table: all certified %s, %.0f s' % (ok, time.time() - t0))
    elif mode == 'all':
        N = int(sys.argv[2])
        with Pool(CPUS) as pool:
            res = certify([(K, Fraction(1), 10**6) for K in range(N, 0, -1)], pool)
        rows = [res[(K, Fraction(1))] for K in range(1, N + 1)]
        bad = [r['K'] for r in rows if not r['certified_positive']]
        worst = min(rows, key=lambda r: r['minlb'] if r['minlb'] is not None else -1)
        json.dump(dict(N=N, failures=bad, rows=rows, secs=time.time() - t0), open('fm_qcert_all.json', 'w'), indent=1)
        print('K = 1..%d: f_K certified positive on the line for every K: %s; smallest certified lower bound of min Q %.5f at K = %d; %.0f s'
              % (N, not bad, worst['minlb'], worst['K'], time.time() - t0))
    elif mode == 'minq':
        # the minimum of Q_{K, 2 beta_K} on the line is Q(0): Q > Q(0) on (0, Y] by cover_min, and for y >= Y >= sqrt(8K+2)
        # every term e_j v^-j He_2j(y) is positive, so Q(y) >= He_0 = 1 > Q(0); Q is even.
        Ks = []
        for part in sys.argv[2].split(','):
            a, _, b = part.partition('-')
            Ks += list(range(int(a), int(b or a) + 1))
        tasks, meta = [], {}
        for K in Ks:
            W = default_prec(K, Fraction(1))
            n0 = math.isqrt((8*K + 2)*4**S0 - 1) + 1
            q0 = Q_data(K, Fraction(1), W).value(Fraction(0))
            assert q0 < 1, (K, q0)
            meta[K] = dict(W=W, Y=n0/2**S0, Q0=q0.str(15), q0=q0)
            for i0 in range(0, n0, 8):
                tasks.append((K, Fraction(1), W, i0, min(n0, i0 + 8)))
        out = {K: dict(K=K, leaves=0, gap=None, at=None, below=None, undecided=[], s_quad=None, b2=None, cpu_secs=0.0)
               for K in Ks}
        with Pool(CPUS) as pool:
            for r in pool.imap_unordered(cover_min, tasks[::-1], chunksize=1):
                o = out[r['K']]
                o['leaves'] += r['leaves']; o['cpu_secs'] += r['secs']; o['undecided'] += r['undecided']
                if r['gap'] is not None and (o['gap'] is None or r['gap'] < o['gap']):
                    o['gap'], o['at'] = r['gap'], [r['at'][0]/2**r['at'][1], (r['at'][0] + 1)/2**r['at'][1]]
                if r['below'] and not o['below']:
                    o['below'] = r['below']
                if r['s_quad'] is not None:
                    o['s_quad'], o['b2'] = r['s_quad'], r['b2']
        rows = []
        for K in Ks:
            o = out[K]; m = meta[K]
            o['W'], o['Y'], o['Q0'] = m['W'], m['Y'], m['Q0']
            o['min_at_0'] = o['below'] is None and not o['undecided'] and o['s_quad'] is not None
            rows.append(o)
            print(json.dumps(o), flush=True)
        ctx.prec = 256
        dec = [(K1, K2) for K1, K2 in zip(Ks, Ks[1:]) if not meta[K2]['q0'] < meta[K1]['q0']]
        ok = all(o['min_at_0'] for o in rows)
        json.dump(dict(ok=ok, Q0_strictly_decreasing=not dec, not_decreasing=dec, rows=rows, secs=time.time() - t0),
                  open('fm_qcert_minq.json', 'w'), indent=1)
        print('minq: min Q = Q(0) certified for every K: %s; Q(0) strictly decreasing along the list: %s; %.0f s'
              % (ok, not dec, time.time() - t0))
