# Chapter 13, the zero completed laws f_K of the Polya law, in the book's units.
#   phi_K(theta) = exp(-beta_K theta^2) prod_{k<=K} (1 - theta^2/gamma_k^2),  beta_K = Var X/2 - sum_{k<=K} gamma_k^-2,
#   f_{K,v}(x) = g_v(x) Q_{K,v}(x/sqrt v),  Q_{K,v}(y) = sum_{j<=K} e_j v^-j He_2j(y),
# with e_j the elementary symmetric functions of gamma_1^-2, ..., gamma_K^-2, so that f_K = f_{K, 2 beta_K}.
# Positivity is monotone in v (f_{K,v'} = f_{K,v} * g_{v'-v}), so there is a threshold v*(K) with f_{K,v} >= 0 iff
# v >= v*; f_K is a density iff v*(K) <= 2 beta_K. The Maxwell construction gives v*(K) <= sum_{k<=K} gamma_k^-2.
# Q has positive coefficients and He_2j(y) > 0 beyond the largest zero of He_2K, which is below sqrt(8K+2), so the
# sign of f_{K,v} is decided on 0 <= y <= sqrt(8K+2). The zeros of He_2K are about pi/sqrt(2K) apart, and Q is
# checked on a grid of step 0.01 over that whole range (at least eight points per half oscillation for K <= 800).
# Q is summed by He_{n+1} = y He_n - n He_{n-1} in arb ball arithmetic (python-flint), so every sign is certified at
# its grid point; the precision is raised until no ball contains 0. Between grid points the check is numerical.
# Usage: python3 fm_hermite.py K1 K2 ...   Reads zeros_1700.txt; writes fm_hermite.json.
#        python3 fm_hermite.py all N       Only the scan at v = 2 beta_K, for every K = 1..N; writes fm_hermite_all.json.
import json, os, sys, time
from multiprocessing import Pool
import numpy as np, mpmath as mp
from flint import arb, ctx

mp.mp.dps = 60
xi = lambda s: s*(s - 1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
VAR = mp.nstr(mp.diff(xi, mp.mpf(1)/2, 2)/xi(mp.mpf(1)/2), 55)     # E X^2 = xi''(1/2)/xi(1/2), no hypothesis
GAM = [l.strip() for l in open('zeros_1700.txt') if l.strip()]


class Family:
    """f_{K,v} for v = t * 2 beta_K, t a dyadic rational, at a working precision that rises when needed."""

    def __init__(self, K):
        self.K = K
        S = sum(mp.mpf(s)**-2 for s in GAM[:K]); v = mp.mpf(VAR) - 2*S
        self.prec = int(3.33*(30 + K*max(0.0, float(mp.log10(32*S/v))))) + 64
        self.setprec(self.prec)

    def setprec(self, prec):
        self.prec = prec; ctx.prec = prec
        w = [1/(arb(s)*arb(s)) for s in GAM[:self.K]]
        e = [arb(1)] + [arb(0)]*self.K
        for wi in w:
            for j in range(self.K, 0, -1):
                e[j] += wi*e[j - 1]
        self.e = e; self.S = sum(w, arb(0)); self.v = arb(VAR) - 2*self.S

    def scan(self, t, ys, full=False):
        """Sign of Q_{K, t v} on the grid: (+1, ...) all certified positive, (-1, ...) some certified negative."""
        while True:
            ctx.prec = self.prec
            vt = self.v*arb(float(t))           # t is dyadic, exact as a float
            c = [self.e[j]/vt**j for j in range(self.K + 1)]
            und, neg, qmin, ymin = False, [], None, None
            for y in ys:
                Y = arb(y); H0, H1, s = arb(1), Y, c[0]
                for n in range(1, 2*self.K):
                    H0, H1 = H1, Y*H1 - n*H0
                    if n & 1:
                        s += c[(n + 1) >> 1]*H1
                if s < 0:
                    neg.append(y)
                    if not full:
                        break
                elif not s > 0:
                    und = True
                    break
                if qmin is None or s.mid() < qmin:
                    qmin, ymin = s.mid(), y
            if neg:
                return -1, dict(neg_y=[min(neg), max(neg)], n_neg=len(neg))
            if not und:
                return 1, dict(Qmin=float(qmin), y_at=ymin)
            self.setprec(int(1.5*self.prec))


def run(K):
    t0 = time.time()
    fam = Family(K)
    top = float(np.sqrt(8*K + 2))
    ys = [float(y) for y in np.arange(0, top + 0.01, 0.01)]
    two_beta = float(fam.v.mid()); S = float(fam.S.mid())
    out = dict(K=K, two_beta=two_beta, maxwell=S, maxwell_ok=bool(two_beta >= S), y_top=ys[-1], npts=len(ys))
    sg, info = fam.scan(1, ys, full=True)
    out['at_two_beta'] = dict(sign=sg, **info)
    if sg > 0:
        lo, hi = mp.mpf(1)/16, mp.mpf(1)
        if fam.scan(lo, ys)[0] > 0:
            out['ratio_below'] = float(lo)
        else:
            for _ in range(14):
                mid = (lo + hi)/2
                if fam.scan(mid, ys)[0] > 0:
                    hi = mid
                else:
                    lo = mid
            out['ratio'] = [float(lo), float(hi)]
            out['vstar'] = [float(lo)*two_beta, float(hi)*two_beta]
            out['at_lo'] = fam.scan(lo, ys, full=True)[1]
    else:                                   # f_K itself fails: how much Gaussian is missing
        lo, hi = mp.mpf(1), mp.mpf(2)
        while fam.scan(hi, ys)[0] < 0:
            lo, hi = hi, 2*hi
        for _ in range(14):
            mid = (lo + hi)/2
            if fam.scan(mid, ys)[0] > 0:
                hi = mid
            else:
                lo = mid
        out['ratio'] = [float(lo), float(hi)]
        out['vstar'] = [float(lo)*two_beta, float(hi)*two_beta]
    out['prec_bits'] = fam.prec
    out['secs'] = round(time.time() - t0, 1)
    print(json.dumps(out), flush=True)
    return out


def tv_table(Ks):
    """Total variation distance of f_K to the Polya density (double precision Fourier inversion)."""
    pi = np.pi
    def Phi(t):
        t = np.abs(t); s = 0.0
        for n in range(1, 8):
            s = s + 2*(2*pi**2*n**4*np.exp(4.5*t) - 3*pi*n**2*np.exp(2.5*t))*np.exp(-pi*n**2*np.exp(2*t))
        return s
    X0 = float(xi(mp.mpf(1)/2))
    xs = np.linspace(-2.4, 2.4, 4801); pp = Phi(xs)/X0
    th = np.arange(0, 700, 0.005)
    gf = np.array([float(g) for g in GAM[:200]])
    out = {}
    for K in Ks:
        beta = (float(VAR) - 2*np.sum(gf[:K]**-2.))/2
        logphi = -beta*th**2
        sgn = np.ones_like(th)
        for g in gf[:K]:
            u = 1 - th**2/g**2
            logphi = logphi + np.log(np.abs(u) + 1e-300); sgn = sgn*np.sign(u)
        phi = sgn*np.exp(np.clip(logphi, -700, 700))
        assert abs(phi[-1]) < 1e-30, (K, phi[-1])
        f = np.array([np.trapezoid(phi*np.cos(x*th), th) for x in xs])/pi
        out[K] = dict(TV=float(0.5*np.trapezoid(np.abs(f - pp), xs)), fmin=float(f.min()))
        print('TV', K, out[K], flush=True)
    return out


def at_two_beta(K):
    """Only the scan at v = 2 beta_K: is f_K itself a probability density?"""
    fam = Family(K)
    ys = [float(y) for y in np.arange(0, float(np.sqrt(8*K + 2)) + 0.01, 0.01)]
    sg, info = fam.scan(1, ys)
    return dict(K=K, sign=sg, **info, prec_bits=fam.prec)


if __name__ == '__main__' and sys.argv[1] == 'all':
    # python3 fm_hermite.py all N: the scan at v = 2 beta_K for every K = 1..N; writes fm_hermite_all.json
    N = int(sys.argv[2])
    with Pool(int(os.environ.get('SLURM_CPUS_PER_TASK', '8'))) as pool:
        res = pool.map(at_two_beta, range(N, 0, -1), chunksize=1)
    res = sorted(res, key=lambda r: r['K'])
    bad = [r['K'] for r in res if r['sign'] <= 0]
    worst = min(res, key=lambda r: r.get('Qmin', -1))
    print('K = 1..%d: f_K certified positive on the grid for all K: %s; smallest min Q %.5f at K = %d'
          % (N, not bad, worst.get('Qmin', float('nan')), worst['K']))
    json.dump(dict(VarX=VAR[:20], N=N, failures=bad, K=res), open('fm_hermite_all.json', 'w'), indent=1)
elif __name__ == '__main__':
    Ks = sorted({int(k) for k in sys.argv[1:]}, reverse=True)
    print('Var X =', VAR[:20], ' cpus', os.cpu_count(), flush=True)
    with Pool(min(len(Ks), int(os.environ.get('SLURM_CPUS_PER_TASK', '8')))) as pool:
        res = pool.map(run, Ks, chunksize=1)
    tv = tv_table(sorted(K for K in Ks if K <= 150))
    for r in res:
        if r['K'] in tv:
            r.update(tv[r['K']])
    json.dump(dict(VarX=VAR[:20], K=sorted(res, key=lambda r: r['K'])), open('fm_hermite.json', 'w'), indent=1)
    print('done')
