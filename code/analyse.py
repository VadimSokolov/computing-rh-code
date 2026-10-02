import numpy as np, json, mpmath as mp
from scipy import integrate
from scipy.interpolate import CubicSpline, CubicHermiteSpline
P1 = json.load(open('part1.json'))
Theta = 200
gl = np.array([np.longdouble(s) for s in P1['ref']])          # first 80 zeros, 45 digits
gall = np.concatenate([np.load(f) for f in ['g_1_400.npy','g_401_1000.npy','g_1001_1700.npy']])
near_hp = gl[:80]
mid = gall[(gall>202) & (gall<300)].astype(np.longdouble)
near = np.concatenate([near_hp, mid])
far = gall[gall>=300]
mp.mp.dps = 30
g1701 = float(mp.zetazero(1701).imag)
Nsm = lambda T: float(mp.siegeltheta(T)/mp.pi + 1)
Tstar = float(mp.findroot(lambda T: mp.siegeltheta(T)/mp.pi + 1 - 1700, 2197.5))
assert gall[-1] <= Tstar < g1701, (gall[-1], Tstar, g1701)
rho = lambda g: (0.5*np.log(g/(2*np.pi)) + 1/(48*g*g))/np.pi
PI = np.pi
def far_part(eps, th):
    c = np.linspace(0, Theta, 401)
    v = []
    for t in c:
        s = np.sum(eps/(eps**2+(t-far)**2) + eps/(eps**2+(t+far)**2))/PI
        tail,_ = integrate.quad(lambda g: (eps/(eps**2+(t-g)**2)+eps/(eps**2+(t+g)**2))/PI*rho(g), Tstar, np.inf, epsabs=1e-16, epsrel=1e-12, limit=200)
        v.append(s+tail)
    return CubicSpline(c, v)(th)
res = dict(Tstar=Tstar, g1701=g1701, N_arg=P1['N_arg'], n_sign=P1['n_sign'], maxdiff=P1['maxdiff'], rows=[])
cells = np.concatenate([[0.0], [float((near_hp[k]+near_hp[k+1])/2) for k in range(78)], [Theta]])
Wt = [0.005, 0.01, 0.02]
Wzero = [float(np.sum(np.exp(-t*gall.astype(float)**2))) for t in Wt]
curves = {}
for eps in [1.0,0.5,0.25,0.1,0.05,0.02,0.01,0.005]:
    d = np.load(f'eps_{eps}.npz'); th = d['theta']; F = d['F']; A = d['arg']
    thl = th.astype(np.longdouble)
    R = np.zeros(len(th), dtype=np.longdouble)
    e = np.longdouble(eps)
    for g in near:
        R += e/(e*e+(thl-g)**2) + e/(e*e+(thl+g)**2)
    R = R/np.longdouble(np.pi)
    R = R + far_part(eps, th)
    D = (F/PI) - R.astype(float)
    tv_m = integrate.trapezoid(np.abs(D), th)
    M = A/PI
    Mc = CubicHermiteSpline(th, M, F/PI)(cells)   # M_eps at the cell boundaries, with dM/dtheta = F/pi (linear interpolation erred by 1e-4 at eps = 1)
    cm = np.diff(Mc)
    tv_cells = float(np.sum(np.abs(cm-1)))
    Nstep = np.searchsorted(near_hp.astype(float), th)   # number of zeros <= theta
    l1 = integrate.trapezoid(np.abs(M-Nstep), th)
    l1z = integrate.trapezoid(np.abs(M-Nstep)*2*th, th)
    W = [float(integrate.trapezoid(np.exp(-t*th**2)*F/PI, th)) for t in Wt]
    row = dict(eps=eps, npts=len(th), mass=float(M[-1]), minF=float(F.min()), tv_matched=float(tv_m),
               sup_matched=float(np.abs(D).max()), tv_raw=float(M[-1]+79), tv_cells=tv_cells,
               cell_min=float(cm.min()), cell_max=float(cm.max()), l1_theta=float(l1), l1_z=float(l1z), W=W)
    res['rows'].append(row); print(row, flush=True)
    step = max(1, len(th)//6000)
    curves[f'th_{eps}'] = th[::step]; curves[f'F_{eps}'] = F[::step]; curves[f'M_{eps}'] = M[::step]
    curves[f'D_{eps}'] = D[::step]
res['Wt'] = Wt; res['Wzero'] = Wzero
json.dump(res, open('analysis.json','w'), indent=1)
np.savez('plotcurves.npz', **curves, zeros=near_hp.astype(float))
print('W zeros', Wzero)
