import numpy as np, json
PI=np.pi
gam = 32.93506158773918969  # fifth ordinate
def K(e,x): return e/(e**2+x**2)
def corr(eps, b, d):
    """(1/pi) Re d/dw log[(1+w^2/a)(1+w^2/abar)/(1+w^2/gam^2)] at w = eps + i b:
    the pair at gam is replaced by the off line quadruple 1/2 +- d + i(+-gam)."""
    w = eps + 1j*b
    a = -(d+1j*gam)**2
    r = 2*w/(a+w**2) + 2*w/(np.conj(a)+w**2) - 2*w/(gam**2+w**2)
    return r.real/PI
out=[]
for d in [0.2, 0.05]:
    for eps in [1.0,0.5,0.3,0.15,0.1,0.07,0.04,0.03,0.02,0.01,0.005,0.002,0.001]:
        db = min(eps,d)/40
        b = np.arange(gam-8, gam+8, db)
        # perturbed minus list measure with gam counted twice: only the local terms differ
        diff = corr(eps,b,d) - (K(eps,b-gam)+K(eps,b+gam))/PI
        tv = np.trapezoid(np.abs(diff), b)
        # negative part of the perturbed density: local model = RH density + corr
        # (RH density near gam is dominated by the gam kernel; include it exactly)
        loc = (K(eps,b-gam)+K(eps,b+gam))/PI + corr(eps,b,d)
        neg = np.trapezoid(np.clip(-loc,0,None), b)
        out.append(dict(d=d,eps=eps,tv=float(tv),neg=float(neg)))
        print(d,eps,tv,neg)
json.dump(out, open('results_sens.json','w'), indent=1)
