import numpy as np, json, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from flint import acb
from polya_flint import setup, xi_pair
plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':0.3,'savefig.dpi':200})
C = np.load('plotcurves.npz'); A = json.load(open('analysis.json')); Z = C['zeros']
cm = plt.cm.plasma
E = [1.0,0.5,0.25,0.1,0.05,0.02,0.01,0.005]
# Fig 1: Thorin density rho_eps(theta)
fig, ax = plt.subplots(figsize=(7.5,3.6))
for i,e in enumerate([0.25,0.1,0.05]):
    th=C[f'th_{e}']; F=C[f'F_{e}']; m=(th>10)&(th<56)
    ax.plot(th[m], F[m]/np.pi, color=cm(0.15+0.3*i), lw=1.1, label=fr'$\varepsilon={e}$')
for g in Z[Z<56]: ax.axvline(g, color='teal', ls='--', lw=0.7, alpha=0.8)
ax.set_yscale('log'); ax.set_xlabel(r'$\theta=\sqrt{z}$'); ax.set_ylabel(r'$\rho_\varepsilon(\theta)$')
ax.set_title('Pólya computed Thorin density at basepoint 1/2+ε (dashed: zero ordinates)')
ax.legend(loc='lower right'); fig.tight_layout(); fig.savefig('fig/density.pdf'); plt.close()
# Fig 2: cumulative Thorin mass vs RH counting function
fig, ax = plt.subplots(1,2,figsize=(8,3.4))
tt=np.linspace(0,200,20001); N=np.searchsorted(Z,tt)
for a,(lo,hi) in zip(ax,[(0,200),(28,52)]):
    a.step(tt, N, where='post', color='black', lw=1.2, label='RH atoms')
    for i,e in enumerate([1.0,0.25,0.05]):
        th=C[f'th_{e}']; M=C[f'M_{e}']; m=(th>=lo)&(th<=hi)
        a.plot(th[m], M[m], color=plt.cm.viridis(0.1+0.35*i), lw=1.1, label=fr'$\varepsilon={e}$')
    a.set_xlim(lo,hi); a.set_xlabel(r'$\theta$')
ax[0].set_ylabel(r'$\mu_\varepsilon((0,\theta^2])$'); ax[1].set_ylim(3,10); ax[0].legend(fontsize=8)
fig.suptitle('Cumulative Thorin mass against the counting measure of the squared ordinates')
fig.tight_layout(); fig.savefig('fig/cumulative.pdf'); plt.close()
# Fig 3: z chart density
fig, ax = plt.subplots(figsize=(7.5,3.4))
for i,e in enumerate([0.5,0.1,0.02]):
    th=C[f'th_{e}']; F=C[f'F_{e}']; m=(th>9)&(th<60)
    ax.plot(th[m]**2, F[m]/(2*np.pi*th[m]), color=plt.cm.cividis(0.1+0.4*i), lw=1, label=fr'$\varepsilon={e}$')
for g in Z[Z<60]: ax.axvline(g**2, color='crimson', ls=':', lw=0.7)
ax.set_yscale('log'); ax.set_xlabel(r'$z$'); ax.set_ylabel(r'$u_\varepsilon(z)$'); ax.legend()
ax.set_title(r'Thorin density in the clock variable $z$ (dotted: $\gamma^2$)')
fig.tight_layout(); fig.savefig('fig/zdensity.pdf'); plt.close()
# Fig 4: rates
R = A['rows']; eps=np.array([r['eps'] for r in R])
fig, ax = plt.subplots(figsize=(6.5,4))
series = [('tv_cells','cell TV',None),('l1_theta',r'$L^1$ of cumulative ($\theta$)',None),
          ('mass','mass defect',lambda r: r['mass']-79),('tv_matched','matched TV',None),
          ('W','heat trace error, t=0.01',lambda r: abs(r['W'][1]-A['Wzero'][1]))]
for i,(k,lab,fn) in enumerate(series):
    y=[fn(r) if fn else r[k] for r in R]
    ax.loglog(eps,y,'o-',color=plt.cm.tab10(i),label=lab,ms=4)
ax.set_xlabel(r'$\varepsilon$'); ax.set_ylabel('distance'); ax.legend(fontsize=8)
ax.set_title('Distances to the RH measure on $z\\leq 200^2$')
fig.tight_layout(); fig.savefig('fig/rates.pdf'); plt.close()
# Fig 5: matched difference
fig, ax = plt.subplots(figsize=(7.5,3))
for i,e in enumerate([0.5,0.05,0.005]):
    th=C[f'th_{e}']; D=C[f'D_{e}']
    ax.semilogy(th, np.abs(D)+1e-18, color=plt.cm.coolwarm(0.1+0.4*i), lw=0.9, label=fr'$\varepsilon={e}$')
ax.set_xlabel(r'$\theta$'); ax.set_ylabel(r'$|\rho_\varepsilon-P_\varepsilon*\widetilde{U}|$'); ax.legend(fontsize=8)
ax.set_title('Pólya density minus the Poisson smoothed model of the RH atoms')
fig.tight_layout(); fig.savefig('fig/matched.pdf'); plt.close()
# Fig 6: sensitivity
P=lambda a,x: a/(np.pi*(a*a+x*x))
x=np.concatenate([-np.logspace(-6,5,40001)[::-1],np.logspace(-6,5,40001)])
r=np.logspace(-2,1.5,120)
tvd=[]; tvo=[]; neg=[]
for q in r:           # q = delta/eps, eps = 1
    f=0.5*(P(1+q,x)+P(1-q,x))-P(1,x); tvd.append(np.trapezoid(np.abs(f),x))
    g=P(1+q,x)+P(1-q,x); neg.append(np.trapezoid(np.clip(-0.5*g,0,None),x))
    tvo.append(2.0 if q<=1 else (8/np.pi)*np.arctan(np.sqrt((q+1)/(q-1)))-2)
fig, ax = plt.subplots(figsize=(6.5,3.8))
ax.loglog(r,tvd,color='darkorange',lw=1.6,label='zero displaced (same count)')
ax.loglog(r,0.4135*r**2,'--',color='darkorange',lw=0.8,label=r'$0.4135(\delta/\varepsilon)^2$')
ax.loglog(r,tvo,color='royalblue',lw=1.6,label='zero omitted from the list')
ax.loglog(r,np.maximum(neg,1e-12),color='seagreen',lw=1.6,label='negative mass of the pair')
ax.set_xlabel(r'$\delta/\varepsilon$'); ax.set_ylabel('total variation'); ax.set_ylim(1e-4,4); ax.legend(fontsize=8)
ax.set_title('Signature of one off line pair at resolution ε')
fig.tight_layout(); fig.savefig('fig/sensitivity.pdf'); plt.close()
json.dump(dict(r=list(r),tvd=tvd,tvo=tvo,neg=neg), open('sens2.json','w'))
# Fig 7: weak convergence illustration with a hypothetical displaced fifth zero
g5=Z[4]; d=0.1; th=np.linspace(g5-3,g5+3,6001)
fig, ax = plt.subplots(figsize=(7.5,3.6))
for i,e in enumerate([0.3,0.1,0.05,0.02]):
    base=sum(P(e,th-g)+P(e,th+g) for g in Z if g!=g5)
    true=base+P(e,th-g5)+P(e,th+g5)
    pert=base+0.5*(P(e+d,th-g5)+P(e-d,th-g5))+P(e,th+g5)
    ax.plot(th,pert,color=plt.cm.magma(0.15+0.22*i),lw=1.2,label=fr'off line, $\varepsilon={e}$')
    if e==0.3: ax.plot(th,true,'k:',lw=1,label=r'on line, $\varepsilon=0.3$')
ax.axhline(0,color='gray',lw=0.8); ax.set_ylim(-3,6)
ax.set_xlabel(r'$\theta$'); ax.set_ylabel(r'$\rho_\varepsilon$')
ax.set_title(r'Fifth zero moved to $1/2\pm 0.1+i\gamma_5$: the dip vanishes as $\varepsilon\to0$')
ax.legend(fontsize=7,ncol=2); fig.tight_layout(); fig.savefig('fig/weaklimit.pdf'); plt.close()
# Fig 8: HCM derivative signs
H=json.load(open('hcm.json'))
fig, ax = plt.subplots(1,2,figsize=(8,3.4),sharey=True)
for a,al in zip(ax,['0.5','1.0']):
    for j,row in enumerate([rw for rw in H[al] if rw['w'] in (3,20)]):
        a.semilogy(range(1,7),row['d'],'o-',color=plt.cm.tab10(j),ms=4,label=f"u={row['u']}, w={row['w']}")
    a.set_title(fr'basepoint $\alpha={al}$'); a.set_xlabel('k')
ax[0].set_ylabel(r'$(-1)^k \partial_w^k H_u(w)/H_u(w)$'); ax[1].legend(fontsize=7)
fig.suptitle('HCM test: alternating derivatives in $w=v+1/v$')
fig.tight_layout(); fig.savefig('fig/hcm.pdf'); plt.close()
# mass derivative check and prime table
S=setup(420,200,3)
a,b=xi_pair(S,acb(0,200)); XiprimeOverXi=float((-(b.imag))/a.real)  # d/dtheta Xi = -Im xi'
out=dict(XipXi200=XiprimeOverXi, pred_coeff=-XiprimeOverXi/np.pi)
# first order cell constant of ch07 and ch15: Xi'/Xi at the cell boundaries b_0=0, b_k=(g_k+g_{k+1})/2, b_79=200
bc=[0.0]+[float((Z[k]+Z[k+1])/2) for k in range(78)]+[200.0]
Lp=[]
for t in bc:
    xa,xb=xi_pair(S,acb(0,t)); Lp.append(float((-(xb.imag))/xa.real))
leaks=[-(Lp[k+1]-Lp[k])/np.pi for k in range(79)]
out['kappa_cell']=float(np.sum(np.abs(leaks))); out['cell_leaks']=leaks
pr=json.load(open('results_prime.json')); new=[]
for q in pr:
    a,b=xi_pair(S,acb(q['alpha']-0.5,q['b'])); q['polya']=float((b/a).real); new.append(q)
out['prime']=new
json.dump(out,open('extra.json','w'),indent=1); print(out['XipXi200'], out['pred_coeff'], out['kappa_cell'])
