# Added for the book (not part of the authors' package): redraws Figures wr:fig:wiener (fig/bp1_wiener.pdf) and wr:fig:salem (fig/bp1_salem_riesz.pdf) of Chapter ch:strip, which no archived script draws, from the archived wiener.json, salem_riesz.json, riesz_explicit.json and g_1_400.npy.
# The book's first versions of the two figures, the merged paper's basepoint_one/merged/fig_wiener.pdf and fig_salem_riesz.pdf
# (kept as misc/repro/reconstruct/part6_bp1_figures/old_bp1_wiener.pdf and old_bp1_salem_riesz.pdf), were drawn by the authors
# with a script that is not in the archive. This script draws the same panels from the same data, in the style of the authors'
# figure scripts (xi/figs.py, scale_mixture/figs2.py, riesz_program/code/figs_riesz.py) and with the axis ranges, line widths
# and legend positions of those files, and writes them in the notation of the book:
#   fig_wiener.pdf       alpha = 1/2 and alpha = 1 where the authors' legend has 0.5 and 1.0, and an italic N;
#   fig_salem_riesz.pdf  theta for the ordinate, which the authors' panel calls t (the book writes |eta(sigma+i theta)|),
#                        sigma = 1/2 and sigma = 1 where the legend has 0.5 and 1.0, Mobius spelt with its umlaut (the
#                        authors' legend prints the TeX escape M\"obius literally), an italic x, and the title of the right
#                        panel states Riesz's criterion as the book does, R(x) = O(x^(1/4+delta)), instead of
#                        "R(x)/x^(1/4) bounded".
# No JSON file stores the explicit formula curve of the right panel of fig_salem_riesz.pdf, so it is evaluated here exactly as
# in riesz_program/code/figs_riesz.py: with rho_n = 1/2 + i gamma_n and gamma_n the first 60 ordinates in g_1_400.npy,
#   R(x)/x^(1/4) = sum_{n<=60} 2 Re[Gamma(1-rho_n/2)/(2 zeta'(rho_n)) x^(i gamma_n/2)] + x^(-1/4) sum_{k<=3} k!/(2 zeta'(-2k)) x^(-k),
# with the trivial coefficients k!/(2 zeta'(-2k)) read from riesz_explicit.json, on 3000 points of [1e3, 1e10]. As a check the
# script recomputes the three numbers maxdev, signal and amp1 of riesz_explicit.json at the 36 points of salem_riesz.json in
# [1e3, 1e10], and it prints the values that Sections wr:sec:wiener and wr:sec:salem read off the data of the two figures.
# Run from basepoint_one/wiener:  python3 figs_wiener.py  -> fig_wiener.pdf, fig_salem_riesz.pdf
# (the book's copies are fig/bp1_wiener.pdf and fig/bp1_salem_riesz.pdf). Requirements: numpy, matplotlib, mpmath.
import json, numpy as np, mpmath as mp, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':0.25,'axes.spines.top':False,'axes.spines.right':False})
cols={'0.5':'#27ae60','0.6':'#c0392b','0.75':'#d68910','1.0':'#1b4f72'}
lab={'0.5':'1/2','0.6':'0.6','0.75':'0.75','1.0':'1'}
W=json.load(open('wiener.json')); S=json.load(open('salem_riesz.json')); E=json.load(open('riesz_explicit.json'))

# Figure wr:fig:wiener. Left: best L^1 error against N, with the barrier |g^(gamma_1)| at the centre (statement (2) of
# Proposition wr:prop:wiener). Right: |h^(gamma_1)| of the best approximant against the target g^(gamma_1).
fig,ax=plt.subplots(1,2,figsize=(9.8,3.6))
print('Figure wr:fig:wiener (wiener.json): barrier |g^(gamma_1)| = exp(-kappa^2 gamma_1^2/2) = %.5f'%W['lb'])
for a in ['0.5','0.6','0.75','1.0']:
    rows=sorted([r for r in W['rows'] if r['alpha']==float(a)],key=lambda r:r['N'])
    N=[r['N'] for r in rows]
    ax[0].plot(N,[r['err'] for r in rows],'o-',color=cols[a],label=r'$\alpha=%s$'%lab[a])
    ax[1].plot(N,[r['hhat_g1'] for r in rows],'o-',color=cols[a],label=r'$\alpha=%s$'%lab[a])
    print('  alpha %-4s N %s  L1 error %s  |h^(gamma_1)| %s  max coefficient sum %.2e'%(a,N,' '.join('%.4f'%r['err'] for r in rows),
          ' '.join('%.3g'%r['hhat_g1'] for r in rows),max(r['coefnorm'] for r in rows)))
ax[0].axhline(W['lb'],color='k',ls=':',lw=1.0,label=r'barrier $|\hat g(\gamma_1)|$ at the centre')
ax[0].set_ylim(0,0.5); ax[0].set_xlabel(r'number of translates $N$'); ax[0].set_ylabel(r'best $L^1$ error')
ax[0].set_title(r'Approximating a Gaussian by translates of $\Phi_\alpha$',fontsize=10); ax[0].legend(fontsize=7,loc='upper right')
ax[1].axhline(W['lb'],color='k',ls=':',lw=1.0,label=r'$\hat g(\gamma_1)$')
ax[1].set_xlabel(r'number of translates $N$'); ax[1].set_ylabel(r'$|\hat h(\gamma_1)|$')
ax[1].set_title('The approximant at the first ordinate',fontsize=10); ax[1].legend(fontsize=7,loc='upper right')
fig.tight_layout(); fig.savefig('fig_wiener.pdf')
print('  axis limits: left x %s y %s, right x %s y %s'%(ax[0].get_xlim(),ax[0].get_ylim(),ax[1].get_xlim(),ax[1].get_ylim()))
plt.close(fig)

# Figure wr:fig:salem. Left: |eta(sigma+i theta)| = |k^_sigma(theta)/Gamma(sigma+i theta)| on the grid of salem_riesz.py.
# Right: R(x)/x^(1/4) from the Mobius series (salem_riesz.json) and from the explicit formula over 60 zeros.
fig,ax=plt.subplots(1,2,figsize=(10,3.6))
print('Figure wr:fig:salem, left panel (salem_riesz.json): local minima of |eta(sigma+i theta)| on (0,60]')
for a in ['0.5','0.6','0.75','1.0']:
    d=S['salem'][a]; t=np.array(d['t']); e=np.array(d['eta'])
    ax[0].semilogy(t,e,color=cols[a],lw=1.0,label=r'$\sigma=%s$'%lab[a])
    i=np.where((e[1:-1]<e[:-2])&(e[1:-1]<e[2:]))[0]+1; i=i[np.argsort(e[i])]
    print('  sigma %-4s grid %d points on [%.2f, %.2f]; %d local minima, the deepest (theta, |eta|): %s'%(a,len(t),t[0],t[-1],len(i),
          ', '.join('(%.3f, %.3g)'%(t[j],e[j]) for j in i[:6])))
ax[0].set_ylim(1e-3,5); ax[0].set_xlabel(r'$\theta$'); ax[0].legend(fontsize=7,loc='lower right')
ax[0].set_title(r'Salem kernel: $|\eta(\sigma+i\theta)|$ on the lines $\mathrm{Re}\,s=\sigma$',fontsize=9)
print('  zeros of 1-2^(1-s) on sigma = 1: theta = 2 pi k/log 2 =',', '.join('%.4f'%(2*np.pi*k/np.log(2)) for k in (1,2,3,4)))
mp.mp.dps=20
g=np.load('g_1_400.npy')[:60]
coef=[complex(mp.gamma(1-mp.mpc(0.5,t)/2)/(2*mp.zeta(mp.mpc(0.5,t),derivative=1))) for t in g]
triv=E['triv']
pred=lambda x: sum(2*(c*x**(0.5j*t)).real for c,t in zip(coef,g))+sum(r*x**(-k) for k,r in triv)/x**0.25
xs=np.array(S['riesz']['x']); R=np.array(S['riesz']['R']); m=(xs>=1e3)&(xs<=1e10)
xf=np.geomspace(1e3,1e10,3000)
ax[1].semilogx(xf,[pred(x) for x in xf],color='#c0392b',lw=0.8,label='explicit formula over 60 zeros')
ax[1].semilogx(xs[m],R[m]/xs[m]**0.25,'o',ms=3,color='#1b4f72',label='$R(x)/x^{1/4}$ from the Möbius series')
ax[1].set_xlabel(r'$x$'); ax[1].legend(fontsize=7,loc='lower right')
ax[1].set_title(r"Riesz's criterion: $R(x)=O(x^{1/4+\delta})$",fontsize=9)
fig.tight_layout(); fig.savefig('fig_salem_riesz.pdf')
print('  axis limits: left x %s y %s, right x %s y %s'%(ax[0].get_xlim(),ax[0].get_ylim(),ax[1].get_xlim(),ax[1].get_ylim()))
plt.close(fig)
ratio=R[m]/xs[m]**0.25; pr=np.array([pred(x) for x in xs[m]]); dev=np.abs(ratio-pr)
chk=dict(maxdev=float(dev.max()),signal=float(np.abs(ratio).max()),amp1=float(2*abs(coef[0])))
print('Figure wr:fig:salem, right panel: %d points of salem_riesz.json in [1e3, 1e10], x from %.4g to %.4g'%(m.sum(),xs[m][0],xs[m][-1]))
for k in ['maxdev','signal','amp1']: print('  %-6s recomputed %.16g  riesz_explicit.json %.16g  rel diff %.1e'%(k,chk[k],E[k],abs(chk[k]/E[k]-1)))
tr=[float(mp.gamma(1+k)/(2*mp.zeta(-2*k,derivative=1))) for k in (1,2,3)]
print('  trivial coefficients k!/(2 zeta\'(-2k)) recomputed %s, riesz_explicit.json %s'%(tr,[r for k,r in triv]))
print('wrote fig_wiener.pdf and fig_salem_riesz.pdf')
