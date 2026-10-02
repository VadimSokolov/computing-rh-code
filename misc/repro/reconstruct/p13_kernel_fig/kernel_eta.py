# kernel_eta.py (round 2 item kernel_fig, not part of the authors' package). A copy of the block
# "# 1 Polya kernel on lines in the strip" of code/plots3.py (its lines 1 to 22), which draws fig/kernel.pdf
# (Figure fig:ch2:kernel of ch/ch02.tex); its output kernel.pdf replaces that figure.
# Only the labels of the right panel change, to follow notation decision 15 of misc/reviews/round2/notation-audit.md
# and item 5 of misc/reviews/round2/notation-decisions.md: each curve is labelled by the height eta = pi/4 - d of its
# line (0.1000, 0.3927, 0.6854, 0.7654) instead of by the distance d below the edge of the strip (0.6854, 0.3927, 0.1,
# 0.02, printed as delta), and the panel title reads |Phi(x+i eta)| instead of |Phi(x+i(pi/4-delta))|.
# The grids, the precision, the series, the colours and the sizes are those of code/plots3.py. Lines 7 and 8 of
# plots3.py (the zeros g and the Poisson kernel P, used only by its later figures) are left out, and the output goes
# to the run folder instead of book/fig/. Lines changed against plots3.py are marked "# changed".
# With the argument "old" the script draws the archived labels instead (a control, kernel_oldlabels.pdf), so that a
# comparison of the three PDFs separates the change of labels from the change of matplotlib version (3.10.8 to 3.11.2).
# Both modes also write the plotted values and the axes geometry to kernel[_oldlabels]_curves.json.
# Run: bash misc/tools/hopper_run.sh -g kernel.pdf kernel_eta.py code/tilted.py   (add "-- old" for the control)
import sys  # changed (mode switch)
import json, numpy as np, mpmath as mp, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from flint import acb, arb, ctx
import tilted
plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':0.25,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':200,'axes.prop_cycle':plt.cycler(color=['#1b4f72','#c0392b','#27ae60','#d68910','#7d3c98','#17a589','#566573'])})
F=''  # changed (was 'book/fig/')
OLD = len(sys.argv) > 1 and sys.argv[1] == 'old'  # changed (mode switch)
name = 'kernel_oldlabels' if OLD else 'kernel'  # changed (mode switch)
# 1 Polya kernel on lines in the strip
ctx.prec=100
fig,ax=plt.subplots(1,2,figsize=(9.5,3.5))
u=np.linspace(-1.2,1.2,601)
Phi=[float(tilted.Phi_c(acb(float(x),0),120).real) for x in u]
ax[0].plot(u,Phi,lw=2); ax[0].fill_between(u,0,Phi,alpha=0.15)
ax[0].set_title('Pólya kernel $\Phi(u)$'); ax[0].set_xlabel('u')
xx=np.linspace(-3,3,601)
curves = {}  # changed (record of the plotted values)
for i,d in enumerate([0.6854,0.3927,0.1,0.02]):
    eta=arb.pi()/4-arb(d)
    v=[abs(complex(tilted.Phi_c(acb(float(x),eta),120))) for x in xx]
    lab = r'$\delta=%s$'%d if OLD else r'$\eta=%.4f$'%(np.pi/4-d)  # changed (label of the decision)
    ax[1].semilogy(xx,np.maximum(v,1e-30),lw=1.4,label=lab)  # changed (label=lab)
    curves[str(d)] = dict(eta=np.pi/4-d, label=lab, v=v)  # changed (record)
ttl = r'$|\Phi(x+i(\pi/4-\delta))|$ inside the strip' if OLD else r'$|\Phi(x+i\eta)|$ inside the strip'  # changed (title of the decision)
ax[1].set_ylim(1e-12,1e5); ax[1].set_title(ttl); ax[1].set_xlabel('x'); ax[1].legend(fontsize=8)  # changed (set_title(ttl))
fig.tight_layout(); fig.savefig(F+name+'.pdf')  # changed (file name; plt.close moved below)
# record of the plotted values and of the axes geometry in PDF points (for pdf_compare.py)
W, H = fig.get_size_inches() * 72
geo = [dict(box=[p.x0 * W, p.y0 * H, p.x1 * W, p.y1 * H], xlim=list(a.get_xlim()), ylim=list(a.get_ylim()), yscale=a.get_yscale(), title=a.get_title(), legend=[t.get_text() for t in a.get_legend().get_texts()] if a.get_legend() else [])
       for a in ax for p in [a.get_position()]]
json.dump(dict(mode='old' if OLD else 'new', matplotlib=matplotlib.__version__, u=list(u), Phi=Phi, x=list(xx), curves=curves, axes=geo), open(name + '_curves.json', 'w'))
plt.close()
for d, c in curves.items():
    v = np.array(c['v'])
    print('%s  d=%-6s eta=%.7f  max|Phi|=%.6g at x=%+.2f  |Phi(i eta)|=%.6g' % (c['label'], d, c['eta'], v.max(), xx[v.argmax()], v[300]))
print('Phi(0) =', Phi[300], ' axes:', json.dumps(geo))
print('matplotlib', matplotlib.__version__, 'wrote', name + '.pdf')
