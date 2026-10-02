# Helper of the reconstruction in misc/repro/reconstruct/p13_perturb_fig/ (not the authors' code): section 8 of
# code/plots3.py (lines 92 to 103), verbatim apart from the output name and the imports it needs (lines 1, 2, 5, 7 of
# plots3.py), run to show what it draws in place of the archived fig/perturb.pdf. Writes plots3_section8.pdf and
# prints how many of the 300 values of |W_F-W|/W are exactly 0 in double precision.
import json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':0.25,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':200,'axes.prop_cycle':plt.cycler(color=['#1b4f72','#c0392b','#27ae60','#d68910','#7d3c98','#17a589','#566573'])})
g=np.concatenate([np.load(f) for f in ['g_1_400.npy','g_401_1000.npy','g_1001_1700.npy']])
# 8 perturbation invisibility (Theorem D type)
gp=g[(g>1000)&(g<2000)]; k=int(np.argmin(np.diff(gp))); A_,B_=gp[k],gp[k+1]; gg=(A_+B_)/2
tt=np.logspace(-3,-0.5,300)
W=np.array([np.sum(np.exp(-g**2*x)) for x in tt])
a=gg**2-1/16-1j*gg/2
dW=np.array([2*np.real(np.exp(-a*x))-np.exp(-A_**2*x)-np.exp(-B_**2*x) for x in tt])
fig,ax=plt.subplots(figsize=(6.5,3.5))
ax.loglog(tt,np.abs(dW)/W+1e-300,lw=1.6)
ax.set_xlabel('t'); ax.set_ylabel(r'$|W_F-W|/W$'); ax.set_ylim(1e-40,1)
ax.set_title('Moving a pair near height %.0f off the line changes the heat trace by'%gg+'\nat most this relative amount',fontsize=9)
fig.tight_layout(); fig.savefig('plots3_section8.pdf'); plt.close()
print('pair',A_,B_)
print('zero values of dW:', int(np.sum(dW==0)), 'of', len(dW), '; max |dW|/W = %g' % np.max(np.abs(dW)/W))
