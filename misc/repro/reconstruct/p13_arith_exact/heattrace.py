# Section 2 of code/plots3.py (item arith_exact): draws fig/heattrace.pdf (Figure fig:ch14:W of ch/ch14.tex) from
# arith.json. The rcParams and the plotting lines are those of code/plots3.py; only the input and output file
# names are arguments (code/plots3.py reads arith.json and writes book/fig/heattrace.pdf).
# Usage: python3 heattrace.py [input json] [output file]   (defaults arith.json, heattrace.pdf; .png also works)
import sys, json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':0.25,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':200,'axes.prop_cycle':plt.cycler(color=['#1b4f72','#c0392b','#27ae60','#d68910','#7d3c98','#17a589','#566573'])})
src = sys.argv[1] if len(sys.argv) > 1 else 'arith.json'
dst = sys.argv[2] if len(sys.argv) > 2 else 'heattrace.pdf'
g=np.concatenate([np.load(f) for f in ['g_1_400.npy','g_401_1000.npy','g_1001_1700.npy']])
# 2 heat trace and arithmetic components
A=json.load(open(src)); t=np.array([r[0] for r in A])
Pi=np.array([float(r[1]) for r in A]); Ar=np.array([float(r[2]) for r in A]); Pr=np.array([float(r[3]) for r in A]); Wz=np.array([float(r[5]) for r in A]); dig=np.array([r[6] for r in A])
tt=np.linspace(0.002,0.3,300); Wt=[np.sum(np.exp(-g**2*x)) for x in tt]
fig,ax=plt.subplots(1,2,figsize=(9.5,3.5))
ax[0].semilogy(tt,Wt,lw=2,label='W(t) from 1700 zeros')
ax[0].semilogy(tt,np.exp(-g[0]**2*tt),'--',lw=1,label=r'$e^{-\gamma_1^2 t}$')
ax[0].semilogy(t,Pi,'o',ms=4,label='polar'); ax[0].semilogy(t,np.abs(Ar),'s',ms=4,label='|archimedean|'); ax[0].semilogy(t[t>0.004],Pr[t>0.004],'^',ms=4,label='prime')
ax[0].set_xlabel('t'); ax[0].legend(fontsize=7); ax[0].set_title('Heat trace and its arithmetic pieces')
ax[1].plot(t,dig,'o-',lw=1.8); ax[1].plot(tt,tt*g[0]**2/np.log(10),'--',color='gray',label=r'$\gamma_1^2t/\log 10$')
ax[1].set_xlabel('t'); ax[1].set_ylabel('decimal digits lost'); ax[1].legend(fontsize=8); ax[1].set_title('Cancellation in the explicit formula')
fig.tight_layout(); fig.savefig(dst); plt.close()
print('wrote', dst, 'from', src)
