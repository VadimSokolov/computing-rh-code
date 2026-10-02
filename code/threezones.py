# Figure fig:ch19:zones (fig/threezones.pdf) from volterra.json, rsboundary.json and mc2.json: the boundary computed
# from the Volterra equation, the tangent boundary and the tail asymptote over the three zones, and the crossing time errors.
import json, os, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':0.3})
V=json.load(open('volterra.json')); R=json.load(open('rsboundary.json')); M=json.load(open('mc2.json'))
tv=np.array(V['t']); bv=np.array(V['b']); tr=np.array(R['t']); br=np.array(R['b']); m=tr<=0.3
a,c,beta=M['tail']['a'],M['tail']['c'],M['tail']['beta']
fig,ax=plt.subplots(1,2,figsize=(10,3.6))
for (lo,hi),col,lab,x in [((0,0.011),'tab:green','small times',0.0055),((0.011,0.05),'tab:orange','bulk',0.0305),((0.05,0.3),'tab:red','right tail',0.175)]:
    ax[0].axvspan(lo,hi,color=col,alpha=0.12,lw=0)
    ax[0].text(x,-6.15,lab,ha='center',va='bottom',fontsize=8,rotation=90 if hi-lo<0.02 else 0)
ax[0].plot(tv,bv,color='#1f4e79',lw=1.8,label='Volterra equation')
ax[0].plot(tr[m],br[m],'--',color='#c0392b',lw=1.1,label='Roberts Shortland tangent')
tt=np.linspace(0.05,0.3,300)
ax[0].plot(tt,beta*tt+a+c*np.log(tt),':',color='black',lw=1.3,label=r'tail $\beta t+a+\frac{3}{2|\beta|}\log t$')
ax[0].set_xlim(0,0.3); ax[0].set_ylim(-6.3,0.9); ax[0].set_xlabel('t'); ax[0].set_ylabel('boundary b(t)')
ax[0].set_title('Brownian boundary of the clock',fontsize=10); ax[0].legend(loc='upper right',fontsize=8)
g=np.array(M['grid']); Ft=np.array(M['Ft']); n=200000  # paths per boundary in mc2.py
ax[1].axhspan(-1/np.sqrt(n),1/np.sqrt(n),color='gray',alpha=0.18,lw=0,label=r'Monte Carlo noise, $\pm1/\sqrt{n}$')
ax[1].plot(g,np.array(M['tangent']['Fe'])-Ft,color='#c0392b',lw=1.1,label=f"tangent, KS {M['tangent']['ks']:.4f}")
ax[1].plot(g,np.array(M['volterra']['Fe'])-Ft,color='#1f4e79',lw=1.1,label=f"Volterra, KS {M['volterra']['ks']:.4f}")
ax[1].axhline(0,color='black',lw=0.6)
ax[1].set_xlabel('t'); ax[1].set_ylabel('simulated minus clock CDF'); ax[1].set_title(f'Crossing time, {n} paths',fontsize=10)
ax[1].legend(loc='upper right',fontsize=8)
os.makedirs('fig',exist_ok=True); fig.tight_layout(); fig.savefig('fig/threezones.pdf'); plt.close()
