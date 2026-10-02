# dbnplot.py: Figure fig:ch21:dbn of Section ch:dbn (ch/ch21.tex), added for the book in September 2026. The archive held
# the figure but no script, and the archived figure plotted the midpoints of the grid cells holding the sign changes.
# This script reads dbn.json, the output of dbn.py, and plots the zeros refined by Newton's method that
# Table tab:ch21:traj prints (key z of each entry of traj), with the collision time t_* = -0.2880 that the text quotes
# (key collision, found by bisection), in the layout of the archived figure.
# Run from the folder of dbn.json (data/, or the run folder after dbn.py); writes fig/dbn.pdf below it, creating fig/:
#   python3 dbnplot.py
import json, os, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':0.3})
D=json.load(open('dbn.json')); tc=D['collision']
fig,ax=plt.subplots(figsize=(6,3.4))
for e in D['traj']:
    ax.plot([e['t']]*len(e['z']),e['z'],'o',ms=4,color='navy' if e['t']>=0 else 'crimson')
ax.axvline(tc,color='gray',ls='--',lw=0.8,label=('collision t = %.3f'%tc).replace('-','\u2212'))
ax.set_xlabel('t (de Bruijn Newman time)'); ax.set_ylabel(r'zeros $\theta$'); ax.legend(fontsize=8,loc='center right')
ax.set_title('The closest pair below height 200 under the heat flow')
os.makedirs('fig',exist_ok=True); fig.tight_layout(); fig.savefig('fig/dbn.pdf'); print('ok', tc)
