# Comparison tool for the item ag_figures (not part of the authors' package): the final run, in one job, of
# riesz_program/code/figs_almost.py, followed by compare_figs.py and render_png.py on the figures it draws.
import subprocess, sys
for s in ('figs_almost.py', 'compare_figs.py', 'render_png.py'):
    print('==== running', s, flush=True)
    r = subprocess.run([sys.executable, '-u', s])
    if r.returncode:
        sys.exit(r.returncode)
