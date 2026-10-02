#!/bin/bash
# Final run of the reconstruction of fig/perturb.pdf (misc/repro/reconstruct/p13_perturb_fig/), on an Intel node:
# perturb.py computes the curve and draws the figures with Matplotlib 3.10.8 (installed with pip --no-deps in
# ../p13_perturb_fig_mpl/mpl3108, the version that drew the archived figure), then compare_fig.py compares the
# curves with the archived figure, old_perturb.pdf.
#SBATCH -J rhb-p13-perturb
#SBATCH -p normal
#SBATCH --constraint=intel
#SBATCH -c 16
#SBATCH --mem=16G
#SBATCH -t 60
#SBATCH -o slurm-%j.log
cd /scratch/vsokolov/rh_book_repro/agents/p13_perturb_fig
echo "[job] $SLURM_JOB_ID on $(hostname), $(grep -m1 'model name' /proc/cpuinfo)"
export PYTHONPATH=/scratch/vsokolov/rh_book_repro/agents/p13_perturb_fig_mpl/mpl3108
export MPLCONFIGDIR=/scratch/vsokolov/rh_book_repro/agents/p13_perturb_fig/mplconfig
PY=/scratch/vsokolov/rh_book_repro/venv/bin/python3
$PY -u perturb.py 16
$PY -u compare_fig.py old_perturb.pdf perturb.json perturb_w1700.pdf perturb_exact.pdf
