#!/bin/bash
# Test run of the proposed code/perturb.py (misc/repro/reconstruct/p13_perturb_fig/code_proposal/) on an Intel node,
# with Matplotlib 3.10.8 (the version that drew the archived fig/perturb.pdf, installed with pip --no-deps in
# ../p13_perturb_fig_mpl/mpl3108), then a comparison of its fig/perturb.pdf with the archived figure.
#SBATCH -J rhb-p13-perturb-code
#SBATCH -p normal
#SBATCH --constraint=intel
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 30
#SBATCH -o slurm-%j.log
cd /scratch/vsokolov/rh_book_repro/agents/p13_perturb_fig_code
echo "[job] $SLURM_JOB_ID on $(hostname), $(grep -m1 'model name' /proc/cpuinfo)"
export PYTHONPATH=/scratch/vsokolov/rh_book_repro/agents/p13_perturb_fig_mpl/mpl3108
export MPLCONFIGDIR=/scratch/vsokolov/rh_book_repro/agents/p13_perturb_fig_code/mplconfig
PY=/scratch/vsokolov/rh_book_repro/venv/bin/python3
start=$SECONDS; $PY -u perturb.py; echo "[job] perturb.py took $((SECONDS-start)) s"
cp /scratch/vsokolov/rh_book_repro/agents/p13_perturb_fig/perturb.json recon.json
$PY -u compare_fig.py old_perturb.pdf recon.json fig/perturb.pdf | grep -E "^==|archived curve"
