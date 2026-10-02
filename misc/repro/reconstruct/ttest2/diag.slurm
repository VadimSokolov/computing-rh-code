#!/bin/bash
# Reconstruction job, not the authors' script: runs the diagnostics floor_check.py (split of the errors of
# Table tab:tilt, ch/ch04.tex) and horner100.py (Horner experiment of Section sec:ch4:balls) on the archived
# code/tilted.py, copied beside them.
#SBATCH -J rhb-ttest2-diag
#SBATCH -p normal
#SBATCH --constraint=intel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 60
#SBATCH -o /scratch/vsokolov/rh_book_repro/agents/reconstruct_ttest2/out/diag_%j.log
cd /scratch/vsokolov/rh_book_repro/agents/reconstruct_ttest2
PY=/scratch/vsokolov/rh_book_repro/venv/bin/python3
echo "job $SLURM_JOB_ID on $(hostname)"
echo "== floor_check.py"; $PY -u floor_check.py
echo "== horner100.py"; $PY -u horner100.py
echo "== done"
