#!/bin/bash
# Reconstruction job, not the authors' script: parameter search for the step h behind Table tab:tilt
# (ch/ch04.tex). One array task per candidate h; each runs the archived code/ttest2.py unchanged
# through search_task.py. Submit with: sbatch --array=0-115 search.slurm
#SBATCH -J rhb-ttest2-search
#SBATCH -p normal
#SBATCH --constraint=intel
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 30
#SBATCH -o /scratch/vsokolov/rh_book_repro/agents/reconstruct_ttest2/out/search_%A_%a.log
cd /scratch/vsokolov/rh_book_repro/agents/reconstruct_ttest2
/scratch/vsokolov/rh_book_repro/venv/bin/python3 -u search_task.py $SLURM_ARRAY_TASK_ID
