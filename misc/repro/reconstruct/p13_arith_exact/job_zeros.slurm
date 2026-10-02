#!/bin/bash
# Item arith_exact: the first 1700 zero ordinates to 70 digits (zeros_hp.py).
#SBATCH -J p13-arith-zeros
#SBATCH -p normal
#SBATCH --constraint=intel
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH -t 60
#SBATCH -o /scratch/vsokolov/rh_book_repro/agents/p13_arith_exact/slurm-%j.log
cd /scratch/vsokolov/rh_book_repro/agents/p13_arith_exact
echo "[p13_arith_exact] job $SLURM_JOB_ID on $(hostname), $(lscpu | grep 'Model name' | sed 's/  */ /g')" >&2
/scratch/vsokolov/rh_book_repro/venv/bin/python3 -u zeros_hp.py
