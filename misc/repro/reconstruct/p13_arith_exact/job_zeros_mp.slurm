#!/bin/bash
# Item arith_exact: check zeros_1700.txt against mpmath.zetazero at 75 digits, 50 zeros per task (zeros_mp_check.py).
#SBATCH -J p13-arith-zmp
#SBATCH -p normal
#SBATCH --constraint=intel
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -t 120
#SBATCH --array=0-33
#SBATCH -o /scratch/vsokolov/rh_book_repro/agents/p13_arith_exact/zmp/slurm-%A_%a.log
cd /scratch/vsokolov/rh_book_repro/agents/p13_arith_exact
echo "[p13_arith_exact] job ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} on $(hostname), $(lscpu | grep 'Model name' | sed 's/  */ /g')" >&2
/scratch/vsokolov/rh_book_repro/venv/bin/python3 -u zeros_mp_check.py $SLURM_ARRAY_TASK_ID
