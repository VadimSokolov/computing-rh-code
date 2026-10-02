#!/bin/bash
# Item arith_exact: the patched arith.py at 40, 60 and 80 digits, and the decomposition of the change (decompose.py).
#SBATCH -J p13-arith-run
#SBATCH -p normal
#SBATCH --constraint=intel
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 180
#SBATCH --array=0-3
#SBATCH -o /scratch/vsokolov/rh_book_repro/agents/p13_arith_exact/slurm-%A_%a.log
cd /scratch/vsokolov/rh_book_repro/agents/p13_arith_exact
echo "[p13_arith_exact] job ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} on $(hostname), $(lscpu | grep 'Model name' | sed 's/  */ /g')" >&2
case $SLURM_ARRAY_TASK_ID in
  0) /scratch/vsokolov/rh_book_repro/venv/bin/python3 -u arith.py 40 arith_dps40.json ;;
  1) /scratch/vsokolov/rh_book_repro/venv/bin/python3 -u arith.py 60 arith.json ;;
  2) /scratch/vsokolov/rh_book_repro/venv/bin/python3 -u arith.py 80 arith_dps80.json ;;
  3) /scratch/vsokolov/rh_book_repro/venv/bin/python3 -u decompose.py ;;
esac
