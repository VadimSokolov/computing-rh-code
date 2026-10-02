#!/bin/bash
#SBATCH --job-name=rhbook_repro
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --exclude=amd[021-060]
#SBATCH --output=/scratch/vsokolov/rh_book_repro/logs/repro_%j.log
# Reproducibility run for the book "Computing the Riemann Hypothesis": regenerate the zero lists,
# run every script behind the tables in dependency order, and compare with the archived data/*.json.
echo "job $SLURM_JOB_ID on $(hostname) cpus $SLURM_CPUS_PER_TASK start $(date)"
# The base folder is $RH_REPRO_BASE if it is set (a rerun in a fresh folder), else /scratch/vsokolov/rh_book_repro;
# run_all.py and compare.py read the same variable.
cd ${RH_REPRO_BASE:-/scratch/vsokolov/rh_book_repro}
/scratch/vsokolov/rh_book_repro/venv/bin/python3 -u run_all.py
/scratch/vsokolov/rh_book_repro/venv/bin/python3 -u compare.py > /dev/null
echo "end $(date)"
