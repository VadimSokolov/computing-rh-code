#!/bin/bash
#SBATCH --job-name=rhbook_ks_ak
#SBATCH --partition=normal
#SBATCH --constraint=intel
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=0-01:00:00
#SBATCH --output=/scratch/vsokolov/rh_book_repro/riesz/ks_ak/ks_ak_%j.log
# Recomputes the arithmetic factors that ks.py hard-codes (misc/repro/ks_ak.py).
echo "job $SLURM_JOB_ID on $(hostname) start $(date)"
cd /scratch/vsokolov/rh_book_repro/riesz/ks_ak
/scratch/vsokolov/rh_book_repro/venv/bin/python3 -u ks_ak.py
echo "end $(date)"
