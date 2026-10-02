#!/bin/bash
#SBATCH --job-name=rhbook_v1min
#SBATCH --partition=normal
#SBATCH --constraint=intel
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=0-06:00:00
#SBATCH --output=/scratch/vsokolov/rh_book_repro/riesz/v1min/v1min_%j.log
# Certifies that the Thorin density of T_1 is minimal at theta = 0 (riesz_program/code/v1min.py, added for the book).
echo "job $SLURM_JOB_ID on $(hostname) start $(date)"
cd /scratch/vsokolov/rh_book_repro/riesz/v1min
/scratch/vsokolov/rh_book_repro/venv/bin/python3 -u v1min.py
echo "end $(date)"
