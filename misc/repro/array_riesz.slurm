#!/bin/bash
#SBATCH --job-name=rhbook_rzarr
#SBATCH --partition=normal
#SBATCH --constraint=intel
#SBATCH --array=0-10
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --output=/scratch/vsokolov/rh_book_repro/riesz/logs/array_%A_%a.log
# Riesz package reproducibility: one independent chain of scripts per array task (see GROUPS in run_riesz.py).
echo "task $SLURM_ARRAY_TASK_ID of job $SLURM_ARRAY_JOB_ID on $(hostname) start $(date)"
cd /scratch/vsokolov/rh_book_repro/riesz
/scratch/vsokolov/rh_book_repro/venv/bin/python3 -u run_riesz.py group $SLURM_ARRAY_TASK_ID
echo "end $(date)"
