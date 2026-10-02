#!/bin/bash
#SBATCH --job-name=rhbook_p6arr
#SBATCH --partition=normal
#SBATCH --constraint=intel
#SBATCH --array=0-17
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --output=/scratch/vsokolov/rh_book_repro/part6/logs/array6_%A_%a.log
# Part VI reproducibility: one independent chain of scripts per array task (see GROUPS in run6.py).
echo "task $SLURM_ARRAY_TASK_ID of job $SLURM_ARRAY_JOB_ID on $(hostname) start $(date)"
cd /scratch/vsokolov/rh_book_repro/part6
/scratch/vsokolov/rh_book_repro/venv/bin/python3 -u run6.py group $SLURM_ARRAY_TASK_ID
echo "end $(date)"
