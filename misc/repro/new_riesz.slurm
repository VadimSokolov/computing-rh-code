#!/bin/bash
#SBATCH --job-name=rhbook_rznew
#SBATCH --partition=normal
#SBATCH --constraint=intel
#SBATCH --array=9-10
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=0-12:00:00
#SBATCH --output=/scratch/vsokolov/rh_book_repro/riesz/logs/new_%A_%a.log
# Adds scripts of a revised Riesz package to an existing Riesz run without redoing the other tasks:
# copies them (and their archived outputs, for comparison) into the run tree and runs the given groups of run_riesz.py.
echo "task $SLURM_ARRAY_TASK_ID of job $SLURM_ARRAY_JOB_ID on $(hostname) start $(date)"
cd /scratch/vsokolov/rh_book_repro/riesz
for f in killing.py killing2.py ks.py; do cp src/code/$f run/code/; done
for f in killing.json ks.json; do cp src/data/$f ref/code/; done
/scratch/vsokolov/rh_book_repro/venv/bin/python3 -u run_riesz.py group $SLURM_ARRAY_TASK_ID
echo "end $(date)"
