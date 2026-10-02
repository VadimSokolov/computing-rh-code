#!/bin/bash
#SBATCH -J rhb-reconstruct_primeside_main
#SBATCH -p normal
#SBATCH --constraint=intel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 120
#SBATCH --array=0-5
#SBATCH -o /scratch/vsokolov/rh_book_repro/agents/reconstruct_primeside/logs/main_%A_%a.log
# Attempted reconstruction of data/results_main.json: one resolution per array task (main_compute.py).
A=/scratch/vsokolov/rh_book_repro/agents/reconstruct_primeside
PY=/scratch/vsokolov/rh_book_repro/venv/bin/python3
EPS=(0.5 0.2 0.1 0.05 0.02 0.01)
e=${EPS[$SLURM_ARRAY_TASK_ID]}
echo "job $SLURM_JOB_ID task $SLURM_ARRAY_TASK_ID eps $e on $(hostname) start $(date)"; lscpu | grep -E '^Model name'
cd $A/main && $PY -u main_compute.py $e; echo "rc=$? end $(date)"
