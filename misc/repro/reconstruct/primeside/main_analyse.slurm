#!/bin/bash
#SBATCH -J rhb-reconstruct_primeside_mainan
#SBATCH -p normal
#SBATCH --constraint=intel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 30
#SBATCH -o /scratch/vsokolov/rh_book_repro/agents/reconstruct_primeside/logs/main_analyse_%j.log
A=/scratch/vsokolov/rh_book_repro/agents/reconstruct_primeside
PY=/scratch/vsokolov/rh_book_repro/venv/bin/python3
echo "job $SLURM_JOB_ID on $(hostname) start $(date)"
cd $A/main && $PY -u main_analyse.py $A/book/results_main.json; echo "rc=$? end $(date)"
