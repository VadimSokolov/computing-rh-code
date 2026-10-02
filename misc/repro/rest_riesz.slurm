#!/bin/bash
#SBATCH --job-name=rhbook_rzrest
#SBATCH --partition=normal
#SBATCH --constraint=intel
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --output=/scratch/vsokolov/rh_book_repro/riesz/logs/rest_%j.log
# Riesz package reproducibility: the figure scripts,
# any task the array left undone, then the comparison with the archived outputs.
echo "job $SLURM_JOB_ID on $(hostname) start $(date)"
cd /scratch/vsokolov/rh_book_repro/riesz
/scratch/vsokolov/rh_book_repro/venv/bin/python3 -u run_riesz.py rest
/scratch/vsokolov/rh_book_repro/venv/bin/python3 -u compare6.py /scratch/vsokolov/rh_book_repro/riesz > /dev/null
echo "end $(date)"
