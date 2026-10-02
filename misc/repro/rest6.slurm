#!/bin/bash
#SBATCH --job-name=rhbook_p6rest
#SBATCH --partition=normal
#SBATCH --constraint=intel
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --output=/scratch/vsokolov/rh_book_repro/part6/logs/rest6_%j.log
# Part VI reproducibility: the tasks that depend on several chains (staging, comparisons, figures),
# any task the array left undone, then the comparison with the archived outputs.
echo "job $SLURM_JOB_ID on $(hostname) start $(date)"
cd /scratch/vsokolov/rh_book_repro/part6
/scratch/vsokolov/rh_book_repro/venv/bin/python3 -u run6.py rest
/scratch/vsokolov/rh_book_repro/venv/bin/python3 -u compare6.py > /dev/null
echo "end $(date)"
