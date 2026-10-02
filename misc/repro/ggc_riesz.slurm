#!/bin/bash
#SBATCH --job-name=rhbook_ggc
#SBATCH --partition=normal
#SBATCH --constraint=intel
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=0-12:00:00
#SBATCH --output=/scratch/vsokolov/rh_book_repro/riesz/logs/ggc_%j.log
# Adds ggc_coef.py (revised Riesz package) to an existing Riesz run without redoing the other tasks,
# then compares every archived output again.
echo "job $SLURM_JOB_ID on $(hostname) start $(date)"
cd /scratch/vsokolov/rh_book_repro/riesz
cp src/code/ggc_coef.py run/code/ && cp src/data/ggc_coef.json ref/code/
/scratch/vsokolov/rh_book_repro/venv/bin/python3 -u run_riesz.py group 8
/scratch/vsokolov/rh_book_repro/venv/bin/python3 -c "import json,glob,os; json.dump({os.path.basename(p)[:-5]: json.load(open(p)) for p in glob.glob('status/*.json')}, open('status.json','w'), indent=1)"
/scratch/vsokolov/rh_book_repro/venv/bin/python3 -u compare6.py /scratch/vsokolov/rh_book_repro/riesz > /dev/null
echo "end $(date)"
