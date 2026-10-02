#!/bin/bash
#SBATCH --job-name=rhbook_rzcmp
#SBATCH --partition=normal
#SBATCH --constraint=intel
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=0-01:00:00
#SBATCH --output=/scratch/vsokolov/rh_book_repro/riesz/logs/compare_%j.log
# Collects the task status files and compares every archived output of the Riesz package with the run.
cd /scratch/vsokolov/rh_book_repro/riesz
/scratch/vsokolov/rh_book_repro/venv/bin/python3 -c "import json,glob,os; json.dump({os.path.basename(p)[:-5]: json.load(open(p)) for p in glob.glob('status/*.json')}, open('status.json','w'), indent=1)"
/scratch/vsokolov/rh_book_repro/venv/bin/python3 -u compare6.py /scratch/vsokolov/rh_book_repro/riesz > /dev/null
echo "compared $(date)"
