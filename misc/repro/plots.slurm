#!/bin/bash
#SBATCH --job-name=rhbook_plots
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --exclude=amd[021-060]
#SBATCH --output=/scratch/vsokolov/rh_book_repro/logs/plots_%j.log
# Question: do sens2.json and extra.json (written by code/plots.py) reproduce? primeside.py cannot run
# (it imports a module "polya" missing from the archive), so plots.py reads the archived results_prime.json.
cd /scratch/vsokolov/rh_book_repro/run
cp ../book_data/results_prime.json .
../venv/bin/python3 -u plots.py && ../venv/bin/python3 -u ../compare.py > /dev/null && grep -E "sens2|extra" ../compare_report.md
