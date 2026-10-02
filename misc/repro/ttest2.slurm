#!/bin/bash
#SBATCH --job-name=rhbook_ttest2
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --exclude=amd[021-060]
#SBATCH --output=/scratch/vsokolov/rh_book_repro/logs/ttest2_%j.log
# Question: does code/ttest2.py reproduce the orders of magnitude in Table tab:tilt (ch04)?
# The book does not record h and xmax; these values give the node counts 2001 and 8401 of the table.
cd /scratch/vsokolov/rh_book_repro/run
echo "delta=0.10 prec=200 h=0.004 xmax=4 (2001 nodes)"
../venv/bin/python3 -u ttest2.py 0.10 200 0.004 4 200 500 1000
echo "delta=0.02 prec=220 h=0.001 xmax=4.2 (8401 nodes)"
../venv/bin/python3 -u ttest2.py 0.02 220 0.001 4.2 1000 2000 5000
