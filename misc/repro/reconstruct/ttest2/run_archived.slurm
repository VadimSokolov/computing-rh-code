#!/bin/bash
# Reconstruction job, not the authors' script. It runs the archived code/ttest2.py and code/ttest.py
# unchanged (copied beside this file with code/tilted.py) with the command lines that reproduce
# Table tab:tilt of Chapter ch:polya (ch/ch04.tex). The book records no command line; these are the
# ones the parameter search (search.slurm) singles out.
#SBATCH -J rhb-ttest2-archived
#SBATCH -p normal
#SBATCH --constraint=intel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 60
#SBATCH -o /scratch/vsokolov/rh_book_repro/agents/reconstruct_ttest2/out/archived_%j.log
cd /scratch/vsokolov/rh_book_repro/agents/reconstruct_ttest2
PY=/scratch/vsokolov/rh_book_repro/venv/bin/python3
echo "job $SLURM_JOB_ID on $(hostname), $(lscpu | grep 'Model name' | sed 's/  */ /g')"
$PY -c "import flint, mpmath; print('python-flint', flint.__version__, 'mpmath', mpmath.__version__)"
md5sum ttest.py ttest2.py tilted.py
echo "== ttest2.py 0.10 200 0.004 4 200 500 1000"
$PY -u ttest2.py 0.10 200 0.004 4 200 500 1000
echo "== ttest2.py 0.02 220 0.001 4.2 1000 2000 5000"
$PY -u ttest2.py 0.02 220 0.001 4.2 1000 2000 5000
echo "== ttest.py 0.1 200 0.004 4 (Horner evaluator, heights fixed in the script)"
$PY -u ttest.py 0.1 200 0.004 4
echo "== done"
