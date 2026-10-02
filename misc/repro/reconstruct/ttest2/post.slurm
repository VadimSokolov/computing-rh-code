#!/bin/bash
# Reconstruction job, not the authors' script: collects the parameter search (collect_search.py) and compares
# Table tab:tilt of ch/ch04.tex with the archived code/ttest2.py output (compare_table.py). Submitted with
# --dependency=afterany on the search array, the archived run and the diagnostics; ARCHLOG names the log of
# run_archived.slurm.
#SBATCH -J rhb-ttest2-post
#SBATCH -p normal
#SBATCH --constraint=intel
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -t 10
#SBATCH -o /scratch/vsokolov/rh_book_repro/agents/reconstruct_ttest2/out/post_%j.log
cd /scratch/vsokolov/rh_book_repro/agents/reconstruct_ttest2
PY=/scratch/vsokolov/rh_book_repro/venv/bin/python3
echo "== collect_search.py"; $PY -u collect_search.py
echo "== compare_table.py $ARCHLOG"; $PY -u compare_table.py "$ARCHLOG"
