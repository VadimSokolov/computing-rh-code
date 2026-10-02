#!/bin/bash
#SBATCH -J rhb-reconstruct_primeside_probe
#SBATCH -p normal
#SBATCH --constraint=intel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 60
#SBATCH -o /scratch/vsokolov/rh_book_repro/agents/reconstruct_primeside/logs/probe_%j.log
# Tail model probe for results_main.json (main_tail_probe.py) and the check of the printed tables (book_tables_check.py).
A=/scratch/vsokolov/rh_book_repro/agents/reconstruct_primeside
PY=/scratch/vsokolov/rh_book_repro/venv/bin/python3
echo "job $SLURM_JOB_ID on $(hostname) start $(date)"; lscpu | grep -E '^Model name'
cd $A/main && cp $A/main_tail_probe.py . && $PY -u main_tail_probe.py $A/book/results_main.json; echo "probe rc=$?"
cd $A && echo "== printed tables against the regenerated files" && $PY -u book_tables_check.py $A/run/results_prime.json $A/run/extra.json
echo "== printed tables against the archived files" && $PY -u book_tables_check.py $A/book/results_prime.json $A/book/extra.json
echo "end $(date)"
