#!/bin/bash
#SBATCH -J rhb-reconstruct_primeside_probe2
#SBATCH -p normal
#SBATCH --constraint=intel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 60
#SBATCH -o /scratch/vsokolov/rh_book_repro/agents/reconstruct_primeside/logs/probe2_%j.log
# Second tail probe for results_main.json, and the number by number comparison of the rebuilt file.
A=/scratch/vsokolov/rh_book_repro/agents/reconstruct_primeside
PY=/scratch/vsokolov/rh_book_repro/venv/bin/python3
echo "job $SLURM_JOB_ID on $(hostname) start $(date)"; lscpu | grep -E '^Model name'
cd $A/main && cp $A/main_tail_probe2.py . && $PY -u main_tail_probe2.py $A/book/results_main.json; echo "probe2 rc=$?"
rm -rf $A/cmp_main && mkdir -p $A/cmp_main && cp $A/main/main_reconstructed.json $A/cmp_main/results_main.json
cd $A/cmp_main && $PY -u $A/compare_primeside.py $A/book $A/cmp_main
echo "end $(date)"
