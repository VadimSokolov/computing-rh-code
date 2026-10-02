#!/bin/bash
#SBATCH -J rhb-reconstruct_primeside_v2
#SBATCH -p normal
#SBATCH --constraint=intel
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH -t 60
#SBATCH -o /scratch/vsokolov/rh_book_repro/agents/reconstruct_primeside/logs/variants2_%j.log
# Second forensic pass on the missing module polya (see shim_variants2.py).
A=/scratch/vsokolov/rh_book_repro/agents/reconstruct_primeside
PY=/scratch/vsokolov/rh_book_repro/venv/bin/python3
echo "job $SLURM_JOB_ID on $(hostname) start $(date)"; lscpu | grep -E '^Model name'
rm -rf $A/variants2 && mkdir -p $A/variants2 && cd $A/variants2
cp $A/code/polya_flint.py $A/shim_variants2.py . && cp $A/book/results_prime.json book_results_prime.json
$PY -u shim_variants2.py; echo "shim_variants2 rc=$? end $(date)"
