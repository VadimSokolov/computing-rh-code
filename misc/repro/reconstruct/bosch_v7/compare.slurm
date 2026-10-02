#!/bin/bash
#SBATCH -J rhb-bosch_v7_cmp2
#SBATCH -p normal
#SBATCH --constraint=intel
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 30
#SBATCH -o /scratch/vsokolov/rh_book_repro/agents/reconstruct_bosch_v7/logs/compare_%j.log
# Reruns only the comparison (compare_bosch_v7.py) after a fix to one of its regular expressions; the collector
# output residue/residue_hp.json of job 1243307 is reused. Submit: sbatch --wait compare.slurm
PY=/scratch/vsokolov/rh_book_repro/venv/bin/python3
D=/scratch/vsokolov/rh_book_repro/agents/reconstruct_bosch_v7
echo "job $SLURM_JOB_ID on $(hostname): $(grep -m1 'model name' /proc/cpuinfo) start $(date)"
cd $D && md5sum src/compare_bosch_v7.py && $PY -u src/compare_bosch_v7.py; echo "compare rc=$? end $(date)"
