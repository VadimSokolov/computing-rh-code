#!/bin/bash
#SBATCH -J rhb-bosch_v7_cmp
#SBATCH -p normal
#SBATCH --constraint=intel
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 60
#SBATCH -o /scratch/vsokolov/rh_book_repro/agents/reconstruct_bosch_v7/logs/collect_%j.log
# Collects the 150 digit residue tasks (residue_hp_collect.py) and runs the comparison (compare_bosch_v7.py).
# Submit after main.slurm, residue.slurm and talbot.slurm: sbatch --dependency=afterany:<main job> collect.slurm
PY=/scratch/vsokolov/rh_book_repro/venv/bin/python3
D=/scratch/vsokolov/rh_book_repro/agents/reconstruct_bosch_v7
echo "job $SLURM_JOB_ID on $(hostname): $(grep -m1 'model name' /proc/cpuinfo) start $(date)"
cd $D/residue && ls residue_*.json | wc -l && $PY -u $D/src/residue_hp_collect.py; echo "collect rc=$?"
cd $D && md5sum src/verify_bosch_xi_v7_output.txt fixed/joined_output.txt equiv/joined_output.txt && $PY -u src/compare_bosch_v7.py; echo "compare rc=$? end $(date)"
