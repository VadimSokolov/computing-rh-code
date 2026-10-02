#!/bin/bash
#SBATCH -J rhb-bosch_v7_res
#SBATCH -p normal
#SBATCH --constraint=intel
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 480
#SBATCH --array=0-49
#SBATCH -o /scratch/vsokolov/rh_book_repro/agents/reconstruct_bosch_v7/logs/residue_%A_%a.log
# Block 6 of code_bosch/verify_bosch_xi_v7.py at 150 digits: zeros n = TASK+1, TASK+51, ... up to 250 and their
# residue coefficients (residue_hp_task.py); residue_hp_collect.py sums the series. Submit: sbatch residue.slurm
PY=/scratch/vsokolov/rh_book_repro/venv/bin/python3
D=/scratch/vsokolov/rh_book_repro/agents/reconstruct_bosch_v7
echo "task $SLURM_ARRAY_TASK_ID of job $SLURM_ARRAY_JOB_ID on $(hostname): $(grep -m1 'model name' /proc/cpuinfo) start $(date)"
mkdir -p $D/residue; cd $D/residue
$PY -u $D/src/residue_hp_task.py $SLURM_ARRAY_TASK_ID 50 250; echo "rc=$? end $(date)"
