#!/bin/bash
#SBATCH -J rhb-bosch_v7_tal
#SBATCH -p normal
#SBATCH --constraint=intel
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 720
#SBATCH --array=0-1
#SBATCH -o /scratch/vsokolov/rh_book_repro/agents/reconstruct_bosch_v7/logs/talbot_%A_%a.log
# Talbot inversion of phi at 150 digits (talbot_hp.py), task 0 at t=0.003 and task 1 at t=0.005. Submit: sbatch talbot.slurm
PY=/scratch/vsokolov/rh_book_repro/venv/bin/python3
D=/scratch/vsokolov/rh_book_repro/agents/reconstruct_bosch_v7
T=(0.003 0.005)
echo "task $SLURM_ARRAY_TASK_ID of job $SLURM_ARRAY_JOB_ID on $(hostname): $(grep -m1 'model name' /proc/cpuinfo) start $(date)"
mkdir -p $D/talbot; cd $D/talbot
$PY -u $D/src/talbot_hp.py ${T[$SLURM_ARRAY_TASK_ID]}; echo "rc=$? end $(date)"
