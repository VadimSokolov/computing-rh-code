#!/bin/bash
#SBATCH --job-name=rhb_v1chk
#SBATCH --partition=normal
#SBATCH --constraint=intel
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0-02:00:00
#SBATCH --output=/scratch/vsokolov/rh_book_repro/v1check/logs/%x_%A_%a.log
# Array task: MODE=grid (100 tasks, v1_grid_mp.py) or MODE=zeros (37 tasks, zeros_mp.py). See misc/repro/v1_collect.py.
cd /scratch/vsokolov/rh_book_repro/v1check
/scratch/vsokolov/rh_book_repro/venv/bin/python3 -u ${MODE}.py $SLURM_ARRAY_TASK_ID $NTASK
