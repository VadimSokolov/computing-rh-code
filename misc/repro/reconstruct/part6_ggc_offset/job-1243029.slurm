#!/bin/bash
#SBATCH -J rhb-part6_ggc_offset/run2
#SBATCH -p normal
#SBATCH --constraint=intel
#SBATCH -c 16
#SBATCH --mem=8G
#SBATCH -t 60
#SBATCH -o /scratch/vsokolov/rh_book_repro/agents/part6_ggc_offset/run2/slurm-%j.log
cd /scratch/vsokolov/rh_book_repro/agents/part6_ggc_offset/run2
echo "[hopper_run] job $SLURM_JOB_ID on $(hostname)" >&2
exec /scratch/vsokolov/rh_book_repro/venv/bin/python3 -u ggc_offset.py
