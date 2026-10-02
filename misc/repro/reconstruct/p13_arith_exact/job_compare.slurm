#!/bin/bash
# Item arith_exact: redraw fig/heattrace.pdf from the new and from the archived arith.json, then compare (compare.py).
#SBATCH -J p13-arith-cmp
#SBATCH -p normal
#SBATCH --constraint=intel
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 30
#SBATCH -o /scratch/vsokolov/rh_book_repro/agents/p13_arith_exact/slurm-%j.log
cd /scratch/vsokolov/rh_book_repro/agents/p13_arith_exact
echo "[p13_arith_exact] job $SLURM_JOB_ID on $(hostname), $(lscpu | grep 'Model name' | sed 's/  */ /g')" >&2
/scratch/vsokolov/rh_book_repro/venv/bin/python3 -u heattrace.py arith.json heattrace.pdf
/scratch/vsokolov/rh_book_repro/venv/bin/python3 -u heattrace.py arith.json heattrace.png
/scratch/vsokolov/rh_book_repro/venv/bin/python3 -u heattrace.py arith_archived.json heattrace_archived_json.pdf
/scratch/vsokolov/rh_book_repro/venv/bin/python3 -u heattrace.py arith_archived.json heattrace_archived_json.png
/scratch/vsokolov/rh_book_repro/venv/bin/python3 -u compare.py
