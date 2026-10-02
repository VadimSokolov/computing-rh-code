#!/bin/bash
#SBATCH -J rhb-reconstruct_primeside
#SBATCH -p normal
#SBATCH --constraint=intel
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 90
#SBATCH -o /scratch/vsokolov/rh_book_repro/agents/reconstruct_primeside/logs/pipeline_%j.log
# Reconstruction of the gap "primeside" (misc/repro/reconstruct/primeside/REPORT.md).
# Runs the archived code/primeside.py unchanged, with the shim polya.py beside it, then the
# archived code/plots.py on its output, then compares with the archived data files.
A=/scratch/vsokolov/rh_book_repro/agents/reconstruct_primeside
R=/scratch/vsokolov/rh_book_repro/run            # earlier full rerun: read only
PY=/scratch/vsokolov/rh_book_repro/venv/bin/python3
echo "job $SLURM_JOB_ID on $(hostname) start $(date)"
lscpu | grep -E '^Model name' ; $PY -c "import sys,numpy,scipy,mpmath,flint,matplotlib; print(sys.version.split()[0], 'numpy', numpy.__version__, 'scipy', scipy.__version__, 'mpmath', mpmath.__version__, 'flint', flint.__version__, 'matplotlib', matplotlib.__version__)"
echo "md5 of the archived scripts used:"; md5sum $A/code/*.py
# 1. forensic check of the shim convention
rm -rf $A/variants && mkdir -p $A/variants && cd $A/variants
cp $A/code/polya_flint.py $A/shim_variants.py .
cp $A/book/results_prime.json book_results_prime.json; cp $A/book/extra.json book_extra.json
t0=$(date +%s); $PY -u shim_variants.py; echo "shim_variants rc=$? $(( $(date +%s) - t0 ))s"
# 2. primeside.py with the shim, then plots.py
rm -rf $A/run && mkdir -p $A/run/fig && cd $A/run
cp $A/code/primeside.py $A/code/plots.py $A/code/polya_flint.py $A/polya.py .
cp $R/plotcurves.npz $R/analysis.json $R/hcm.json .
t0=$(date +%s); $PY -u primeside.py > primeside.out 2>&1; rc=$?; echo "primeside rc=$rc $(( $(date +%s) - t0 ))s"; cat primeside.out
t0=$(date +%s); $PY -u plots.py > plots.out 2>&1; rc=$?; echo "plots rc=$rc $(( $(date +%s) - t0 ))s"; cat plots.out
ls -la $A/run $A/run/fig
# 3. comparisons
cd $A
$PY -u compare_primeside.py $A/book $A/run
$PY -u archimedean_check.py $A/book/results_prime.json $A/run/results_prime.json
echo "end $(date)"
