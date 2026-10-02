#!/bin/bash
#SBATCH -J rhb-reconstruct_spacing
#SBATCH -p normal
#SBATCH --constraint=intel
#SBATCH -c 16
#SBATCH --mem=16G
#SBATCH -t 60
#SBATCH -o /scratch/vsokolov/rh_book_repro/agents/reconstruct_spacing/slurm-%j.log
# Reconstruction of data/spacing.json (misc/repro/reconstruct/spacing/). spacing.py runs twice: in fresh/ with
# ordinates that the script computes itself (mpmath.zetazero at 30 digits), and in npy/ with the zero lists of the
# earlier full rerun (/scratch/vsokolov/rh_book_repro/run, made by misc/repro/zeros.py). Both outputs are compared
# with the archived file (book_data/spacing.json, a copy of data/spacing.json), and the candidate figure with the
# book's fig/spacing.pdf (copied beside this file as book_fig_spacing.pdf). Submit from the folder below:
# sbatch --wait run_spacing.slurm. Job 1241556 ran an earlier version of this file and of spacing.py, with the step
# "$PY -u variants.py $BOOK fresh" (out/variants.json); variants2.py, variants3.py and pdf_bars.py ran through
# misc/tools/hopper_run.sh. REPORT.md lists the jobs.
set -e
PY=/scratch/vsokolov/rh_book_repro/venv/bin/python3
D=/scratch/vsokolov/rh_book_repro/agents/reconstruct_spacing
BOOK=/scratch/vsokolov/rh_book_repro/book_data/spacing.json
RUN=/scratch/vsokolov/rh_book_repro/run
cd $D
echo "job $SLURM_JOB_ID on $(hostname): $(grep -m1 'model name' /proc/cpuinfo)"
md5sum $BOOK
rm -rf fresh npy; mkdir fresh npy
cp spacing.py fresh/; cp spacing.py npy/
cp $RUN/g_1_400.npy $RUN/g_401_1000.npy $RUN/g_1001_1700.npy npy/
echo "== spacing.py, ordinates computed by the script"; (cd fresh && $PY -u spacing.py)
echo "== spacing.py, zero lists of the earlier rerun"; (cd npy && $PY -u spacing.py)
echo "== zero lists: fresh against the earlier rerun"; for f in g_1_400.npy g_401_1000.npy g_1001_1700.npy; do cmp fresh/$f npy/$f && echo "$f identical"; done
echo "== extras and histogram: fresh against npy"; diff fresh/spacing_extra.json npy/spacing_extra.json || true; cmp fresh/spacing_hist.json npy/spacing_hist.json && echo "spacing_hist.json identical"
echo "== comparison with the archived file"; $PY -u compare_spacing.py $BOOK fresh/spacing.json npy/spacing.json
cmp $BOOK fresh/spacing.json && echo "fresh/spacing.json is byte identical to the archived file"
echo "== candidate figure against the book's fig/spacing.pdf"; md5sum book_fig_spacing.pdf
$PY -u fig_check.py book_fig_spacing.pdf fresh/spacing.pdf fresh/spacing_hist.json
echo "== end"
