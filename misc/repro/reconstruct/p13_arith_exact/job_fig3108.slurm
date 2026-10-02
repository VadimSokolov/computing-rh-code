#!/bin/bash
# Item arith_exact: redraw fig/heattrace.pdf with matplotlib 3.10.8, the version that wrote the archived file (its
# PDF metadata), installed privately in mpl3108/, with the archived creation date, and compare the PDFs (pdfcmp.py).
#SBATCH -J p13-arith-fig
#SBATCH -p normal
#SBATCH --constraint=intel
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 20
#SBATCH -o /scratch/vsokolov/rh_book_repro/agents/p13_arith_exact/slurm-%j.log
cd /scratch/vsokolov/rh_book_repro/agents/p13_arith_exact
echo "[p13_arith_exact] job $SLURM_JOB_ID on $(hostname), $(lscpu | grep 'Model name' | sed 's/  */ /g')" >&2
export PYTHONPATH=/scratch/vsokolov/rh_book_repro/agents/p13_arith_exact/mpl3108 SOURCE_DATE_EPOCH=1790359548
/scratch/vsokolov/rh_book_repro/venv/bin/python3 -c "import matplotlib; print('matplotlib', matplotlib.__version__, matplotlib.__file__)"
/scratch/vsokolov/rh_book_repro/venv/bin/python3 -u heattrace.py arith.json heattrace_mpl3108.pdf
/scratch/vsokolov/rh_book_repro/venv/bin/python3 -u heattrace.py arith_archived.json heattrace_archived_json_mpl3108.pdf
echo "== archived json, 3.10.8, against fig/heattrace.pdf"; /scratch/vsokolov/rh_book_repro/venv/bin/python3 -u pdfcmp.py heattrace_archived_json_mpl3108.pdf old_heattrace.pdf
echo "== new json, 3.10.8, against fig/heattrace.pdf"; /scratch/vsokolov/rh_book_repro/venv/bin/python3 -u pdfcmp.py heattrace_mpl3108.pdf old_heattrace.pdf
echo "== new json against archived json, both 3.10.8"; /scratch/vsokolov/rh_book_repro/venv/bin/python3 -u pdfcmp.py heattrace_mpl3108.pdf heattrace_archived_json_mpl3108.pdf
echo "== new json, 3.11.2 (job 1243544), against fig/heattrace.pdf"; /scratch/vsokolov/rh_book_repro/venv/bin/python3 -u pdfcmp.py heattrace.pdf old_heattrace.pdf | head -4
