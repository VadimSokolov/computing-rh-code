#!/bin/bash
#SBATCH -J rhb-R2-23-final
#SBATCH -p normal
#SBATCH --constraint=intel
#SBATCH -c 16
#SBATCH --mem=16G
#SBATCH -t 60
#SBATCH --array=1,3,5
#SBATCH -o /scratch/vsokolov/rh_book_repro/agents/R2-23/logs/final_%A_%a.log
# Unit R2-23: checks of the scripts added to or changed in the archive (misc/repro/archive-changes.md).
# One task per script, each in a fresh folder, with the files exactly as they are now in the book folder (src/).
PY=/scratch/vsokolov/rh_book_repro/venv/bin/python3
D=/scratch/vsokolov/rh_book_repro/agents/R2-23
S=$D/src; H=$S/helpers; RUN=/scratch/vsokolov/rh_book_repro/run
echo "task $SLURM_ARRAY_TASK_ID of job $SLURM_ARRAY_JOB_ID on $(hostname): $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2) start $(date)"
$PY -c "import sys,numpy,scipy,mpmath,flint,matplotlib; print('python', sys.version.split()[0], 'numpy', numpy.__version__, 'scipy', scipy.__version__, 'mpmath', mpmath.__version__, 'flint', flint.__version__, 'matplotlib', matplotlib.__version__)"
t0=$(date +%s)
case $SLURM_ARRAY_TASK_ID in
0)  # (a) code/hankel.py beside the zero lists of riesz_program/zeros
    W=$D/final_t0_hankel; rm -rf $W; mkdir -p $W; cd $W
    cp $S/code/hankel.py $S/zeros/g_*.npy .; md5sum hankel.py g_*.npy
    $PY -u hankel.py; echo "rc=$?"
    $PY -u $H/cmpjson.py $S/data/hankel.json hankel.json 0
    cmp $S/data/hankel.json hankel.json && echo "hankel.json is byte identical to data/hankel.json"; md5sum hankel.json ;;
1)  # (b) code/spacing.py, once with the zero lists and once computing them, then the figure check
    W=$D/final_t1_spacing; rm -rf $W; mkdir -p $W/npy $W/fresh; cd $W
    cp $S/code/spacing.py npy/; cp $S/zeros/g_*.npy npy/; cp $S/code/spacing.py fresh/; md5sum npy/spacing.py
    echo "== with the zero lists"; (cd npy && $PY -u spacing.py; echo "rc=$?")
    echo "== without the zero lists"; (cd fresh && $PY -u spacing.py; echo "rc=$?")
    for f in g_1_400.npy g_401_1000.npy g_1001_1700.npy; do cmp fresh/$f $S/zeros/$f && echo "computed $f is byte identical to riesz_program/zeros/$f"; done
    for v in npy fresh; do echo "== $v/spacing.json against data/spacing.json"; $PY -u $H/cmpjson.py $S/data/spacing.json $v/spacing.json 0
        cmp $S/data/spacing.json $v/spacing.json && echo "$v/spacing.json is byte identical to data/spacing.json"; done
    md5sum npy/spacing.json fresh/spacing.json $S/data/spacing.json
    cmp npy/spacing.pdf fresh/spacing.pdf && echo "the two spacing.pdf are byte identical"
    echo "== bars of the book's fig/spacing.pdf and of the new spacing.pdf against the new histogram"
    cp $H/pdf_bars.py $H/fig_check.py .; $PY -u fig_check.py $S/fig/spacing.pdf npy/spacing.pdf npy/spacing_hist.json ;;
2)  # (c) code/sinh_bench.py, default (K = 400000) and closed
    W=$D/final_t2_sinh_bench; rm -rf $W; mkdir -p $W; cd $W
    cp $S/code/sinh_bench.py .; md5sum sinh_bench.py
    $PY -u sinh_bench.py; echo "rc=$?"; $PY -u sinh_bench.py closed; echo "rc=$?"
    echo "== default run against data/sinh_bench.json (tolerance 1.1e-12)"; $PY -u $H/cmpjson.py $S/data/sinh_bench.json sinh_bench.json 1.1e-12
    echo "== closed against the reconstruction's sinh_bench_closed.json"; $PY -u $H/cmpjson.py $H/sinh_bench_closed.json sinh_bench_closed.json 0 ;;
3)  # (d) basepoint_one/wiener/riesz_explicit.py: from the archived salem_riesz.json, and after salem_riesz.py
    W=$D/final_t3_riesz_explicit; rm -rf $W; mkdir -p $W/direct $W/chain; cd $W
    cp $S/wiener/riesz_explicit.py $S/wiener/salem_riesz.json $S/wiener/g_1_400.npy direct/
    cp $S/wiener/riesz_explicit.py $S/wiener/salem_riesz.py $S/wiener/g_1_400.npy chain/; md5sum direct/riesz_explicit.py
    echo "== riesz_explicit.py with the archived salem_riesz.json"; (cd direct && $PY -u riesz_explicit.py; echo "rc=$?")
    echo "== salem_riesz.py; riesz_explicit.py"; (cd chain && $PY -u salem_riesz.py > salem.out; echo "salem rc=$?"; $PY -u riesz_explicit.py > explicit.out; echo "explicit rc=$?"; tail -3 explicit.out)
    cmp chain/salem_riesz.json $S/wiener/salem_riesz.json && echo "recomputed salem_riesz.json is byte identical to the archived file"
    for v in direct chain; do echo "== $v/riesz_explicit.json against basepoint_one/wiener/riesz_explicit.json"; $PY -u $H/cmpjson.py $S/wiener/riesz_explicit.json $v/riesz_explicit.json 0
        cmp $S/wiener/riesz_explicit.json $v/riesz_explicit.json && echo "$v/riesz_explicit.json is byte identical to the archived file"; done ;;
4)  # (e) code/ttest2.py and code/ttest.py with the command lines now recorded at their tops
    W=$D/final_t4_ttest; rm -rf $W; mkdir -p $W; cd $W
    cp $S/code/ttest.py $S/code/ttest2.py $S/code/tilted.py .; md5sum ttest.py ttest2.py tilted.py; head -5 ttest2.py; head -5 ttest.py
    { echo "== ttest2.py 0.10 200 0.004 4 200 500 1000"; $PY -u ttest2.py 0.10 200 0.004 4 200 500 1000
      echo "== ttest2.py 0.02 220 0.001 4.2 1000 2000 5000"; $PY -u ttest2.py 0.02 220 0.001 4.2 1000 2000 5000
      echo "== ttest.py 0.1 200 0.004 4 (Horner evaluator, heights fixed in the script)"; $PY -u ttest.py 0.1 200 0.004 4
      echo "== done"; } > new.log 2>&1
    cat new.log; $PY -u $H/cmp_ttest.py $H/archived_1241681.log new.log ;;
5)  # (f) code/primeside.py with code/polya.py, then the HEAD version of code/plots.py on its output
    W=$D/final_t5_primeside; rm -rf $W; mkdir -p $W/fig; cd $W
    cp $S/code/primeside.py $S/code/polya.py $S/code/polya_flint.py .; cp $S/code/plots_HEAD.py plots.py
    cp $RUN/plotcurves.npz $RUN/analysis.json $RUN/hcm.json .; md5sum *.py
    $PY -u primeside.py > primeside.out 2>&1; echo "primeside rc=$?"; cat primeside.out
    echo "== results_prime.json against data/results_prime.json"; $PY -u $H/cmpjson.py $S/data/results_prime.json results_prime.json 1e-9
    $PY -u plots.py > plots.out 2>&1; echo "plots rc=$?"; tail -3 plots.out
    for f in extra.json sens2.json; do cmp $S/data/$f $f && echo "$f is byte identical to data/$f"; done ;;
6)  # (g) verification_reciprocal_zeta/spot_checks.py with the command now in its docstring
    W=$D/final_t6_spot_checks; rm -rf $W; mkdir -p $W; cd $W
    cp $S/verif/spot_checks.py .; sha256sum spot_checks.py
    $PY spot_checks.py > spot_checks_report.md 2> stderr.txt; echo "rc=$? stderr bytes $(wc -c < stderr.txt)"
    sha256sum spot_checks_report.md $S/verif/spot_checks_report.md; wc -l -c spot_checks_report.md
    cmp spot_checks_report.md $S/verif/spot_checks_report.md && echo "the report is byte identical to verification_reciprocal_zeta/spot_checks_report.md" ;;
7)  # (h) code_bosch/verify_bosch_xi_v7.py with the gaptest.py lines fixed, in an empty folder; outputs joined in block order
    W=$D/final_t7_bosch; rm -rf $W; mkdir -p $W; cd $W
    cmp $S/bosch/verify_bosch_xi_v7.py $S/bosch/verify_bosch_xi_v7_hcm.py && echo "code_bosch and bosch_hcm_xi copies are byte identical"
    cp $S/bosch/verify_bosch_xi_v7.py .; md5sum verify_bosch_xi_v7.py; wc -l verify_bosch_xi_v7.py
    $PY -u verify_bosch_xi_v7.py > stdout.txt 2> stderr.txt; echo "rc=$?"; cat stderr.txt
    JOIN="verify_v2_output.txt verify_v3_output.txt thorin_p_output.txt open_problems_output.txt op5_output.txt gW_output.txt gpos_output.txt gaptest_output.txt gaptest3_output.txt"
    cat $JOIN > joined_output.txt; wc -l $JOIN joined_output.txt
    $PY -u $H/cmp_bosch.py joined_output.txt $S/bosch/reconstruct_fixed_joined.txt $S/bosch/archived_output.txt ;;
esac
echo "seconds=$(( $(date +%s) - t0 )) end $(date)"
