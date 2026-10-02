#!/bin/bash
#SBATCH -J rhb-bosch_v7
#SBATCH -p normal
#SBATCH --constraint=intel
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 120
#SBATCH --array=0-4
#SBATCH -o /scratch/vsokolov/rh_book_repro/agents/reconstruct_bosch_v7/logs/main_%A_%a.log
# Reconstruction of the bosch_v7 gap (misc/repro/reconstruct/bosch_v7). One independent run per array task, each in
# an empty folder: 0 the fixed copy verify_bosch_xi_v7_fixed.py; 1 the archived script unchanged, with the
# reconstructed gaptest.py beside it; 2 supp_lines37_39.py; 3 supp_line109.py; 4 supp_gaps_table.py.
# Tasks 0 and 1 join the block output files in block order, as the printed output does (Appendix app:bosch).
# Submit: sbatch main.slurm (scripts in src/).
PY=/scratch/vsokolov/rh_book_repro/venv/bin/python3
D=/scratch/vsokolov/rh_book_repro/agents/reconstruct_bosch_v7
echo "task $SLURM_ARRAY_TASK_ID of job $SLURM_ARRAY_JOB_ID on $(hostname): $(grep -m1 'model name' /proc/cpuinfo) start $(date)"
JOIN="verify_v2_output.txt verify_v3_output.txt thorin_p_output.txt open_problems_output.txt op5_output.txt gW_output.txt gpos_output.txt gaptest_output.txt gaptest3_output.txt"
case $SLURM_ARRAY_TASK_ID in
0) W=$D/fixed; S=verify_bosch_xi_v7_fixed.py; FILES="$D/src/$S";;
1) W=$D/equiv; S=verify_bosch_xi_v7.py; FILES="$D/src/$S $D/src/gaptest.py";;
2) W=$D/supp37; S=supp_lines37_39.py; FILES="$D/src/$S";;
3) W=$D/supp109; S=supp_line109.py; FILES="$D/src/$S";;
4) W=$D/suppgaps; S=supp_gaps_table.py; FILES="$D/src/$S";;
esac
rm -rf $W; mkdir -p $W; cp $FILES $W/; cd $W; md5sum *.py
t0=$(date +%s)
$PY -u $S > stdout.txt 2>&1; rc=$?
echo "rc=$rc seconds=$(( $(date +%s) - t0 ))"
cat stdout.txt
if [ $SLURM_ARRAY_TASK_ID -le 1 ]; then cat $JOIN > joined_output.txt; wc -l $JOIN joined_output.txt; fi
echo "end $(date)"
