#!/bin/bash
# Submit an sbatch array on Hopper (Intel nodes): one task per entry of a tasks JSON file.
# Usage: bash submit_array.sh NAME SCRIPT TASKS.json MINUTES MEM [MAXPAR] -- files...
# Runs "python3 SCRIPT $SLURM_ARRAY_TASK_ID TASKS.json" per task, waits for the whole array,
# then copies result_*.json and the logs back into ./out_NAME/.
set -e
name=$1; script=$2; tasks=$3; mins=$4; mem=$5; shift 5
maxpar=400
if [ "$1" != "--" ]; then maxpar=$1; shift; fi
[ "$1" = "--" ] && shift
dir=/scratch/vsokolov/rh_book_repro/agents/N3-$name
N=$(jq length "$tasks")
ssh -o BatchMode=yes hopper "mkdir -p $dir"
scp -q -o BatchMode=yes "$script" "$tasks" "$@" "hopper:$dir/"
tb=$(basename "$tasks"); sb=$(basename "$script")
job="#!/bin/bash
#SBATCH -J N3-$name
#SBATCH -p normal
#SBATCH --constraint=intel
#SBATCH -c 1
#SBATCH --mem=$mem
#SBATCH -t $mins
#SBATCH --array=0-$((N-1))%$maxpar
#SBATCH -o $dir/slurm-%A_%a.log
cd $dir
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
exec /scratch/vsokolov/rh_book_repro/venv/bin/python3 -u $sb \$SLURM_ARRAY_TASK_ID $tb"
printf '%s\n' "$job" | ssh -o BatchMode=yes hopper "cat > $dir/job.slurm && sbatch --wait -Q $dir/job.slurm >/dev/null; echo done"
mkdir -p out_$name
scp -q -o BatchMode=yes "hopper:$dir/result_*.json" out_$name/ 2>/dev/null || true
ssh -o BatchMode=yes hopper "cd $dir && tar czf logs.tgz slurm-*.log" && scp -q -o BatchMode=yes "hopper:$dir/logs.tgz" out_$name/ || true
ls out_$name | wc -l
