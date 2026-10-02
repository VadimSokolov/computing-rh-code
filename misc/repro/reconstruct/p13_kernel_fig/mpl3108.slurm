#!/bin/bash
# mpl3108.slurm (round 2 item kernel_fig): runs kernel_eta.py in both modes with matplotlib 3.10.8, the version that
# drew the archived fig/kernel.pdf, instead of the 3.11.2 of the shared venv. The wheel was unpacked with
#   pip install --no-deps --target=mpl3108_pkg matplotlib-3.10.8-cp311-cp311-manylinux2014_x86_64.whl
# into the agent folder; every other package (python-flint 0.9, numpy, pillow, ...) comes from the shared venv.
#SBATCH -J rhb-p13_kernel_fig-mpl3108
#SBATCH -p normal
#SBATCH --constraint=intel
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -t 20
#SBATCH -o /scratch/vsokolov/rh_book_repro/agents/p13_kernel_fig/mpl3108/slurm-%j.log
cd /scratch/vsokolov/rh_book_repro/agents/p13_kernel_fig/mpl3108
echo "job $SLURM_JOB_ID on $(hostname), $(grep -m1 'model name' /proc/cpuinfo)"
export PYTHONPATH=/scratch/vsokolov/rh_book_repro/agents/p13_kernel_fig/mpl3108_pkg
/scratch/vsokolov/rh_book_repro/venv/bin/python3 -u kernel_eta.py old
/scratch/vsokolov/rh_book_repro/venv/bin/python3 -u kernel_eta.py
