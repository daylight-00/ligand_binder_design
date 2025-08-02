#!/bin/bash
#SBATCH -J boltz
#SBATCH -c 16
#SBATCH -q normal
#SBATCH -p gpu-farm
#SBATCH --gres=gpu:4
source /home/hwjang/miniforge3/bin/activate boltz2
source /home/hwjang/.slurmrc
slurm_start
###############################################
export BOLTZ_CACHE=/home/hwjang/aipd/boltz/cache
boltz predict \
    yaml/0 \
    --out_dir output \
    --num_workers 8 \
    --devices 4 \
###############################################
slurm_end