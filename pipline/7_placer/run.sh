#!/bin/bash
#SBATCH -J placer
#SBATCH -c 4
##SBATCH -q normal
#SBATCH -p gpu-farm
#SBATCH --gres=gpu:1
#SBATCH -w node25
#SBATCH --array=0-3
source /home/hwjang/miniforge3/bin/activate rfdaa
source /home/hwjang/.slurmrc
slurm_start
###############################################
export PYTHONPATH="${PYTHONPATH}:/home/hwjang/aipd/PLACER/"
python placer.py placer_input/split_$SLURM_ARRAY_TASK_ID.parquet
###############################################
slurm_end