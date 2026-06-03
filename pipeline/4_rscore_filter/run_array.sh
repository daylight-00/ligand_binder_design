#!/bin/bash
#SBATCH -J rscore
#SBATCH -c 64
#SBATCH -q high_cpu_users
#SBATCH -p cpu-farm
#SBATCH --array=0-3
source /home/hwjang/miniforge3/bin/activate pyrosetta
source /home/hwjang/.slurmrc
slurm_start
###############################################
export PYTHONPATH="${PYTHONPATH}:./relax_score_simple"
python 1_relax_sc_simple_hb_mp.py \
    "rscore_input/split_$SLURM_ARRAY_TASK_ID.parquet"
###############################################
slurm_end