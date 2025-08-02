#!/bin/bash
#SBATCH -J dssp
#SBATCH -c 64
#SBATCH -q high_cpu_users
#SBATCH -p cpu-farm
source /home/hwjang/miniforge3/bin/activate pyrosetta
source /home/hwjang/.slurmrc
slurm_start
###############################################
python 1_pyrosetta_dssp.py \
    pdb_list.csv \
    ../0_params/PHT.params
###############################################
slurm_end