#!/bin/bash
#SBATCH -J sasa
#SBATCH -c 64
#SBATCH -q normal
#SBATCH -p cpu-farm
source /home/hwjang/miniforge3/bin/activate pyrosetta
source /home/hwjang/.slurmrc
slurm_start
###############################################
python 2_pyrosetta_lig_atm_sasa.py \
    pdb_list.csv \
    ../0_params/PHT.params
###############################################
slurm_end