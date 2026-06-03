#!/bin/bash
#SBATCH -J params
#SBATCH -c 1
#SBATCH -q normal
#SBATCH -p cpu-farm
source /home/hwjang/miniforge3/bin/activate rfdaa
source /home/hwjang/.slurmrc
slurm_start
###############################################
export ROSETTA_PATH="/home/hwjang/aipd/rosetta"
export PYTHONPATH="$PYTHONPATH:$ROSETTA_PATH/source/scripts/python/public"

export INPUT="Conformer3D_COMPOUND_CID_1775.sdf"
export RES_NAME="PHT"

pymol -cq pymol_script.py -- $INPUT $RES_NAME.pdb $RES_NAME.mol2
(cat << EOF
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N   
ATOM      2  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C   
ATOM      3  C   ALA A   1       0.000   0.000   0.000  1.00  0.00           C   
ATOM      4  O   ALA A   1       0.000   0.000   0.000  1.00  0.00           O   
ATOM      5  CB  ALA A   1       0.000   0.000   0.000  1.00  0.00           C   
TER
EOF
cat $RES_NAME.pdb) > tmp && mv tmp $RES_NAME.pdb

$ROSETTA_PATH/source/scripts/python/public/generic_potential/mol2genparams.py --input $RES_NAME.mol2 --resname $RES_NAME

###############################################
slurm_end