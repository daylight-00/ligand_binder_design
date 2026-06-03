#!/bin/bash
#SBATCH -J lmpnn
#SBATCH -c 4
#SBATCH -q normal
#SBATCH -p gpu-farm
#SBATCH --gres=gpu:1
#SBATCH --array=0-3
source /home/hwjang/miniforge3/bin/activate rfdaa
source /home/hwjang/.slurmrc
slurm_start
###############################################
export LIGAND_MPNN_PATH="/home/hwjang/aipd/LigandMPNN"
python $LIGAND_MPNN_PATH/run.py \
    --seed 111 \
    --model_type "ligand_mpnn" \
    --out_folder "output" \
    --pdb_path_multi ../2_backbone_filter/lmpnn_input/split_$SLURM_ARRAY_TASK_ID.json \
    --number_of_batches 1 \
    --batch_size 8 \
    --pack_side_chains 1 \
    --number_of_packs_per_design 1 \
    --pack_with_ligand_context 1 \
    --temperature 0.2
###############################################
slurm_end