#!/bin/bash
#SBATCH -J rfdaa
#SBATCH -c 4
#SBATCH -q normal
#SBATCH -p gpu-farm
#SBATCH --gres=gpu:1
#SBATCH --array=0-7
source /home/hwjang/miniforge3/bin/activate rfdaa
source /home/hwjang/.slurmrc
slurm_start
###############################################
export RFDIFFUSION_PATH="/home/hwjang/aipd/rf_diffusion_all_atom"
export N_PER_ARRAY=8
python $RFDIFFUSION_PATH/run_inference.py \
    inference.ckpt_path=$RFDIFFUSION_PATH/RFDiffusionAA_paper_weights.pt \
    inference.deterministic=False \
    diffuser.T=30 \
    inference.output_prefix=output/result \
    inference.input_pdb=../0_params/PHT.pdb \
    contigmap.contigs=[\'100-130\'] \
    inference.ligand=PHT \
    inference.num_designs=$N_PER_ARRAY \
    inference.design_startnum=$((SLURM_ARRAY_TASK_ID*N_PER_ARRAY))
###############################################
slurm_end