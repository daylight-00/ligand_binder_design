#!/bin/bash
#SBATCH -J af3_inf
#SBATCH -c 4
#SBATCH -q normal
#SBATCH -p gpu-farm
#SBATCH --gres=gpu:1
#SBATCH --array=0-3
source /home/hwjang/miniforge3/bin/activate af3
source /home/hwjang/.slurmrc
slurm_start
###############################################
export XLA_FLAGS="--xla_gpu_enable_triton_gemm=false"
export XLA_PYTHON_CLIENT_PREALLOCATE=true
export XLA_CLIENT_MEM_FRACTION=0.95
export ALPHAFOLD3_PATH="/home/hwjang/project/alphafold3"

python $ALPHAFOLD3_PATH/run_alphafold.py \
    --run_inference=true \
    --run_data_pipeline=false \
    --input_dir "json/$SLURM_ARRAY_TASK_ID" \
    --output_dir "output" \
    # --save_embeddings=true \
###############################################
slurm_end