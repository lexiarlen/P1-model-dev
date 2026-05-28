#!/usr/bin/env bash
#SBATCH --job-name=get_last_dump
#SBATCH --output=/home/groups/earlew/arlenlex/P1-model-dev/simple_shear/size_effect/out_files/%x_%A.out
#SBATCH --error=/home/groups/earlew/arlenlex/P1-model-dev/simple_shear/size_effect/out_files/%x_%A.err
#SBATCH --time=00:10:00
#SBATCH -p serc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12

ml purge
source /home/groups/earlew/arlenlex/miniconda3-local/etc/profile.d/conda.sh
conda activate postprocess_npy
export PATH="$CONDA_PREFIX/bin:$PATH"

python3 /home/groups/earlew/arlenlex/P1-model-dev/simple_shear/size_effect/analysis/get_last_dump.py 
