#!/usr/bin/env bash
#SBATCH --job-name=meta_analyze_job
#SBATCH --output=/home/groups/earlew/arlenlex/P1-model-dev/nares/out_files/analysis/%x_%A.out
#SBATCH --error=/home/groups/earlew/arlenlex/P1-model-dev/nares/out_files/analysis/%x_%A.err
#SBATCH --time=00:12:00
#SBATCH -p serc

echo "Starting batch script analyze_ensemble.sbatch"

set -euo pipefail

ml purge
source /home/groups/earlew/arlenlex/miniconda3-local/etc/profile.d/conda.sh
conda activate ovito
export PATH="$CONDA_PREFIX/bin:$PATH"

echo "analyze running"


python -u /home/groups/earlew/arlenlex/P1-model-dev/nares/analysis_code/analyze_ensemble.py \
  --J "200" \
  --max_wind_speed "20" \
  --parent_dir_name "test" \
  --rho_a "1" \
  --C_a "0.0012" \
  --T_max "36" \
  --radius "500" \
  --rho_i "920" \

# echo "running meta analyze"

# ml purge
# source /home/groups/earlew/arlenlex/miniconda3-local/etc/profile.d/conda.sh
# conda activate nextsim
# export PATH="$CONDA_PREFIX/bin:$PATH"

# python -u /home/groups/earlew/arlenlex/P1-EKI-nares/EKI/meta_analyze_ensemble_robust.py --J "200" --parent_dir_name "domain_scalings/quadratic" --power "2" --variable "domain_scaling" --natoms "2936"



