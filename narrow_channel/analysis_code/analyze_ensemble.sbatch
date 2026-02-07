#!/usr/bin/env bash
#SBATCH --job-name=analyze_job
#SBATCH --output=/scratch/groups/earlew/arlenlex/P1-model-dev/nares/out_files/analysis/%x_%A.out
#SBATCH --error=/scratch/groups/earlew/arlenlex/P1-model-dev/nares/out_files/analysis/%x_%A.err
#SBATCH --time=01:00:00
#SBATCH -p serc
#SBATCH -c 12


echo "Starting batch script analyze_ensemble.sbatch"

set -euo pipefail

: "${J:?J is required}"
: "${max_wind_speed:?max_wind_speed is required}"
: "${parent_dir_name:?parent_dir_name is required}"
: "${rho_a:?rho_a is required}"
: "${C_a:?C_a is required}"
: "${T_max:?T_max is required}"
: "${radius:?radius is required}"
: "${rho_i:?rho_i is required}"

ml purge
source /home/groups/earlew/arlenlex/miniconda3-local/etc/profile.d/conda.sh
conda activate ovito
export PATH="$CONDA_PREFIX/bin:$PATH"

python -u /home/groups/earlew/arlenlex/P1-model-dev/nares/analysis_code/analyze_ensemble.py \
  --J "${J}" \
  --max_wind_speed "${max_wind_speed}" \
  --parent_dir_name "${parent_dir_name}" \
  --rho_a "${rho_a}" \
  --C_a "${C_a}" \
  --T_max "${T_max}" \
  --radius "${radius}" \
  --rho_i "${rho_i}" \
