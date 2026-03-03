#!/usr/bin/env bash
#SBATCH --job-name=analyze_job
#SBATCH --output=/scratch/groups/earlew/arlenlex/P1-model-dev/narrow_channel/out_files/analysis/%x_%A.out
#SBATCH --error=/scratch/groups/earlew/arlenlex/P1-model-dev/narrow_channel/out_files/analysis/%x_%A.err
#SBATCH --time=01:00:00
#SBATCH -p serc
#SBATCH -c 12


echo "Starting batch script analyze_ensemble.sbatch"

set -euo pipefail

: "${J:?J is required}"
: "${base_dir:?base_dir is required}"

# 1) analyze individual runs
ml purge
source /home/groups/earlew/arlenlex/miniconda3-local/etc/profile.d/conda.sh
conda activate ovito
export PATH="$CONDA_PREFIX/bin:$PATH"

python -u /home/groups/earlew/arlenlex/P1-model-dev/narrow_channel/analysis/analyze_ensemble.py --base_dir "${base_dir}" --J "${J}"

# 2) npy to parquet
ml purge
source /home/groups/earlew/arlenlex/miniconda3-local/etc/profile.d/conda.sh
conda activate nextsim
export PATH="$CONDA_PREFIX/bin:$PATH"

python /home/groups/earlew/arlenlex/P1-model-dev/ensemble_runs/analysis/np_to_pd.py --J "${J}" --base_dir "${base_dir}"

