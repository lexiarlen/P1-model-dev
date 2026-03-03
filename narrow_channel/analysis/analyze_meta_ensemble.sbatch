#!/usr/bin/env bash
#SBATCH --job-name=meta_analyze_job
#SBATCH --output=/home/groups/earlew/arlenlex/P1-EKI-nares/EKI/out_files/%x_%A.out
#SBATCH --error=/home/groups/earlew/arlenlex/P1-EKI-nares/EKI/out_files/%x_%A.err
#SBATCH --time=00:20:00        # quick analysis
#SBATCH -p serc

echo "Starting batch script analyze_meta_ensemble.sbatch" 

set -euo pipefail

: "${J:?}"
: "${run_name:?}"
: "${max_wind_speed:?}"
: "${parent_dir_name:?parent_dir_name is required}"
: "${power:?power is required}"
: "${var:?var is required}"
: "${natoms:?natoms is required}"

ml purge
source /home/groups/earlew/arlenlex/miniconda3-local/etc/profile.d/conda.sh
conda activate nextsim
export PATH="$CONDA_PREFIX/bin:$PATH"

# Compute y, update params,
python /home/groups/earlew/arlenlex/P1-model-dev/ensemble_runs/analysis_code/meta_analyze_ensemble4.py --J "${J}" --parent_dir_name "${parent_dir_name}" --power="${power}" --variable="${var}" --natoms="${natoms}"
