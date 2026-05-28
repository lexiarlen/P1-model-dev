#!/usr/bin/env bash
#SBATCH --job-name=analyze_job
#SBATCH --output=../out_files/%x_%A.out
#SBATCH --error=../out_files/%x_%A.err
#SBATCH --time=01:10:00
#SBATCH -p serc
#SBATCH -c 16


echo "Starting batch script analyze_ensemble.sbatch"

set -euo pipefail

# Resolve paths without hard-coding a user-specific checkout location.
find_project_dir() {
    local dir="$SCRIPT_DIR"
    while [[ "$dir" != "/" ]]; do
        if [[ -d "$dir/narrow_channel" && -d "$dir/packing" && -d "$dir/simple_shear" ]]; then
            printf '%s\n' "$dir"
            return 0
        fi
        dir="$(dirname "$dir")"
    done
    printf '%s\n' "$SCRIPT_DIR"
}
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${PROJECT_DIR:-$(find_project_dir)}"
SCRATCH_PROJECT_DIR="${SCRATCH_PROJECT_DIR:-${PROJECT_DIR}}"
MODULEFILES_DIR="${MODULEFILES_DIR:-/home/groups/earlew/arlenlex/.local/modulefiles}"
CONDA_SH="${CONDA_SH:-/home/groups/earlew/arlenlex/miniconda3-local/etc/profile.d/conda.sh}"


: "${J:?J is required}"
: "${base_dir:?base_dir is required}"
: "${state_dir:?state_dir is required}"

# 1) analyze individual runs
ml purge
source ${CONDA_SH}
conda activate ovito
export PATH="$CONDA_PREFIX/bin:$PATH"

python -u ${PROJECT_DIR}/narrow_channel/parameter_sweep/analysis/analyze_ensemble.py --base_dir "${base_dir}" --J "${J}"

# 2) npy to parquet
ml purge
source ${CONDA_SH}
conda activate postprocess_npy
export PATH="$CONDA_PREFIX/bin:$PATH"

python ${PROJECT_DIR}/narrow_channel/parameter_sweep/analysis/np_to_pd.py  --state_dir "${state_dir}"

