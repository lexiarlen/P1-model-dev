#!/usr/bin/env bash
#SBATCH --job-name=get_last_dump
#SBATCH --output=../out_files/%x_%A.out
#SBATCH --error=../out_files/%x_%A.err
#SBATCH --time=00:10:00
#SBATCH -p serc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12

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
PROJECT_DIR="${PROJECT_DIR:-/home/groups/earlew/arlenlex/P1-model-dev}"
SCRATCH_PROJECT_DIR="${SCRATCH_PROJECT_DIR:-/scratch/groups/earlew/arlenlex/P1-model-dev}"
MODULEFILES_DIR="${MODULEFILES_DIR:-/home/groups/earlew/arlenlex/.local/modulefiles}"
CONDA_SH="${CONDA_SH:-/home/groups/earlew/arlenlex/miniconda3-local/etc/profile.d/conda.sh}"


ml purge
source ${CONDA_SH}
conda activate postprocess_npy
export PATH="$CONDA_PREFIX/bin:$PATH"

export DATA_DIR="${DATA_DIR:-${SCRATCH_PROJECT_DIR}/simple_shear/size_effect/ensemble10_smno}"

python3 ${PROJECT_DIR}/simple_shear/size_effect/analysis/get_last_dump.py 
