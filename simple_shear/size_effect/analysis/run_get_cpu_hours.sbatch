#!/usr/bin/env bash
#SBATCH --job-name=parse_size_effect
#SBATCH --output=../out_files/parse_size_effect_%j.out
#SBATCH --error=../out_files/parse_size_effect_%j.err
#SBATCH --time=00:5:00
#SBATCH -p serc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

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


# ---------- environment ----------
ml purge
source ${CONDA_SH}
conda activate ovito

# ---------- paths ----------
script_path="${PROJECT_DIR}/simple_shear/size_effect/analysis/get_cpuhours.py"

base_dir="${PROJECT_DIR}/simple_shear/size_effect/out_files_smno"

out_dir="${SCRATCH_PROJECT_DIR}/simple_shear/size_effect/ensemble10_smno"
out_path="${out_dir}/runs_dict.json"

# ---------- files to parse ----------
files=(
    "size_effect_tests_22983281.out"
)

# ---------- run parser ----------
python "$script_path" \
    --base-dir "$base_dir" \
    --files "${files[@]}" \
    --out-path "$out_path"