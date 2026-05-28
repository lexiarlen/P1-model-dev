#!/usr/bin/env bash
#SBATCH --job-name=kc_ensemble
#SBATCH --output=../out_files/%x_%A_%a.out
#SBATCH --error=../out_files/%x_%A_%a.err
#SBATCH --time=10:00:00            
#SBATCH -p serc
#SBATCH --array=0-4%5
#SBATCH --nodes=1
#SBATCH --ntasks=3
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
ml use ${MODULEFILES_DIR}
ml load lammps_lexi chemistry

export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

# ---------- fixed params ----------
kb=1e9
eps=0.00002

# ---------- kc sweep ----------
kc_vals=(1e5 1e6 1e7 1e8 1e9)

# ---------- dirs ----------
packing_dir="${SCRATCH_PROJECT_DIR}/packing/narrow_channel/atomfiles_r400/"
base_dir="${SCRATCH_PROJECT_DIR}/narrow_channel_smno/kc_ensemble_small"
state_dir="${base_dir}/state"
output_dir="${base_dir}/lammps_output"

mkdir -p "$state_dir"
mkdir -p "$output_dir"

# ---------- current array task ----------
idx="${SLURM_ARRAY_TASK_ID}"
kc="${kc_vals[$idx]}"
MEMBER=$((idx + 1))

# ---------- sanity checks ----------
: "${kc:?missing kc for idx=${idx}}"

echo "RUNNING idx=${idx} member=${MEMBER} ranks=${SLURM_NTASKS}"
echo "kb=${kb}"
echo "eps=${eps}"
echo "kc=${kc}"
echo "output_dir=${output_dir}"
echo "packing_dir=${packing_dir}"

# ---------- optional parameter table ----------
PARAMS="${state_dir}/theta.tsv"

if [[ "$idx" -eq 0 ]]; then
    echo -e "member\tkb\teps\tkc" > "$PARAMS"

    for i in "${!kc_vals[@]}"; do
        m=$((i + 1))
        echo -e "${m}\t${kb}\t${eps}\t${kc_vals[$i]}" >> "$PARAMS"
    done
fi

# ---------- run LAMMPS ----------
mpirun -oversubscribe -np "${SLURM_NTASKS}" lmp \
    -in ${PROJECT_DIR}/narrow_channel/lammps_scripts/in.strait \
    -var knb "${kb}" \
    -var eps_crit "${eps}" \
    -var knc "${kc}" \
    -var member "${MEMBER}" \
    -var outdir "${output_dir}/" \
    -var atomfile_dir "${packing_dir}" \
    -log none