#!/usr/bin/env bash
#SBATCH --job-name=kc_checks
#SBATCH --output=../out_files/%x_%A_%a.out
#SBATCH --error=../out_files/%x_%A_%a.err
#SBATCH --time=36:00:00
#SBATCH -p serc
#SBATCH --array=0-3%4
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
PROJECT_DIR="${PROJECT_DIR:-/home/groups/earlew/arlenlex/P1-model-dev}"
SCRATCH_PROJECT_DIR="${SCRATCH_PROJECT_DIR:-/scratch/groups/earlew/arlenlex/P1-model-dev}"
MODULEFILES_DIR="${MODULEFILES_DIR:-/home/groups/earlew/arlenlex/.local/modulefiles}"
CONDA_SH="${CONDA_SH:-/home/groups/earlew/arlenlex/miniconda3-local/etc/profile.d/conda.sh}"


# ---------- environment ----------
ml purge
ml use ${MODULEFILES_DIR}
ml load lammps_lexi chemistry

export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

# ---------- parameter arrays ----------
kbs=(175340990 998969920 5615634000 3217549500)
eps_vals=(5.9494822e-7 1.1519924e-5 7.19289e-6 2.3013252e-5)

# ---------- fixed kc ----------
kc_val=1e9

# ---------- dirs ----------
packing_dir="${SCRATCH_PROJECT_DIR}/packing/narrow_channel/atomfiles_r400/"
base_dir="${SCRATCH_PROJECT_DIR}/narrow_channel_smno/kc1e9_regime_checks_higher_damp/"
state_dir="${base_dir}/state"
output_dir="${base_dir}/lammps_output"

mkdir -p "$state_dir"
mkdir -p "$output_dir"

# ---------- current array task ----------
idx="${SLURM_ARRAY_TASK_ID}"

kb="${kbs[$idx]}"
eps="${eps_vals[$idx]}"
kc="${kc_val}"

MEMBER=$((idx + 1))

# ---------- sanity checks ----------
: "${kb:?missing kb for idx=${idx}}"
: "${eps:?missing eps for idx=${idx}}"
: "${kc:?missing kc for idx=${idx}}"

echo "RUNNING idx=${idx} member=${MEMBER} ranks=${SLURM_NTASKS}"
echo "kb=${kb}"
echo "eps=${eps}"
echo "kc=${kc}"
echo "output_dir=${output_dir}"
echo "packing_dir=${packing_dir}"


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