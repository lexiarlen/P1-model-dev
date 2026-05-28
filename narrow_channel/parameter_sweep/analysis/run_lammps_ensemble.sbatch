#!/usr/bin/env bash
#SBATCH --job-name=run_lammps
#SBATCH --output=../out_files/%x_%A_%a.out
#SBATCH --error=../out_files/%x_%A_%a.err
#SBATCH --time=44:00:00            
#SBATCH -p serc
#SBATCH --array=1-1%50              # launcher will override with --array=1-${J}%${J}
#SBATCH --nodes=1
#SBATCH --ntasks=3                 # 3 MPI ranks per ensemble member
#SBATCH --cpus-per-task=1          # 1 CPU per MPI rank

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


# ---------- required inputs via --export from launcher ----------
: "${PARAMS:?PARAMS (theta.tsv) path is required}"
: "${J:?J (ensemble size) is required}"
: "${output_dir:?output_dir is required}"
: "${packing_dir:?packing_dir is required}"

# ---------- environment ----------
ml purge
ml use ${MODULEFILES_DIR}
ml load lammps_lexi chemistry

# pure MPI: no OpenMP threads
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

# ---------- current ensemble member (1-based) ----------
MEMBER="${SLURM_ARRAY_TASK_ID}"

# ---------- read the MEMBER-th nonblank, noncomment line as "kb eps" ----------
read -r kb eps < <(awk 'NF && $1 !~ /^#/' "${PARAMS}" | sed -n "${MEMBER}p")

# ---------- sanity checks ----------
: "${kb:?missing kb for MEMBER=${MEMBER}}"
: "${eps:?missing eps for MEMBER=${MEMBER}}"

echo "RUNNING member=${MEMBER} ranks=${SLURM_NTASKS} kb=${kb} eps=${eps}"

# ---------- run LAMMPS on 3 MPI ranks ----------
mpirun -oversubscribe -np "${SLURM_NTASKS}" lmp \
    -in ${PROJECT_DIR}/narrow_channel/lammps_scripts/in.strait \
    -var knb "${kb}" \
    -var eps_crit "${eps}" \
    -var member "${MEMBER}" \
    -var outdir "${output_dir}" \
    -var atomfile_dir "${packing_dir}" \
    -log none
