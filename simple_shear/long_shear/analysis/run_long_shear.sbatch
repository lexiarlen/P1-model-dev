#!/usr/bin/env bash
#SBATCH --job-name=longer_shear
#SBATCH --output=../out_files/%x_%A.out
#SBATCH --error=../out_files/%x_%A.err
#SBATCH --time=40:00:00
#SBATCH -p serc
#SBATCH -c 4

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


module purge
module use ${MODULEFILES_DIR}
module load lammps_lexi chemistry

s=8095
L=300
N=031475

mpirun -oversubscribe -n 4 lmp \
        -var L_pack "$L" \
        -var natoms  "$N" \
        -var atomfile_dir "${SCRATCH_PROJECT_DIR}/packing/size-effect/ensemble10/L${L}/atomfiles_s${s}/" \
        -var out_dir "${SCRATCH_PROJECT_DIR}/simple_shear/long_shear/" \
        -var outfile_name "N${N}Dump_seed${s}_1000hours_smno.lammps" \
        -var run_time "3600000" \
        -var num_dumps "1000" \
        -var dt "0.025" \
        -var knb "1e9" \
        -var eps_crit "0.00005" \
        -in  ${PROJECT_DIR}/simple_shear/in.shear \
        -log none