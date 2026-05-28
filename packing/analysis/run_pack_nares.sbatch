#!/usr/bin/env bash
#SBATCH --job-name=pack_nares
#SBATCH --output=../out_files/%x_%A.out
#SBATCH --error=../out_files/%x_%A.err
#SBATCH --time=00:30:00
#SBATCH -p serc
#SBATCH -c 1

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


# 1) create packing for narrow channel simulations with radius = 400 m => granularity <= 0.033
module purge
module use ${MODULEFILES_DIR}
module load lammps_lexi chemistry

EXEC="lmp"

mpirun -oversubscribe -np 1 \
      $EXEC -var H "600" -var W "350" \
      -var outdir "${SCRATCH_PROJECT_DIR}/packing/narrow_channel/" \
      -var outfile_name "datafile_r400.lammps" \
      -in ${PROJECT_DIR}/packing/lammps_scripts/in.pack -log none

# 2) Switch to your conda environment and generate atomfiles
ml purge
source ${CONDA_SH}
conda activate nextsim
export PATH="$CONDA_PREFIX/bin:$PATH"

radius_scaling=400

outdir="${SCRATCH_PROJECT_DIR}/packing/narrow_channel/atomfiles_r400/"
mkdir -p "$outdir"

python3 create_atomfiles.py \
    "${SCRATCH_PROJECT_DIR}/packing/narrow_channel/datafile_r400.lammps" \
    "$radius_scaling" \
    "$outdir"
