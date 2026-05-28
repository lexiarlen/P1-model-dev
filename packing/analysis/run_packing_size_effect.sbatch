#!/usr/bin/env bash
#SBATCH --job-name=making_packings_for_size_effect_ensemble_runs
#SBATCH --output=../out_files/%x_%A.out
#SBATCH --error=../out_files/%x_%A.err
#SBATCH --time=24:00:00
#SBATCH -p serc
#SBATCH -c 2

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


module purge
module use ${MODULEFILES_DIR}
module load lammps_lexi chemistry

# 1) create 10 differently seeded packings for each domain size = 100e3/L
Ls=(50 60 80 100 200 300 400 500 600 700 800 900)
seeds1=(9834 1825 2893 2283 1923 1289 8923 6573 9023 1098)               
seeds2=(8302 3094 3894 1092 4999 5673 1927 2380 8095 8920)
for L in "${Ls[@]}"; do
    for idx in "${!seeds1[@]}"; do
        s1=${seeds1[$idx]}
        s2=${seeds2[$idx]}
        mpirun -oversubscribe -np 2 \
        lmp -var H "${L}" -var W "${L}" \
        -var seed1 "${s1}" -var seed2 "${s2}" \
        -var outdir "${SCRATCH_PROJECT_DIR}/packing/size-effect/ensemble10/L${L}/" \
        -var outfile_name "datafile_s${s2}.lammps" \
        -in ${PROJECT_DIR}/packing/lammps_scripts/in.pack \
        -log none
    done
done

# 2) generate atomfiles with the python script
ml purge
source ${CONDA_SH}
conda activate nextsim
export PATH="$CONDA_PREFIX/bin:$PATH"

true_L=100000
for L in "${Ls[@]}"; do
    for s2 in "${seeds2[@]}"; do
        radius_scaling=$(bc -l <<< "$true_L / $L")

        outdir="${SCRATCH_PROJECT_DIR}/packing/size-effect/ensemble10/L${L}/atomfiles_s${s2}"
        mkdir -p "$outdir"

        python3 create_atomfiles.py \
        "${SCRATCH_PROJECT_DIR}/packing/size-effect/ensemble10/L${L}/datafile_s${s2}.lammps" \
        "$radius_scaling" \
        "$outdir"
    done
done
