#!/usr/bin/env bash
#SBATCH --job-name=size_effect_tests
#SBATCH --output=../out_files/%x_%A.out
#SBATCH --error=../out_files/%x_%A.err
#SBATCH --time=72:00:00
#SBATCH -p serc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=9

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


# 1) Run LAMMPS for each domain size and granularity
Ls=(50 60 80 100 200 300 400 500 600 700 800 900) 
num_atoms=(000873 001257 002237 003496 013987 031475 055954 087431 125900 171366 223824 283279) 
ncores=(1 1 1 2 2 4 4 9 9 9 9 9)  
seeds2=(8302 3094 3894 1092 4999 5673 1927 2380 8095 8920)
dts=(0.1 0.1 0.1 0.075 0.025 0.01 0.01 0.01 0.005 0.005 0.005 0.005)
run_time=1800
num_dumps=500

for idx in "${!Ls[@]}"; do
    L=${Ls[$idx]}
    N=${num_atoms[$idx]}
    NP=${ncores[$idx]}
    DT=${dts[$idx]}
    for s in "${seeds2[@]}"; do
        mpirun -oversubscribe -n $NP lmp \
            -var L_pack "$L" \
            -var natoms  "$N" \
            -var dt "$DT" \
            -var atomfile_dir "${SCRATCH_PROJECT_DIR}/packing/size-effect/ensemble10/L${L}/atomfiles_s${s}/" \
            -var out_dir "${SCRATCH_PROJECT_DIR}/simple_shear/size_effect/ensemble10_smno/L${L}/" \
            -var outfile_name "N${N}Dump_seed${s}.lammps" \
            -var run_time "$run_time" \
            -var num_dumps "$num_dumps" \
            -var knb "1e9" \
            -var eps_crit "0.00005" \
            -in  ${PROJECT_DIR}/simple_shear/in.shear \
            -log none
    done
done

# 2) Switch to ovito conda environment to postprocess
ml purge
source ${CONDA_SH}
conda activate ovito
export PATH="$CONDA_PREFIX/bin:$PATH"

# Create numpy arrays of stresses using Ovito (pandas incompatible with ovito on Sherlock)
python3 ${PROJECT_DIR}/simple_shear/size_effect/analysis/compute_macroscopic_stress.py \
    "${SCRATCH_PROJECT_DIR}/simple_shear/size_effect/ensemble10_smno/" 

# 3) Turn numpy arrays into parquet for visualization
ml purge
source ${CONDA_SH}
conda activate postprocess_npy
export PATH="$CONDA_PREFIX/bin:$PATH"

python3 ${PROJECT_DIR}/simple_shear/size_effect/analysis/numpy_to_nc_of_stresses.py \
    "${SCRATCH_PROJECT_DIR}/simple_shear/size_effect/ensemble10_smno/" 
