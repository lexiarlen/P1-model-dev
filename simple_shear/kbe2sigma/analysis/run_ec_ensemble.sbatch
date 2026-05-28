#!/usr/bin/env bash
#SBATCH --job-name=kb2elastic_ensemble
#SBATCH --output=../out_files/%x_%A.out
#SBATCH --error=../out_files/%x_%A.err
#SBATCH --time=04:00:00
#SBATCH -p serc
#SBATCH --ntasks=2
#SBATCH --nodes=1

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


module purge
module use ${MODULEFILES_DIR}
module load lammps_lexi chemistry

# --- ensemble config ---
L=300
num_atoms=031475
s=8302

kb=1e9
dt=0.025
eps=(5.00e-06 6.00e-05 1.15e-04 1.70e-04 2.25e-04 2.80e-04 3.35e-04 3.90e-04 4.45e-04 5.00e-04) 
base_out="${SCRATCH_PROJECT_DIR}/simple_shear/kbe2sigma_smno/L${L}"
atomfile_dir="${SCRATCH_PROJECT_DIR}/packing/size-effect/ensemble10/L${L}/atomfiles_s${s}/"
run_time="5400"
num_dumps="500"

mkdir -p "${base_out}"

for idx in "${!eps[@]}"; do
    ep="${eps[$idx]}"

    case_dir="${base_out}/eps_${ep}"
    mkdir -p "${case_dir}"

    outfile="N${num_atoms}_s${s}_eps${ep}.lammps"

    echo "[RUN] kb=${kb} eps=${ep} dt=${dt} -> ${case_dir}/${outfile}"

    mpirun -n 2 lmp \
        -var L_pack "${L}" \
        -var natoms  "${num_atoms}" \
        -var knb "${kb}" \
        -var eps_crit "${ep}" \
        -var atomfile_dir "${atomfile_dir}" \
        -var out_dir "${case_dir}/" \
        -var outfile_name "${outfile}" \
        -var run_time "${run_time}" \
        -var num_dumps "${num_dumps}" \
        -var dt "${dt}" \
        -in ${PROJECT_DIR}/simple_shear/in.shear \
        -log none
done

# 1) create numpy files for each ensemble member
ml purge
source ${CONDA_SH}
conda activate ovito
export PATH="$CONDA_PREFIX/bin:$PATH"

python3 ${PROJECT_DIR}/simple_shear/kbe2sigma/analysis/compute_max_strengths.py --base-dir "${base_out}" --qc-plots

# 2) convert numpy files to parquet
ml purge
source ${CONDA_SH}
conda activate postprocess_npy
export PATH="$CONDA_PREFIX/bin:$PATH"

python3 ${PROJECT_DIR}/simple_shear/kbe2sigma/analysis/numpy_to_pd.py --base-dir "${base_out}"
