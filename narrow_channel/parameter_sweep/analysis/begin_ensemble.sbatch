#!/usr/bin/env bash
#SBATCH --job-name=begin_ensemble
#SBATCH --output=../out_files/%x_%A.out
#SBATCH --error=../out_files/%x_%A.err
#SBATCH --time=00:2:00
#SBATCH -p serc
#SBATCH -c 1

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


OUTDIR="${SCRATCH_PROJECT_DIR}/narrow_channel/out_files/lammps_jobs"

if [[ -d "$OUTDIR" ]]; then
    echo "Cleaning $OUTDIR"
    rm -rf -- "${OUTDIR:?}/"*
fi

ml purge 
source ${CONDA_SH}
conda activate ovito
export PATH="$CONDA_PREFIX/bin:$PATH"

J=200               # number of samples
natoms=73442
radius=400

packing_dir="${SCRATCH_PROJECT_DIR}/packing/narrow_channel/atomfiles_r400/" 
base_dir="${SCRATCH_PROJECT_DIR}/narrow_channel_smno"
state_dir="${base_dir}/state"
output_dir="${base_dir}/lammps_output/"

mkdir -p "$state_dir"
mkdir -p "$output_dir"

# sample parameters
python ${PROJECT_DIR}/narrow_channel/parameter_sweep/analysis/loghypercube_parameter_sampling.py \
    --kb-min 1e8  --kb-max 1e10 \
    --eps-min 1e-7 --eps-max 5e-5 \
    --J "${J}" --seed 123 \
    --out "${state_dir}/theta.tsv"

PARAMS="${state_dir}/theta.tsv"

#submit the model array for this u
ARRAY_JOBID=$(sbatch --parsable \
    --array=1-"${J}"%"${J}" \
    --export=ALL,PARAMS="${PARAMS}",J="${J}",packing_dir="${packing_dir}",output_dir="${output_dir}" \
    ${PROJECT_DIR}/narrow_channel/parameter_sweep/analysis/run_lammps_ensemble.sbatch)

#submit the analysis job for this u
ANALYZE_JOBID=$(sbatch --parsable \
    --dependency=afterok:${ARRAY_JOBID} \
    --export=ALL,J="${J}",ARRAY_JOBID="${ARRAY_JOBID}",base_dir="${base_dir}",state_dir="${state_dir}" \
    ${PROJECT_DIR}/narrow_channel/parameter_sweep/analysis/analyze_ensemble.sbatch)