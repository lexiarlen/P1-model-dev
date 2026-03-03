#!/usr/bin/env bash
#SBATCH --job-name=begin_ensemble
#SBATCH --output=/home/groups/earlew/arlenlex/P1-model-dev/narrow_channel/out_files/%x_%A.out
#SBATCH --error=/home/groups/earlew/arlenlex/P1-model-dev/narrow_channel/out_files/%x_%A.err
#SBATCH --time=00:10:00
#SBATCH -p serc
#SBATCH -c 1

set -euo pipefail

ml purge 
source /home/groups/earlew/arlenlex/miniconda3-local/etc/profile.d/conda.sh
conda activate ovito
export PATH="$CONDA_PREFIX/bin:$PATH"

J=200               # number of samples
natoms=52457
radius=400

packing_dir="/scratch/groups/earlew/arlenlex/P1-model-dev/packing/narrow_channel/atomfiles_r400/" 
base_dir="/scratch/groups/earlew/arlenlex/P1-model-dev/narrow_channel"
state_dir="${base_dir}/state"
output_dir="${base_dir}/lammps_output/"

mkdir -p "$state_dir"
mkdir -p "$output_dir"

# sample parameters
python /home/groups/earlew/arlenlex/P1-model-dev/narrow_channel/analysis/loghypercube_parameter_sampling.py \
    --kb-min 1e8  --kb-max 1e10 \
    --eps-min 1e-6 --eps-max 1e-4 \
    --J "${J}" --seed 123 \
    --out "${state_dir}/theta.tsv"

PARAMS="${state_dir}/theta.tsv"

#submit the model array for this u
ARRAY_JOBID=$(sbatch --parsable \
    --array=1-"${J}"%"${J}" \
    --export=ALL,PARAMS="${PARAMS}",J="${J}",packing_dir="${packing_dir}",output_dir="${output_dir}" \
    /home/groups/earlew/arlenlex/P1-model-dev/narrow_channel/analysis/run_lammps_ensemble.sbatch)

#submit the analysis job for this u
ANALYZE_JOBID=$(sbatch --parsable \
    --dependency=afterok:${ARRAY_JOBID} \
    --export=ALL,J="${J}",ARRAY_JOBID="${ARRAY_JOBID}",base_dir="${base_dir}" \
    /home/groups/earlew/arlenlex/P1-model-dev/narrow_channel/analysis/analyze_ensemble.sbatch)


