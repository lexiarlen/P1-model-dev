#!/usr/bin/env bash
#SBATCH --job-name=begin_ensemble
#SBATCH --output=/home/groups/earlew/arlenlex/P1-model-dev/nares/out_files/%x_%A.out
#SBATCH --error=/home/groups/earlew/arlenlex/P1-model-dev/nares/out_files/%x_%A.err
#SBATCH --time=00:10:00
#SBATCH -p serc
#SBATCH -c 1

set -euo pipefail

ml purge 
source /home/groups/earlew/arlenlex/miniconda3-local/etc/profile.d/conda.sh
conda activate ovito
export PATH="$CONDA_PREFIX/bin:$PATH"

J=200
ell=1
rho_i=920
rho_a=1
power=1
rho_o=1027
parent_dir_name="baseline_correct_geometry"
u=20
T_max=36
C_a=0.0012
natoms=47001
radius=500
packing_dir="/scratch/groups/earlew/arlenlex/P1-model-dev/packing/atomfiles_test_r${radius}/" # move this...
damping_ratio=1
C=40e3

ANALYZE_IDS=()

RUN_STATE="/scratch/groups/earlew/arlenlex/P1-model-dev/nares/${parent_dir_name}/state"
# mkdir -p "$RUN_STATE"

# theta for this run
# python /home/groups/earlew/arlenlex/P1-model-dev/nares/analysis_code/make_theta0_loghypercube2.py \
#     --kb-min 1e8  --kb-max 1e10 \
#     --eps-min 1e-6 --eps-max 1e-4 \
#     --J "${J}" --seed 123 \
#     --out "${RUN_STATE}/theta.tsv"

PARAMS="${RUN_STATE}/theta.tsv"

# #submit the model array for this u
# ARRAY_JOBID=$(sbatch --parsable \
#     --array=1-"${J}"%"${J}" \
#     --export=ALL,PARAMS="${PARAMS}",J="${J}",T_max="${T_max}",C="${C}",natoms="${natoms}",radius="${radius}",packing_dir="${packing_dir}",domain_scaling="${ell}",damping_ratio="${damping_ratio}",power="${power}",max_wind_speed="${u}",rho_i="${rho_i}",rho_a="${rho_a}",rho_o="${rho_o}",C_a="${C_a}",parent_dir_name="${parent_dir_name}" \
#     /home/groups/earlew/arlenlex/P1-model-dev/nares/analysis_code/run_lammps_general.sbatch)

# #submit the analysis job for this u
# ANALYZE_JOBID=$(sbatch --parsable \
#     --dependency=afterok:${ARRAY_JOBID} \
#     --export=ALL,J="${J}",T_max="${T_max}",ARRAY_JOBID="${ARRAY_JOBID}",radius="${radius}",max_wind_speed="${u}",rho_i="${rho_i}",rho_a="${rho_a}",domain_scaling="${ell}",C_a="${C_a}",rho_o="${rho_o}",parent_dir_name="${parent_dir_name}" \
#     /home/groups/earlew/arlenlex/P1-model-dev/nares/analysis_code/analyze_ensemble.sbatch)

ANALYZE_JOBID=$(sbatch --parsable \
    --export=ALL,J="${J}",T_max="${T_max}",radius="${radius}",max_wind_speed="${u}",rho_i="${rho_i}",rho_a="${rho_a}",domain_scaling="${ell}",C_a="${C_a}",rho_o="${rho_o}",parent_dir_name="${parent_dir_name}" \
    /home/groups/earlew/arlenlex/P1-model-dev/nares/analysis_code/analyze_ensemble.sbatch)


