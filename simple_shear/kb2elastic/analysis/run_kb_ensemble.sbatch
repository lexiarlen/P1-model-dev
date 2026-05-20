#!/usr/bin/env bash
#SBATCH --job-name=kb2elastic_ensemble
#SBATCH --output=/home/groups/earlew/arlenlex/P1-model-dev/simple_shear/kb2elastic_smno/out_files/%x_%A.out
#SBATCH --error=/home/groups/earlew/arlenlex/P1-model-dev/simple_shear/kb2elastic_smno/out_files/%x_%A.err
#SBATCH --time=04:00:00
#SBATCH -p serc
#SBATCH --ntasks=2
#SBATCH --nodes=1

set -euo pipefail

module purge
module use /home/groups/earlew/arlenlex/.local/modulefiles
module load lammps_lexi chemistry

# --- ensemble config ---
L=300
num_atoms=031475
s=8302

kbs=(1e8 1.2e9 2.3e9 3.4e9 4.5e9 5.6e9 6.7e9 7.8e9 8.9e9 1.0e10)
eps=(0.01 0.00083333 0.00043478 0.00029412 0.00022222 0.00017857 0.00014925 0.00012821 0.00011236 0.0001) # set high strength b/c want pre-failure physics
dts=(0.1 0.025 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.0075)
base_out="/scratch/groups/earlew/arlenlex/P1-model-dev/simple_shear/kb2elastic_smno/L${L}"
atomfile_dir="/scratch/groups/earlew/arlenlex/P1-model-dev/packing/size-effect/ensemble10/L${L}/atomfiles_s${s}/"
run_time="1800"
num_dumps="500"

# mkdir -p "${base_out}"

# for idx in "${!kbs[@]}"; do
#     kb="${kbs[$idx]}"
#     ep="${eps[$idx]}"
#     dt="${dts[$idx]}"

#     case_dir="${base_out}/kb_${kb}_eps_${ep}"
#     mkdir -p "${case_dir}"

#     outfile="N${num_atoms}_s${s}_kb${kb}_eps${ep}.lammps"

#     echo "[RUN] kb=${kb} eps=${ep} dt=${dt} -> ${case_dir}/${outfile}"

#     mpirun -n 2 lmp \
#         -var L_pack "${L}" \
#         -var natoms  "${num_atoms}" \
#         -var knb "${kb}" \
#         -var eps_crit "${ep}" \
#         -var atomfile_dir "${atomfile_dir}" \
#         -var out_dir "${case_dir}/" \
#         -var outfile_name "${outfile}" \
#         -var run_time "${run_time}" \
#         -var num_dumps "${num_dumps}" \
#         -var dt "${dt}" \
#         -in /home/groups/earlew/arlenlex/P1-model-dev/simple_shear/in.shear \
#         -log none
# done

# 1) create numpy files for each ensemble member
ml purge
source /home/groups/earlew/arlenlex/miniconda3-local/etc/profile.d/conda.sh
conda activate ovito
export PATH="$CONDA_PREFIX/bin:$PATH"

python3 /home/groups/earlew/arlenlex/P1-model-dev/simple_shear/kb2elastic/analysis/compute_elastic_moduli.py --base-dir "${base_out}" 

# 2) convert numpy files to parquet
ml purge
source /home/groups/earlew/arlenlex/miniconda3-local/etc/profile.d/conda.sh
conda activate postprocess_npy
export PATH="$CONDA_PREFIX/bin:$PATH"

python3 /home/groups/earlew/arlenlex/P1-model-dev/simple_shear/kb2elastic/analysis/numpy_to_pd.py --base-dir "${base_out}"
