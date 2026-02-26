#!/usr/bin/env bash
#SBATCH --job-name=kb2elastic_ensemble
#SBATCH --output=/home/groups/earlew/arlenlex/P1-model-dev/simple_shear/kb2elastic/out_files/%x_%A.out
#SBATCH --error=/home/groups/earlew/arlenlex/P1-model-dev/simple_shear/kb2elastic/out_files/%x_%A.err
#SBATCH --time=00:20:00
#SBATCH -p serc
#SBATCH --ntasks=4
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
eps=(0.0005 0.0000416667 0.0000217391 0.0000147059 0.0000111111 0.0000089286 0.0000074627 0.0000064103 0.000005618 0.000005)
dts=(0.1 0.025 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.0075)
base_out="/scratch/groups/earlew/arlenlex/P1-model-dev/simple_shear/kb2elastic/L${L}"
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

#     # filename pattern MUST match N*.lammps for compute_elastic_moduli.py
#     outfile="N${num_atoms}_s${s}_kb${kb}_eps${ep}.lammps"

#     echo "[RUN] kb=${kb} eps=${ep} dt=${dt} -> ${case_dir}/${outfile}"

#     mpirun -n 4 lmp \
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

# --- postprocess part 1 (ovito env) ---
ml purge
source /home/groups/earlew/arlenlex/miniconda3-local/etc/profile.d/conda.sh
conda activate ovito
export PATH="$CONDA_PREFIX/bin:$PATH"

python3 /home/groups/earlew/arlenlex/P1-model-dev/simple_shear/kb2elastic/analysis/compute_elastic_moduli_qc.py \
    --base-dir "${base_out}" --qc-plots

# --- postprocess part 2 (pandas env) ---
ml purge
source /home/groups/earlew/arlenlex/miniconda3-local/etc/profile.d/conda.sh
conda activate postprocess_npy
export PATH="$CONDA_PREFIX/bin:$PATH"

python3 /home/groups/earlew/arlenlex/P1-model-dev/simple_shear/kb2elastic/analysis/numpy_to_pd.py \
    --base-dir "${base_out}"
