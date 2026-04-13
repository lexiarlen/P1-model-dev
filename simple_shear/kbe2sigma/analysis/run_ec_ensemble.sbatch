#!/usr/bin/env bash
#SBATCH --job-name=kb2elastic_ensemble
#SBATCH --output=/home/groups/earlew/arlenlex/P1-model-dev/simple_shear/kbe2sigma/out_files/%x_%A.out
#SBATCH --error=/home/groups/earlew/arlenlex/P1-model-dev/simple_shear/kbe2sigma/out_files/%x_%A.err
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

kb=1e9
dt=0.025
eps=(5.00e-06 6.00e-05 1.15e-04 1.70e-04 2.25e-04 2.80e-04 3.35e-04 3.90e-04 4.45e-04 5.00e-04) 
base_out="/scratch/groups/earlew/arlenlex/P1-model-dev/simple_shear/kbe2sigma/L${L}"
atomfile_dir="/scratch/groups/earlew/arlenlex/P1-model-dev/packing/size-effect/ensemble10/L${L}/atomfiles_s${s}/"
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
        -in /home/groups/earlew/arlenlex/P1-model-dev/simple_shear/in.shear \
        -log none
done

# 1) create numpy files for each ensemble member
ml purge
source /home/groups/earlew/arlenlex/miniconda3-local/etc/profile.d/conda.sh
conda activate ovito
export PATH="$CONDA_PREFIX/bin:$PATH"

python3 /home/groups/earlew/arlenlex/P1-model-dev/simple_shear/kbe2sigma/analysis/compute_max_strengths.py --base-dir "${base_out}" --qc-plots

# 2) convert numpy files to parquet
ml purge
source /home/groups/earlew/arlenlex/miniconda3-local/etc/profile.d/conda.sh
conda activate postprocess_npy
export PATH="$CONDA_PREFIX/bin:$PATH"

python3 /home/groups/earlew/arlenlex/P1-model-dev/simple_shear/kbe2sigma/analysis/numpy_to_pd.py --base-dir "${base_out}"
