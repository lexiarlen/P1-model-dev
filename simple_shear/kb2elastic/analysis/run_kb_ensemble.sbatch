#!/usr/bin/env bash
#SBATCH --job-name=size_effect_tests_ensemble_plastic
#SBATCH --output=/home/groups/earlew/arlenlex/P1-EKI-nares/size_effect_tests/out_files/%x_%A.out
#SBATCH --error=/home/groups/earlew/arlenlex/P1-EKI-nares/size_effect_tests/out_files/%x_%A.err
#SBATCH --time=24:00:00
#SBATCH -p serc
#SBATCH --ntasks=2
#SBATCH --nodes=1


module purge
module use /home/groups/earlew/arlenlex/.local/modulefiles
module load lammps_lexi chemistry

# 1) Run LAMMPS for each domain size (serially); small test run
L=200
num_atoms=13496
s=9834
kbs=(1e8 1.2e9 2.3e9 3.4e9 4.5e9 5.6e9 6.7e9 7.8e9 8.9e9 1.0e10) 
eps=(0.0005 0.0000416667 0.0000217391 0.0000147059 0.0000111111 0.0000089286 0.0000074627 0.0000064103 0.000005618 0.000005) # corresponds to sigma = 50,000 Pa
dts=(0.1 0.025 0.025 0.025 0.01 0.01 0.01 0.01 0.01 0.01)
s=8302

for idx in "${!kbs[@]}"; do
    kb=${kbs[$idx]}
    ep=${eps[$idx]}
    dt=${dts[$idx]}
    mpirun -oversubscribe -n 2 lmp \
        -var L_pack "$L" \
        -var natoms  "$num_atoms" \
        -var knb "$kb" \
        -var eps_crit "$ep" \
        -var atomfile_dir "/oak/stanford/groups/earlew/arlenlex/P1-model-dev/packing/size-effect/ensemble/L${L}/atomfiles_s${s}/" \
        -var out_dir "/scratch/groups/earlew/arlenlex/P1-model-dev/simple_shear/kb2elastic/L${L}/" \
        -var outfile_name "N${N}Dump_seed${s}.lammps" \
        -var run_time "1800" \
        -var num_dumps "500" \
        -var dt "$dt" \
        -in  /home/groups/earlew/arlenlex/P1-model-dev/simple_shear/in.shear \
        -log none 
done

# 2) Switch to conda environment w/ ovito and postprocess part 1
ml purge
source /home/groups/earlew/arlenlex/miniconda3-local/etc/profile.d/conda.sh
conda activate ovito

export PATH="$CONDA_PREFIX/bin:$PATH"

python3 /home/groups/earlew/arlenlex/P1-model-dev/simple_shear/kb2elastic/analysis/compute_elastic_moduli.py \
    "/scratch/groups/earlew/arlenlex/P1-model-dev/simple_shear/kb2elastic/L${L}/"

# 3) Switch to conda environment w/ xarray and postprocess part 2
ml purge
source /home/groups/earlew/arlenlex/miniconda3-local/etc/profile.d/conda.sh
conda activate nextsim

export PATH="$CONDA_PREFIX/bin:$PATH"

python3 /home/groups/earlew/arlenlex/P1-model-dev/simple_shear/kb2elastic/analysis/nnumpy_to_kb2e.py \
    "/scratch/groups/earlew/arlenlex/P1-model-dev/simple_shear/kb2elastic/L${L}/"