#!/usr/bin/env bash
#SBATCH --job-name=size_effect_tests
#SBATCH --output=/home/groups/earlew/arlenlex/P1-model-dev/size_effect/out_files/%x_%A.out
#SBATCH --error=/home/groups/earlew/arlenlex/P1-model-dev/size_effect/out_files/%x_%A.err
#SBATCH --time=01:00:00
#SBATCH -p serc
#SBATCH --ntasks=9
#SBATCH --nodes=1


module purge
module use /home/groups/earlew/arlenlex/.local/modulefiles
module load lammps_lexi chemistry


# 1) Run LAMMPS for each domain size and granularity
Ls=(100 200 300 400 500 600) 
num_atoms=(003496 013987 031473 055952 087426 125900) 
ncores=(1 2 2 4 4 9 9)  
seeds2=(8302 3094 3894 1092 4999)
dts=(0.1 0.05 0.025 0.01 0.01 0.01)
run_time=1800
num_dumps=500

for idx in "${!Ls[@]}"; do
    L=${Ls[$idx]}
    N=${num_atoms[$idx]}
    NP=${ncores[$idx]}
    for s in "${seeds2[@]}"; do
        mpirun -oversubscribe -n $NP lmp \
            -var L_pack "$L" \
            -var natoms  "$N" \
            -var atomfile_dir "/oak/stanford/groups/earlew/arlenlex/P1-model-dev/packing/size-effect/ensemble/L${L}/atomfiles_s${s}/" \
            -var out_dir "/scratch/groups/earlew/arlenlex/P1-model-dev/size_effect/data/ensemble/L${L}/" \
            -var outfile_name "N${N}Dump_seed${s}.lammps" \
            -var run_time "$run_time" \
            -var num_dumps "$num_dumps" \
            -var knb "1e9" \
            -var eps_crit "0.00005" \ 
            -in  /home/groups/earlew/arlenlex/P1-model-dev/simple_shear/in.shear \
            -log none
    done
done

# 2) Switch to ovito conda environment to postprocess
ml purge
source /home/groups/earlew/arlenlex/miniconda3-local/etc/profile.d/conda.sh
conda activate ovito
export PATH="$CONDA_PREFIX/bin:$PATH"

# Create numpy arrays of stresses using Ovito (xarray incompatible with ovito2 on Sherlock
python3 /home/groups/earlew/arlenlex/P1-model-dev/simple_shear/size_effect/analysis_code/compute_macroscopic_stress_ensemble.py \
    "/scratch/groups/earlew/arlenlex/P1-model-dev/simple_shear/size_effect/data/ensemble/" 

# 3) Turn numpy arrays into a dataset for visualization
ml purge
source /home/groups/earlew/arlenlex/miniconda3-local/etc/profile.d/conda.sh
conda activate nextsim
export PATH="$CONDA_PREFIX/bin:$PATH"

python3 /home/groups/earlew/arlenlex/P1-model-dev/simple_shear/size_effect/analysis_code/numpy_to_nc_of_stresses.py \
    "/scratch/groups/earlew/arlenlex/P1-model-dev/simple_shear/size_effect/data/ensemble/" 