#!/usr/bin/env bash
#SBATCH --job-name=size_effect_tests
#SBATCH --output=/home/groups/earlew/arlenlex/P1-model-dev/simple_shear/size_effect/out_files_smno/%x_%A.out
#SBATCH --error=/home/groups/earlew/arlenlex/P1-model-dev/simple_shear/size_effect/out_files_smno/%x_%A.err
#SBATCH --time=72:00:00
#SBATCH -p serc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=9

module purge
module use /home/groups/earlew/arlenlex/.local/modulefiles
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
            -var atomfile_dir "/scratch/groups/earlew/arlenlex/P1-model-dev/packing/size-effect/ensemble10/L${L}/atomfiles_s${s}/" \
            -var out_dir "/scratch/groups/earlew/arlenlex/P1-model-dev/simple_shear/size_effect/ensemble10_smno/L${L}/" \
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

# Create numpy arrays of stresses using Ovito (pandas incompatible with ovito on Sherlock)
python3 /home/groups/earlew/arlenlex/P1-model-dev/simple_shear/size_effect/analysis/compute_macroscopic_stress.py \
    "/scratch/groups/earlew/arlenlex/P1-model-dev/simple_shear/size_effect/ensemble10_smno/" 

# 3) Turn numpy arrays into parquet for visualization
ml purge
source /home/groups/earlew/arlenlex/miniconda3-local/etc/profile.d/conda.sh
conda activate postprocess_npy
export PATH="$CONDA_PREFIX/bin:$PATH"

python3 /home/groups/earlew/arlenlex/P1-model-dev/simple_shear/size_effect/analysis/numpy_to_nc_of_stresses.py \
    "/scratch/groups/earlew/arlenlex/P1-model-dev/simple_shear/size_effect/ensemble10_smno/" 
