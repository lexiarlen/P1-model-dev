#!/usr/bin/env bash
#SBATCH --job-name=size_effect_tests
#SBATCH --output=/home/groups/earlew/arlenlex/P1-model-dev/size_effect/out_files/%x_%A.out
#SBATCH --error=/home/groups/earlew/arlenlex/P1-model-dev/size_effect/out_files/%x_%A.err
#SBATCH --time=00:20:00
#SBATCH -p serc
#SBATCH --ntasks=9
#SBATCH --nodes=1


module purge
module use /home/groups/earlew/arlenlex/.local/modulefiles
module load lammps_lexi chemistry

EXEC="lmp"

# 1) Run LAMMPS for each domain size (serially)
Ls=(600) #(50 100 200 300 400 500 600) 
num_atoms=(125900) #(000873 003496 013987 031473 055952 087426 105785)
ncores=(9) #(1 2 2 4 4 9 9)  
seeds1=(9834 0825 2893 2283 1923)
seeds2=(8302 3094 3894 1092 4999)

for idx in "${!Ls[@]}"; do
    L=${Ls[$idx]}
    N=${num_atoms[$idx]}
    NP=${ncores[$idx]}
    for s in "${seeds2[@]}"; do
        mpirun -oversubscribe -n $NP $EXEC \
            -var L_pack "$L" \
            -var natoms  "$N" \
            -var atomfile_dir "/oak/stanford/groups/earlew/arlenlex/P1-model-dev/packing/size-effect/ensemble/L${L}/atomfiles_s${s}/" \
            -var out_dir "/scratch/groups/earlew/arlenlex/P1-model-dev/size_effect/data/big_natoms2/L${L}/" \
            -var outfile_name "N${N}Dump_seed${s}.lammps" \
            -in  /home/groups/earlew/arlenlex/P1-model-dev/size_effect/lammps_scripts/in.shear \
            -log none
    done
done

# 2) Switch to conda/python environment and post-process
ml purge
source /home/groups/earlew/arlenlex/miniconda3-local/etc/profile.d/conda.sh
conda activate ovito
export PATH="$CONDA_PREFIX/bin:$PATH"

# Create numpy arrays of stresses using Ovito (xarray incompatible)
python3 /home/groups/earlew/arlenlex/P1-model-dev/size_effect/analysis_code/compute_macroscopic_stress_ensemble.py \
    "/scratch/groups/earlew/arlenlex/P1-model-dev/size_effect/data/ensemble_final/"
