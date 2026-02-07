#!/usr/bin/env bash
#SBATCH --job-name=size_effect_tests_ensemble_plastic
#SBATCH --output=/home/groups/earlew/arlenlex/P1-EKI-nares/size_effect_tests/out_files/%x_%A.out
#SBATCH --error=/home/groups/earlew/arlenlex/P1-EKI-nares/size_effect_tests/out_files/%x_%A.err
#SBATCH --time=66:00:00
#SBATCH -p serc
#SBATCH --ntasks=9
#SBATCH --nodes=1


module purge
module use /home/groups/earlew/arlenlex/.local/modulefiles
module load lammps_lexi chemistry

EXEC="lmp"

# 1) Run LAMMPS for each domain size (serially); small test run
L=50 #300
num_atoms=873 #031473
ncores=1 #4  
dt=0.1 #0.001 
s=9834
kbs=(1e8 5e8 1e9 2.5e9 5e9 7.5e9 1e10) # should probs latin hypercube?
eps=(0.000500 0.000100 0.000050 0.000020 0.000010 0.000006667 0.000005) # corresponds to sigma = 50,000 Pa
#seeds1=(9834 0825 2893 2283 1923)

for idx in "${!kbs[@]}"; do
    kb=${kbs[$idx]}
    ep=${eps[$idx]}
    mpirun -oversubscribe -n $NP $EXEC \
        -var L_pack "$L" \
        -var natoms  "$num_atoms" \
        -var dt "$dt" \
        -var knb "$kb" \
        -var eps_crit "$ep" \
        -var atomfile_dir "/oak/stanford/groups/earlew/arlenlex/P1-EKI-nares/packing/size-effect/ensemble/L${L}/atomfiles_s${s}/" \
        -var out_dir "/scratch/groups/earlew/arlenlex/P1-EKI-nares/kb2elastic/ensemble_kb/L${L}/" \
        -var outfile_name "N${N}Dump_seed${s}.lammps" \
        -in  /home/groups/earlew/arlenlex/P1-EKI-nares/size_effect_tests/in.shear
done

# 2) Switch to conda/python environment and post-process
ml purge
source /home/groups/earlew/arlenlex/miniconda3-local/etc/profile.d/conda.sh
conda activate ovito

# Ensure the env’s bin is at the front of PATH:
export PATH="$CONDA_PREFIX/bin:$PATH"

# make sure the script path is correct
python3 /home/groups/earlew/arlenlex/P1-EKI-nares/size_effect_tests/compute_macroscopic_stress_ensemble.py \
    "/scratch/groups/earlew/arlenlex/P1-EKI-nares/size_effect_tests/data/ensemble_anharmonic/"
