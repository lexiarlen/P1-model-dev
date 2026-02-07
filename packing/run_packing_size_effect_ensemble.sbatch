#!/usr/bin/env bash
#SBATCH --job-name=making_packings_for_size_effect_ensemble_runs
#SBATCH --output=/home/groups/earlew/arlenlex/P1-model-dev/packing/out_files/%x_%A.out
#SBATCH --error=/home/groups/earlew/arlenlex/P1-model-deve/packing/out_files/%x_%A.err
#SBATCH --time=12:00:00
#SBATCH -p serc
#SBATCH -c 2

module purge
module use /home/groups/earlew/arlenlex/.local/modulefiles
module load lammps_lexi chemistry

EXEC="lmp"

# 1) create 6 differently seeded packings for each domain size = 100e3/L
Ls=(100 200 300 400 500 600)
seeds1=(9834 1825 2893 2283 1923)               
seeds2=(8302 3094 3894 1092 4999)
for L in "${Ls[@]}"; do
    for idx in "${!seeds1[@]}"; do
        s1=${seeds1[$idx]}
        s2=${seeds2[$idx]}
        mpirun -oversubscribe -np 2 \
        $EXEC -var H "${L}" -var W "${L}" \
        -var seed1 "${s1}" -var seed2 "${s2}" \
        -var outdir "/scratch/groups/earlew/arlenlex/P1-model-dev/packing/size-effect/ensemble/L${L}/" \
        -var outfile_name "datafile_s${s2}.lammps" \
        -in /home/groups/earlew/arlenlex/P1-model-dev/packing/in.pack \
        -log none
    done
done

# 2) generate atomfiles with the python script
ml purge
source /home/groups/earlew/arlenlex/miniconda3-local/etc/profile.d/conda.sh
conda activate nextsim
export PATH="$CONDA_PREFIX/bin:$PATH"

true_L=100000
for L in "${Ls[@]}"; do
    for s2 in "${seeds2[@]}"; do
        radius_scaling=$(bc <<< "$true_L / $L")

        outdir="/scratch/groups/earlew/arlenlex/P1-model-dev/packing/size-effect/ensemble/L${L}/atomfiles_s${s2}"
        mkdir -p "$outdir"

        python3 create_atomfiles.py \
        "/scratch/groups/earlew/arlenlex/P1-model-dev/packing/size-effect/ensemble/L${L}/datafile_s${s2}.lammps" \
        "$radius_scaling" \
        "$outdir"
    done
done
