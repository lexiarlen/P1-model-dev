#!/usr/bin/env bash
#SBATCH --job-name=pack_nares
#SBATCH --output=/home/groups/earlew/arlenlex/P1-model-dev/packing/out_files/%x_%A.out
#SBATCH --error=/home/groups/earlew/arlenlex/P1-model-dev/packing/out_files/%x_%A.err
#SBATCH --time=00:10:00
#SBATCH -p serc
#SBATCH -c 1

# 1) create packing for narrow channel simulations with radius = 400 m => granularity <= 0.033
module purge
module use /home/groups/earlew/arlenlex/.local/modulefiles
module load lammps_lexi chemistry

EXEC="lmp"

mpirun -oversubscribe -np 1 \
      $EXEC -var H "600" -var W "250" \
      -var outdir "/scratch/groups/earlew/arlenlex/P1-model-dev/packing/narrow_channel/" \
      -var outfile_name "datafile_r400.lammps" \
      -in /home/groups/earlew/arlenlex/P1-model-dev/packing/lammps_scripts/in.pack -log none

# 2) Switch to your conda environment and generate atomfiles
ml purge
source /home/groups/earlew/arlenlex/miniconda3-local/etc/profile.d/conda.sh
conda activate nextsim
export PATH="$CONDA_PREFIX/bin:$PATH"

radius_scaling=400

outdir="/scratch/groups/earlew/arlenlex/P1-model-dev/packing/narrow_channel/atomfiles_r400/"
mkdir -p "$outdir"

python3 create_atomfiles.py \
    "/scratch/groups/earlew/arlenlex/P1-model-dev/packing/narrow_channel/datafile_r400.lammps" \
    "$radius_scaling" \
    "$outdir"
