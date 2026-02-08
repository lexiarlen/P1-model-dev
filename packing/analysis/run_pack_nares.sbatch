#!/usr/bin/env bash
#SBATCH --job-name=pack_nares
#SBATCH --output=/home/groups/earlew/arlenlex/P1-model-dev/packing/out_files/%x_%A.out
#SBATCH --error=/home/groups/earlew/arlenlex/P1-model-deve/packing/out_files/%x_%A.err
#SBATCH --time=00:20:00
#SBATCH -p serc
#SBATCH -c 1

# 1) create packing for nares simulations with radius = 500 m => granularity < 0.05
module purge
module use /home/groups/earlew/arlenlex/.local/modulefiles
module load lammps_lexi chemistry

EXEC="lmp"

mpirun -oversubscribe -np 1 \
      $EXEC -var H "480" -var W "280" \
      -var outdir "/scratch/groups/earlew/arlenlex/P1-model-dev/packing/nares/" \
      -var outfile_name "atomfiles_r500.lammps" \
      -in /home/groups/earlew/arlenlex/P1-model-dev/packing/in.pack -log none

# 2) Switch to your conda environment and generate atomfiles
ml purge
source /home/groups/earlew/arlenlex/miniconda3-local/etc/profile.d/conda.sh
conda activate nextsim
export PATH="$CONDA_PREFIX/bin:$PATH"

radius_scaling=500

outdir="/scratch/groups/earlew/arlenlex/P1-model-dev/packing/nares/atomfiles_r500"
mkdir -p "$outdir"

python3 create_atomfiles.py \
    "/scratch/groups/earlew/arlenlex/P1-model-dev/packing/nares/datafile_r500.lammps" \
    "$radius_scaling" \
    "$outdir"
