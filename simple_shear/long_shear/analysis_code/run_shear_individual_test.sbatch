#!/usr/bin/env bash
#SBATCH --job-name=floe_formation_test
#SBATCH --output=/home/groups/earlew/arlenlex/P1-model-dev/shear/out_files/test_sh_%j.out
#SBATCH --error=/home/groups/earlew/arlenlex/P1-model-dev/shear/out_files/test_sh_%j.err
#SBATCH --time=44:00:00
#SBATCH -p serc
#SBATCH -c 4

module purge
module use /home/groups/earlew/arlenlex/.local/modulefiles
module load lammps_lexi chemistry

EXEC="lmp"

SCRIPT="/home/groups/earlew/arlenlex/P1-model-dev/shear/lammps_scripts/in.shear"

mpirun -oversubscribe -n 4 $EXEC \
        -var L_pack "400" \
        -var natoms  "55952" \
        -var dt "0.025" \
        -log none \
        -in  $SCRIPT