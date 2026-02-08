#!/usr/bin/env bash
#SBATCH --job-name=floe_formation_test
#SBATCH --output=/home/groups/earlew/arlenlex/P1-model-dev/simple_shear/long_shear/out_files/test_sh_%j.out
#SBATCH --error=/home/groups/earlew/arlenlex/P1-model-dev/simple_shear/long_shear/out_files/test_sh_%j.err
#SBATCH --time=24:00:00
#SBATCH -p serc
#SBATCH -c 2

module purge
module use /home/groups/earlew/arlenlex/.local/modulefiles
module load lammps_lexi chemistry

s=8302

mpirun -oversubscribe -n 2 lmp \
        -var L_pack "300" \
        -var natoms  "31473" \
        -var atomfile_dir "/oak/stanford/groups/earlew/arlenlex/P1-model-dev/packing/size-effect/ensemble/L${L}/atomfiles_s${s}/" \
        -var out_dir "/scratch/groups/earlew/arlenlex/P1-model-dev/size_effect/data/long_shear/" \
        -var outfile_name "N${N}Dump_seed${s}_1000hours.lammps" \
        -var run_time "3600000" \
        -var num_dumps "1000" \
        -var dt "0.025" \
        -var knb "1e9" \
        -var eps_crit "0.00005" \ 
        -in  /home/groups/earlew/arlenlex/P1-model-dev/simple_shear/in.shear \
        -log none