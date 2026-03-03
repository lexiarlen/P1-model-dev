#!/usr/bin/env bash
#SBATCH --job-name=run_model
#SBATCH --output=/scratch/groups/earlew/arlenlex/P1-model-dev/nares/out_files/lammps_jobs/%x_%A_%a.out
#SBATCH --error=/scratch/groups/earlew/arlenlex/P1-model-dev/nares/out_files/lammps_jobs/%x_%A_%a.err
#SBATCH --time=15:00:00            # modify for big sims!
#SBATCH -p serc
#SBATCH --array=1-1%1              # launcher will override with --array=1-${J}%${J}
#SBATCH --nodes=1
#SBATCH --ntasks=3                 # 4 MPI ranks per ensemble member
#SBATCH --cpus-per-task=1          # 1 CPU per MPI rank

set -euo pipefail

# ---------- required inputs via --export from launcher ----------
: "${PARAMS:?PARAMS (theta.tsv) path is required}"
: "${J:?J (ensemble size) is required}"
: "${max_wind_speed:?max_wind_speed is required}"
: "${rho_i:?rho_i is required}"
: "${rho_a:?rho_a is required}"
: "${C_a:?C_a is required}"
: "${T_max:?T_max is required}"
: "${rho_o:?rho_o is required}"
: "${damping_ratio:?damping_ratio is required}"
: "${parent_dir_name:?parent_dir_name is required}"
: "${power:?power is required}"
: "${packing_dir:?packing_dir is required}"
: "${natoms:?natoms is required}"
: "${radius:?radius is required}"
: "${domain_scaling:?domain_scaling is required}"

# ---------- environment ----------
ml purge
ml use /home/groups/earlew/arlenlex/.local/modulefiles
ml load lammps_lexi chemistry

# pure MPI: no OpenMP threads
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

# ---------- current ensemble member (1-based) ----------
MEMBER="${SLURM_ARRAY_TASK_ID}"

# ---------- read the MEMBER-th nonblank, noncomment line as "kb eps" ----------
read -r kb eps < <(awk 'NF && $1 !~ /^#/' "${PARAMS}" | sed -n "${MEMBER}p")

# ---------- sanity checks ----------
: "${kb:?missing kb for MEMBER=${MEMBER}}"
: "${eps:?missing eps for MEMBER=${MEMBER}}"

echo "MEMBER=${MEMBER} ranks=${SLURM_NTASKS} kb=${kb} eps=${eps} max_wind_speed=${max_wind_speed} rho_i=${rho_i} rho_a=${rho_a} power=${power}"

# ---------- run LAMMPS on 4 MPI ranks ----------
mpirun -oversubscribe -np "${SLURM_NTASKS}" lmp \
    -in /home/groups/earlew/arlenlex/P1-model-dev/nares/lammps_scripts/in.strait \
    -var knb "${kb}" \
    -var eps_crit "${eps}" \
    -var member "${MEMBER}" \
    -var max_wind_speed "${max_wind_speed}" \
    -var rho_i "${rho_i}" \
    -var rho_a "${rho_a}" \
    -var rho_o "${rho_o}" \
    -var drag_a "${C_a}" \
    -var T_max "${T_max}" \
    -var damping_ratio "${damping_ratio}" \
    -var parent_dir_name "${parent_dir_name}" \
    -var atomfile_dir "${packing_dir}" \
    -var radius "${radius}" \
    -var natoms "${natoms}" \
    -var power "${power}" \
    -var domain_scaling "${domain_scaling}" \
    -log none
