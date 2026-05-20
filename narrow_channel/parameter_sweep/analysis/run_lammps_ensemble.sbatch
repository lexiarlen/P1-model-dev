#!/usr/bin/env bash
#SBATCH --job-name=run_lammps
#SBATCH --output=/scratch/groups/earlew/arlenlex/P1-model-dev/narrow_channel/out_files/lammps_jobs/%x_%A_%a.out
#SBATCH --error=/scratch/groups/earlew/arlenlex/P1-model-dev/narrow_channel/out_files/lammps_jobs/%x_%A_%a.err
#SBATCH --time=44:00:00            
#SBATCH -p serc
#SBATCH --array=1-1%50              # launcher will override with --array=1-${J}%${J}
#SBATCH --nodes=1
#SBATCH --ntasks=3                 # 4 MPI ranks per ensemble member
#SBATCH --cpus-per-task=1          # 1 CPU per MPI rank

set -euo pipefail

# ---------- required inputs via --export from launcher ----------
: "${PARAMS:?PARAMS (theta.tsv) path is required}"
: "${J:?J (ensemble size) is required}"
: "${output_dir:?output_dir is required}"
: "${packing_dir:?packing_dir is required}"

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

echo "RUNNING member=${MEMBER} ranks=${SLURM_NTASKS} kb=${kb} eps=${eps}"

# ---------- run LAMMPS on 3 MPI ranks ----------
mpirun -oversubscribe -np "${SLURM_NTASKS}" lmp \
    -in /home/groups/earlew/arlenlex/P1-model-dev/narrow_channel/lammps_scripts/in.strait \
    -var knb "${kb}" \
    -var eps_crit "${eps}" \
    -var member "${MEMBER}" \
    -var outdir "${output_dir}" \
    -var atomfile_dir "${packing_dir}" \
    -log none
