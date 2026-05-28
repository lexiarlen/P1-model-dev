#!/usr/bin/env bash
#SBATCH --job-name=kc_ensemble
#SBATCH --output=/scratch/groups/earlew/arlenlex/P1-model-dev/narrow_channel/kc_checks/out_files/lammps_jobs/%x_%A_%a.out
#SBATCH --error=/scratch/groups/earlew/arlenlex/P1-model-dev/narrow_channel/kc_checks/out_files/lammps_jobs/%x_%A_%a.err
#SBATCH --time=10:00:00            
#SBATCH -p serc
#SBATCH --array=0-4%5
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=1

set -euo pipefail

# ---------- environment ----------
ml purge
ml use /home/groups/earlew/arlenlex/.local/modulefiles
ml load lammps_lexi chemistry

export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

# ---------- fixed params ----------
kb=1e9
eps=0.00002

# ---------- kc sweep ----------
kc_vals=(1e5 1e6 1e7 1e8 1e9)

# ---------- dirs ----------
packing_dir="/scratch/groups/earlew/arlenlex/P1-model-dev/packing/narrow_channel/atomfiles_r400/"
base_dir="/scratch/groups/earlew/arlenlex/P1-model-dev/narrow_channel_smno/kc_ensemble_small"
state_dir="${base_dir}/state"
output_dir="${base_dir}/lammps_output"

mkdir -p "$state_dir"
mkdir -p "$output_dir"

# ---------- current array task ----------
idx="${SLURM_ARRAY_TASK_ID}"
kc="${kc_vals[$idx]}"
MEMBER=$((idx + 1))

# ---------- sanity checks ----------
: "${kc:?missing kc for idx=${idx}}"

echo "RUNNING idx=${idx} member=${MEMBER} ranks=${SLURM_NTASKS}"
echo "kb=${kb}"
echo "eps=${eps}"
echo "kc=${kc}"
echo "output_dir=${output_dir}"
echo "packing_dir=${packing_dir}"

# ---------- optional parameter table ----------
PARAMS="${state_dir}/theta.tsv"

if [[ "$idx" -eq 0 ]]; then
    echo -e "member\tkb\teps\tkc" > "$PARAMS"

    for i in "${!kc_vals[@]}"; do
        m=$((i + 1))
        echo -e "${m}\t${kb}\t${eps}\t${kc_vals[$i]}" >> "$PARAMS"
    done
fi

# ---------- run LAMMPS ----------
mpirun -oversubscribe -np "${SLURM_NTASKS}" lmp \
    -in /home/groups/earlew/arlenlex/P1-model-dev/narrow_channel/lammps_scripts/in.strait \
    -var knb "${kb}" \
    -var eps_crit "${eps}" \
    -var knc "${kc}" \
    -var member "${MEMBER}" \
    -var outdir "${output_dir}/" \
    -var atomfile_dir "${packing_dir}" \
    -log none