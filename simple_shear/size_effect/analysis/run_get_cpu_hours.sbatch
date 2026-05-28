#!/usr/bin/env bash
#SBATCH --job-name=parse_size_effect
#SBATCH --output=/home/groups/earlew/arlenlex/P1-model-dev/simple_shear/size_effect/out_files_smno/parse_size_effect_%j.out
#SBATCH --error=/home/groups/earlew/arlenlex/P1-model-dev/simple_shear/size_effect/out_files_smno/parse_size_effect_%j.err
#SBATCH --time=00:5:00
#SBATCH -p serc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

set -euo pipefail

# ---------- environment ----------
ml purge
source /home/groups/earlew/arlenlex/miniconda3-local/etc/profile.d/conda.sh
conda activate ovito

# ---------- paths ----------
script_path="/home/groups/earlew/arlenlex/P1-model-dev/simple_shear/size_effect/get_cpuhours.py"

base_dir="/home/groups/earlew/arlenlex/P1-model-dev/simple_shear/size_effect/out_files_smno"

out_dir="/scratch/groups/earlew/arlenlex/P1-model-dev/simple_shear/size_effect/ensemble10_smno"
out_path="${out_dir}/runs_dict.json"

# ---------- files to parse ----------
files=(
    "size_effect_tests_22983281.out"
)

# ---------- run parser ----------
python "$script_path" \
    --base-dir "$base_dir" \
    --files "${files[@]}" \
    --out-path "$out_path"