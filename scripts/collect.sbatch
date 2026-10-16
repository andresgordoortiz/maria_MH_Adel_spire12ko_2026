#!/usr/bin/env bash
#SBATCH --job-name=spire-collect
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --qos=short
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the project root}"
exec scripts/run_R.sh scripts/03_collect.R "$1"
