#!/usr/bin/env bash
#SBATCH --job-name=spire-betAS
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --qos=vlong
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the project root}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
if [[ $# != 1 ]]; then echo 'Expected run directory' >&2; exit 2; fi
scripts/run_R.sh scripts/02_simulate.R "$1"
exec scripts/run_R.sh scripts/03_collect.R "$1"
