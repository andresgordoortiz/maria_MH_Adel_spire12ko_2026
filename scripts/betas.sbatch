#!/usr/bin/env bash
#SBATCH --job-name=spire-betAS
#SBATCH --output=logs/%x.%A_%a.out
#SBATCH --error=logs/%x.%A_%a.err
#SBATCH --qos=vlong
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the project root}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
if [[ $# != 1 ]]; then echo 'Expected run directory' >&2; exit 2; fi
exec scripts/run_R.sh scripts/02_simulate.R "$(( ${SLURM_ARRAY_TASK_ID:?Array index required} + ${SPIRE_CHUNK_OFFSET:-0} + 1 ))" "$1"
