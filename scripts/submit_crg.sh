#!/usr/bin/env bash
# Run on a CRG submission host after restoring renv and filtering.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
if [[ $# != 1 ]]; then echo 'Usage: scripts/submit_crg.sh results/NEW_RUN_NAME' >&2; exit 2; fi
command -v sbatch >/dev/null
mkdir -p logs
run_dir=$1
scripts/run_R.sh scripts/00_init_run.R "$run_dir"
job=$(sbatch --parsable scripts/betas.sbatch "$run_dir")
printf '%s\n' "$job" > "$run_dir/slurm_jobs.txt"
printf 'Submitted single full-table simulation and CSV export job %s.\n' "${job%%;*}"
printf 'Monitor with: squeue -j %s\n' "${job%%;*}"
