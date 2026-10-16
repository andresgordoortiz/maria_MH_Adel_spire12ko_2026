#!/usr/bin/env bash
# Run this on a CRG submission host after restoring renv and filtering.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
if [[ $# != 1 ]]; then echo 'Usage: scripts/submit_crg.sh results/NEW_RUN_NAME' >&2; exit 2; fi
command -v sbatch >/dev/null
mkdir -p logs
run_dir=$1
scripts/run_R.sh scripts/00_init_run.R "$run_dir"
chunks=$(($(wc -l < data/processed/chunks.csv) - 1))
# Submit in groups to accommodate sites with MaxArraySize=1001; IDs stay global.
# Each group has its own dependency so collection waits for every chunk.
window=${ARRAY_WINDOW:-1000}
concurrent=${MAX_CONCURRENT:-20}
[[ "$window" =~ ^[1-9][0-9]*$ && "$concurrent" =~ ^[1-9][0-9]*$ ]] || exit 2
jobs=()
# Global task indexes cannot exceed MaxArraySize. Use an offset per array instead.
for ((offset=0; offset<chunks; offset+=window)); do
  count=$((chunks-offset)); ((count>window)) && count=$window
  job=$(sbatch --parsable --array="0-$((count-1))%${concurrent}" \
    --export="ALL,SPIRE_CHUNK_OFFSET=${offset}" scripts/betas.sbatch "$run_dir")
  jobs+=("${job%%;*}")
  printf '%s\n' "$job" >> "$run_dir/slurm_jobs.txt"
  printf 'Submitted simulation array %s (%s chunks).\n' "${job%%;*}" "$count"
done
dependency=$(IFS=:; echo "${jobs[*]}")
collection=$(sbatch --parsable --dependency="afterok:${dependency}" scripts/collect.sbatch "$run_dir")
printf '%s\n' "$collection" >> "$run_dir/slurm_jobs.txt"
printf 'Submitted collection job %s, after all simulations succeed.\n' "${collection%%;*}"
printf 'Monitor with: squeue -j %s\n' "$(IFS=,; echo "${jobs[*]},${collection%%;*}")"
