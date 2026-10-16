#!/usr/bin/env bash
# Run under srun: install/restore the environment, then prepare filtered inputs.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
export PATH="$HOME/.local/bin:$PATH"
if ! command -v rig >/dev/null 2>&1; then
  echo 'rig must be installed and available on PATH before running setup.' >&2
  exit 1
fi
scripts/setup_env.sh
scripts/run_R.sh scripts/01_filter.R
echo 'Preparation complete. Submit from the login node with:'
echo 'scripts/submit_crg.sh results/spire_N10_nsim1000_run01'
