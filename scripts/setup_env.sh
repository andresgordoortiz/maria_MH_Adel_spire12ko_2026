#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.."

export PATH="$HOME/.local/bin:$PATH"

# Load cluster dependency required by igraph
module load GLPK/5.0-GCCcore-14.3.0

if ! command -v rig >/dev/null 2>&1; then
  bash scripts/install_rig.sh
fi

version=$(tr -d '[:space:]' < .R-version)

rig add "$version" --user --without-pak --without-fonts

exec rig run --user -r "$version" \
  --rscript -f scripts/setup.R