#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
export PATH="$HOME/.local/bin:$PATH"
# An explicitly supplied Rscript can be a CRG module or container wrapper.
# Each R script checks the exact version and pinned betAS revision.
if [[ -n "${RSCRIPT:-}" ]]; then
  exec "$RSCRIPT" "$@"
fi
version=$(tr -d '[:space:]' < .R-version)
if command -v rig >/dev/null 2>&1; then
  # Resolve the rig-managed executable, then pass Rscript arguments verbatim.
  # rig 0.10's -f/--rscript path adds a literal --args to script arguments.
  r_bin=$(rig run --user -r "$version" --rscript -e 'cat(R.home("bin"))' -- --vanilla)
  exec "$r_bin/Rscript" "$@"
fi
exec Rscript "$@"
