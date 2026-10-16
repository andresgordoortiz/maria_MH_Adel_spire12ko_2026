#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
version=$(tr -d '[:space:]' < .R-version)
rig add "$version" --user --without-pak --without-fonts
exec rig run --user -r "$version" --rscript -f scripts/setup.R
