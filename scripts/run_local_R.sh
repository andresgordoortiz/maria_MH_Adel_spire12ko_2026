#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
# Reuse the installed R/packages read-only; matching host libraries preserve their ABI.
runtime=${CONTAINER_RUNTIME:-podman}
r_root=$(readlink -f "$HOME/.local/share/rig/r/4.5.2")
betas=$(find "$HOME/.cache/R/renv/library" -path '*/R-4.5/*/betAS' -print -quit)
[[ -x "$r_root/bin/Rscript" && -n "$betas" ]]
mounts=()
for path in /usr /lib /lib64 "$HOME/.local/share/rig" "$HOME/.cache/R/renv" "$PWD"; do
  mounts+=(-v "$path:$path:ro")
done
for path in results reports data/processed/spire_fastq metadata/spire_fastq audit; do
  mkdir -p "$path"
  mounts+=(-v "$PWD/$path:$PWD/$path:rw")
done
mounts+=(-v "$PWD/data/inputs/spire_fastq:$PWD/data/inputs/spire_fastq:ro")
exec "$runtime" run --rm --network=none "${mounts[@]}" -w "$PWD" \
  -e R_HOME="$r_root/lib/R" -e R_LIBS_USER="${betas%/betAS}" \
  -e FONTCONFIG_FILE="$HOME/.local/share/rig/fontconfig/fonts.conf" \
  -e RENV_CONFIG_AUTOLOADER_ENABLED=FALSE -e BETAS_DATASET="${BETAS_DATASET:-spire_fastq}" \
  -e MATT_OUTPUT_DIR="${MATT_OUTPUT_DIR:-data/processed/spire_fastq/matt}" \
  -e RMAPS_OUTPUT_DIR="${RMAPS_OUTPUT_DIR:-reports/spire_final/rmaps3_outputs}" \
  -e OPENBLAS_NUM_THREADS=1 -e OMP_NUM_THREADS=1 \
  --entrypoint "$r_root/bin/Rscript" \
  docker.io/andresgordoortiz/splicing_analysis_r_crg:v1.5 --vanilla "$@"
