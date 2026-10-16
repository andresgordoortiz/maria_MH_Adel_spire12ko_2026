#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
if [[ $# != 2 ]]; then
  echo 'Usage: bash analysis/run_matt.sh ANALYSIS_DIR NEW_MATT_DIR' >&2
  exit 1
fi
analysis_dir=$1
matt_dir=$2
[[ -f "$analysis_dir/COMPLETE" && ! -e "$matt_dir" ]] || { echo 'Need a completed analysis and a new Matt output directory.' >&2; exit 1; }
for path in "$analysis_dir" "$matt_dir"; do
  [[ "$path" =~ ^[A-Za-z0-9_./-]+$ && "$path" != /* && "/$path/" != *'/../'* ]] || {
    echo 'Use repository-relative paths with letters, digits, /, _, - or . only; upstream Matt interpolates filenames.' >&2; exit 1;
  }
done
image='docker.io/andresgordoortiz/matt-container@sha256:6b729b82e312f06153f8faeb643402df4e1729485fba90359a6cf3928c24c20e'
runtime=${CONTAINER_RUNTIME:-podman}
command -v "$runtime" >/dev/null
[[ -f data/reference/gencode.vM10.annotation.gtf && -f data/reference/GRCm38.primary_assembly.genome.fa ]]
mkdir -p "$matt_dir"
python3 - "$analysis_dir" "$matt_dir" <<'PY'
import csv, pathlib, sys
from analysis.prepare_references import file_hash
source, dest = map(pathlib.Path, sys.argv[1:])
with (dest/'input_checksums.csv').open('w',newline='') as f:
    w=csv.writer(f);w.writerow(['file','sha256'])
    for kind in ('EX','IR'):
        name=f'spire_{kind}_matt.tab'
        w.writerow([name,file_hash(source/'external_inputs'/name)])
PY
extra=()
[[ "$runtime" != podman ]] || extra+=(--userns=keep-id)
for kind in EX IR; do
  command=cmpr_exons
  [[ "$kind" != IR ]] || command=cmpr_introns
  "$runtime" run --rm "${extra[@]}" -v "$PWD:/data" -w "/data/$matt_dir" --entrypoint /bin/bash "$image" -c '
    set -euo pipefail
    export RENV_CONFIG_AUTOLOADER_ENABLED=FALSE
    matt "$1" "/data/$2/external_inputs/spire_$3_matt.tab" START END SCAFFOLD STRAND GENEID_ENSEMBL \
      /data/data/reference/gencode.vM10.annotation.gtf /data/data/reference/GRCm38.primary_assembly.genome.fa \
      Mmus 150 "DATASET[up,down,ndiff]" "$3" -notrbts
  ' _ "$command" "$analysis_dir" "$kind" "$matt_dir" > "$matt_dir/$kind.log" 2>&1
  feature_suffix=efeatures
  [[ "$kind" != IR ]] || feature_suffix=ifeatures
  [[ -s "$matt_dir/$kind/spire_${kind}_matt_with_${feature_suffix}.tab" ]]
done
python3 analysis/validate_matt.py "$analysis_dir" "$matt_dir"
printf '%s\n' "$image" > "$matt_dir/COMPLETE"
echo "Matt complete: $matt_dir"
