> Historical exploratory handoff, not used in the publication figures. The
> SPIRE counts and prepared files below predate the corrected FASTQ analysis.
> For new inputs use dataset `spire_fastq` and a new output directory; the
> container recipe remains unbuilt/unvalidated. Current publication motifs
> come from the completed rMAPS3 analysis.

# RNAmotifs 2.0: checked handoff for Spire and FMN2

Checked on 2026-10-17 against upstream commit
`94d46a76350ed8c7d78cbce111703e2ccbab4fe3` of
[matteocereda/RNAmotifs](https://github.com/matteocereda/RNAmotifs).

**EX/SE is runnable after conversion. IR/RI is experimental and has source-level
region/plotting problems; do not use it as a validated replacement. Alt5/Alt3
are unsupported.** This is motif-cluster discovery, not the same known-RBP
motif test or output format as rMAPS2.

## What is in this repository

`results/{spire,fmn2}_final/external_inputs/` contains custom rMAPS tables,
not rMATS files. These come from VASTDB/betAS on GRCm38/mm10, with 1-based
closed coordinates. `analysis/external.R` writes adjacent event-ID manifests.
The existing groups are higher inclusion/retention in KO (`up`), lower in KO
(`down`), and all other eligible tested events (`ndiff`).

| Study/type | up | down | ndiff |
| --- | ---: | ---: | ---: |
| Spire EX | 250 | 286 | 89,922 |
| Spire IR | 127 | 114 | 81,199 |
| FMN2 EX | 84 | 61 | 19,992 |
| FMN2 IR | 15 | 31 | 13,411 |

`prepare_rnamotifs.py` verifies exported coordinate-file hashes, validates
coordinates and paired manifests, and writes native semicolon-delimited,
headerless input. Numeric IDs map back to VASTDB EVENT IDs in
`event_manifest.csv`. All ndiff events are retained, with no sampling or
reclassification. The native control class therefore means this study's
unselected background, not demonstrated constitutive exons.

Native fields are `row_id;event_id;chr;strand;v5;v6;v7;v8;class[;RI]`.
The class values `2`, `-2`, `0` are category encodings, not delta-PSI values.
For EX, v5-v8 are left exon end, target start minus one, target end, right
exon start minus one. For IR they are left exon start minus one, left exon
end, right exon start minus one, right exon end. Coordinates remain in
ascending genomic order on both strands; RNAmotifs handles strand orientation.
Do not pass `--from-rmats` for these converted inputs.

The following outputs have already been prepared locally (ignored by Git):

```
data/processed/rnamotifs2_spire_inputs/spire_EX.txt
data/processed/rnamotifs2_spire_inputs/spire_IR.txt
data/processed/rnamotifs2_fmn2_inputs/fmn2_EX.txt
data/processed/rnamotifs2_fmn2_inputs/fmn2_IR.txt
```

To prepare them on another machine, from the repository root:

```bash
python3 tests/rnamotifs_checks.py
python3 analysis/prepare_rnamotifs.py reports/spire_final spire_fastq data/processed/rnamotifs2_spire_fastq_inputs
python3 analysis/prepare_rnamotifs.py reports/fmn2_final fmn2 data/processed/rnamotifs2_fmn2_inputs
```

Output directories must be new. Input hashes, row-ID manifests and counts are
saved alongside the input files; `COMPLETE` marks successful conversion.

## Container execution only

The earlier temporary host build, R test library and new renv staging/cache
entries were removed after the user requested environment isolation. See
`audit/rnamotifs2/environment_cleanup.json`. The active libraries and tracked
R environment files were checked. Future evaluations must use uv, pixi or
containers; this mixed Python/R/C++ stack uses a container.

[The Containerfile](rnamotifs.Containerfile) pins the RNAmotifs source revision
and installs its compilers and R dependencies inside the image. **This revised
container recipe has not been built or run.** The successful smoke tests below
are historical host-execution evidence, not container validation. R package
versions are not yet locked; validate and pin the built image before using it
for final analyses.

From the repository root, the following commands would build a local image and
create separate genome/output directories. They have not been executed:

```bash
export RNAMOTIFS_PROJECT="$PWD"
export RNAMOTIFS="$RNAMOTIFS_PROJECT/data/processed/rnamotifs2_container"
mkdir -p "$RNAMOTIFS/build-context" "$RNAMOTIFS/genomes/mm10" "$RNAMOTIFS/results"
podman build -f analysis/rnamotifs.Containerfile -t localhost/rnamotifs2:94d46a7 "$RNAMOTIFS/build-context"
```

## Reuse the existing mm10 genome

The checked upstream genome-download script additionally needs `twoBitToFa`,
which is not included. Avoid another genome download: split this repo's
existing GENCODE GRCm38 primary-assembly FASTA into the chromosome strings
RNAmotifs expects. This writes about 2.65 GB.

```bash
python3 - <<'PY'
import os
from pathlib import Path

project = Path(os.environ['RNAMOTIFS_PROJECT'])
destination = Path(os.environ['RNAMOTIFS']) / 'genomes/mm10'
destination.mkdir(parents=True, exist_ok=True)
chromosomes = {f'chr{i}' for i in range(1, 20)} | {'chrX', 'chrY'}
seen = set()
output = None
with (project / 'data/reference/GRCm38.primary_assembly.genome.fa').open() as fasta:
    for line in fasta:
        if line.startswith('>'):
            if output:
                output.close()
            chrom = line[1:].split()[0]
            output = (destination / f'{chrom}.string').open('x') if chrom in chromosomes else None
            if output:
                seen.add(chrom)
        elif output:
            output.write(line.strip().upper())
if output:
    output.close()
assert seen == chromosomes, f'Missing chromosomes: {chromosomes - seen}'
PY
```

## Run EX/SE

Keep the environment variables above in the shell. Use distinct names for
each study/event type/run: upstream writes intermediate data by analysis name
and can overwrite an existing same-day run with the same settings. Run the
studies sequentially initially.

```bash
podman run --rm --network=none \
  -v "$RNAMOTIFS_PROJECT/data/processed/rnamotifs2_spire_inputs:/inputs:ro" \
  -v "$RNAMOTIFS/genomes:/opt/RNAmotifs/genomes:ro" \
  -v "$RNAMOTIFS/results:/opt/RNAmotifs/results:rw" \
  localhost/rnamotifs2:94d46a7 /inputs/spire_EX.txt \
  -n spire_EX_v1 -g mm10 --event-type SE \
  --in-exon 50 --in-intron 250 --n 30 -e 30 \
  -b 1000 --p-empirical 0.01 -c 2 --top-n 20

podman run --rm --network=none \
  -v "$RNAMOTIFS_PROJECT/data/processed/rnamotifs2_fmn2_inputs:/inputs:ro" \
  -v "$RNAMOTIFS/genomes:/opt/RNAmotifs/genomes:ro" \
  -v "$RNAMOTIFS/results:/opt/RNAmotifs/results:rw" \
  localhost/rnamotifs2:94d46a7 /inputs/fmn2_EX.txt \
  -n fmn2_EX_v1 -g mm10 --event-type SE \
  --in-exon 50 --in-intron 250 --n 30 -e 30 \
  -b 1000 --p-empirical 0.01 -c 2 --top-n 20
```

`--in-exon 50 --in-intron 250` preserves the requested flanking extents,
but does not reproduce rMAPS2's algorithm or region tests. `--n 30` is the
motif-clustering window and `-e 30` is RNAmotifs' enrichment window.
`--top-n 20` limits plots only; bootstrap results retain all 512 candidate motifs.
These are explicit starting parameters, not parameters calibrated for oocytes.

Results are under `$RNAMOTIFS/results/NAME/DATE_mm10_.../`:
`bootstrap_1000.tsv`, `rnamotifs.log`, raw count/position data, and, if motifs
pass selection, `MRMs_*.csv`, positional CSVs and RNA-map PDFs. A successful
run with no selected motifs can produce no PDF; that is not an installation
failure. The current wrapper suppresses output from R, so for debugging run
the `Rscript src/R/selection.R ...` command directly with the paths and
parameters printed by the wrapper/traceback.

**Always specify a compatible empirical cutoff.** In this commit,
`pEmp = (1 + count)/(1 + bootstraps)`. The defaults `-b 10000` and
`--p-empirical 0.00005` cannot select a motif: the minimum attainable p-value
is about 0.0001. With 1,000 bootstraps the minimum is about 0.001; 0.01 is
attainable. Selection also applies an adaptive Fisher/BH criterion.

## Limits that matter here

- **RI source problems:** `rnamotifs_search.cpp::prepare_regions` uses fixed
  exon-skipping search windows for RI too. With 50/250 spans, these fail to
  cover some requested retained-intron sequence for 12,994 Spire and 2,511
  FMN2 IR events (including 60 and 8 differential events respectively), even
  after merging search windows from all events on each chromosome/strand.
  `selection.R` also reads SE-oriented positional windows, inconsistent with
  the RI positions emitted by `rnamotifs_core.cpp`. A separate
  `ir_splicing_map.R` exists but is not invoked by the normal CLI. Successful
  RI execution does not resolve these scientific-output issues.
- **RBP identification:** core RNAmotifs returns enriched short motif
  clusters, not the original panel of named RBPs. MaRs requires external
  eCLIP files and calibrated `Tables/` and `Rdata/` reference panels. Its CLI
  offers HepG2/K562 human-cell panels; selecting mm10 does not establish
  compatibility with mouse/oocyte binding data or provide a liftover.
  No suitable mouse/oocyte reference panel was verified here.
- **Existing figure pipeline:** outputs cannot be passed to
  `RMAPS_OUTPUT_DIR`/`analysis/render_rmaps3.R`. That renderer expects named-RBP
  rMAPS2 p-value tables with eight EX/four IR regions. RNAmotifs has different
  motif IDs, tests and region definitions. Keep its results separate.
- **Scale:** Spire EX has 90,458 rows, far above the 4,368-exon upstream
  benchmark. Full-data runtime, peak memory and scientific findings were
  not measured. Upstream preloads whole used chromosomes; the tiny subset
  memory figures below cannot be extrapolated to all 21 chromosomes. Use
  the workstation cautiously or an HPC allocation for full runs.

## Verification performed

1. Independently cross-checked every EX/IR export in both studies against
   VASTDB target/flanking coordinates, betAS direction thresholds and mm10
   chromosome lengths. No mismatch, invalid ordering or duplicate complete
   coordinate rows was found. Checked start-minus-one conversions and both
   strands with assertions; modified hashes and existing output paths fail.
2. Built the pinned upstream C++ code and loaded the required R packages.
3. Ran 512-motif/1,000-bootstrap end-to-end smoke jobs on real subsets:
   Spire EX 104 rows (13.2 s), Spire IR 104 (10.5 s), FMN2 EX 95 (14.2 s),
   FMN2 IR 84 (8.1 s), two cores. All exited zero, retained the expected
   classes and wrote 512 valid probability rows. These tiny subsets selected
   no motifs and are installation checks, not biological results. A separate
   synthetic enrichment fixture exercises the plot-producing branch.

Audit inputs, logs, source revision and checks are retained locally under
`audit/rnamotifs2/` (ignored by Git). The earlier executable build/genome and
installed R packages were removed; retained audit files are source text and
analysis artifacts. No revised container image has been built or run.
