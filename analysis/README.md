# Spire splicing figures

This is the splicing-only adaptation of the paper's code at commit `e231983fd82736971aef7fc938cf34bed5ce3121`. See [SOURCES.md](SOURCES.md), [panels.tsv](panels.tsv) and the preserved upstream [MIT notice](LICENSE.reference).

`scripts/` runs cluster simulations. `analysis/` consumes results and produces figures. `input_tables/legacy/` preserves historical data; `audit/scripts/` preserves the retired audit helpers locally. Data, downloaded reference code, libraries, logs and figures stay out of Git. Active simulation filenames, configuration and the root lockfile are unchanged.

## Reproduce the current provisional analysis

Run from the repository root with the existing filtered `data/processed/spire_N10.rds` and the three complete legacy CSVs present:

```bash
# Extra GO/plot dependencies; isolated from the simulation library.
scripts/run_R.sh analysis/setup.R

# First pass generates core figures and the full event list.
scripts/run_R.sh analysis/run.R legacy results/spire_legacy_initial

# Download and verify VASTDB + GENCODE M10 resources; never change the input table.
python3 analysis/prepare_references.py results/spire_legacy_initial/tables/all_events.csv

# Add annotations, GO and validated external-tool inputs.
scripts/run_R.sh analysis/run.R legacy results/spire_legacy_annotated

# Same Matt container as the paper. Requires rootless Podman, or CONTAINER_RUNTIME=docker.
bash analysis/run_matt.sh results/spire_legacy_annotated data/processed/spire_matt_legacy

# Final figures, including all significant Matt features and GC gradients.
MATT_OUTPUT_DIR=data/processed/spire_matt_legacy \
  scripts/run_R.sh analysis/run.R legacy results/spire_legacy_paper
```

Each output directory must be new. Open `OUTPUT/index.html` for the figure gallery; `figures/` holds PDFs and 300 dpi PNGs, and `tables/` holds plotted values, event lists, tests, exclusions and provenance. `COMPLETE` means the requested run finished; `tables/analysis_status.csv` records pending external steps. A completed core-only run can still have pending annotations: check that status table. The analysis does not submit a simulation or recompute historical statistics.

The additional analysis library is under `analysis/renv/library/` and has its own `analysis/renv.lock`. It does not change the running betAS library. Reference downloads need roughly 5 GB on disk (compressed files plus the genome/GTF and subset annotations); the Matt container requires additional space. Matt's own reports are retained alongside the feature tables. Matt upstream filenames must be repository-relative without spaces.

## Use newly generated cluster statistics

After the simulation export has its `COMPLETE` marker:

```bash
scripts/run_R.sh analysis/run.R \
  results/spire_N10_nsim1000_run01/csv/spire_all_events.csv \
  results/spire_new_initial
bash analysis/run_matt.sh results/spire_new_initial data/processed/spire_matt_new
MATT_OUTPUT_DIR=data/processed/spire_matt_new \
  scripts/run_R.sh analysis/run.R \
  results/spire_N10_nsim1000_run01/csv/spire_all_events.csv results/spire_new_paper
```

Keep the matching filtered RDS locally. The loader validates every sample PSI and annotation against it, checks complete event order, and rejects unfinished new runs. Matt outputs are reused only if both exact exported input hashes match; regenerate them when event classifications change. Never reuse old RBP outputs with new inputs unless their hashes match. FMN2 simulation support remains intact; this figure analysis intentionally accepts Spire only.

## Methods and limits

- Cont_a/b/c versus DKO_a/b/c, exactly as named in the inclusion table. Positive delta means higher inclusion/site use/retention in DKO. Independent animals, pooling, pairing and batches are not established by these column names.
- N≥10 and complete PSI in all six samples, as filtered by the existing workflow. Paper selection: upstream `FDR <= 0.05` and `abs(deltapsi) >= 0.10`; Pdiff is not an extra selection criterion for these figures. Delta is a fraction: 0.10 is 10 percentage points. The reference scripts use inclusive cutoffs, even where prose says “lower than”.
- betAS calls its empirical per-event false-positive rate `FDR`. It is not a BH-adjusted multiple-testing FDR. Exact zeros have finite simulation resolution, and tiny fitted changes at PSI boundaries can have FPR zero. Raw values are preserved; the volcano floor, compression and jitter affect display only. The effect-size filter matches the paper and does not calibrate the FPR.
- Legacy reconstruction requires identical statistical vectors across the three full `spire_pdiff_*.csv` files, and exact annotation/sample-PSI agreement at every identified row. ANN annotations come from the inclusion table at verified positions. The historical R/betAS revision, simulation settings and seed remain unverified. Historical “differential” and “combination” exports are not pooled into this analysis.
- Main EX counts include all `MmuEX` events, including ANN, matching the paper scripts. Counts by COMPLEX and the missing-annotation list are also exported. ANN is an annotation module, not proof that every event is constitutive.
- Miniheatmaps use median centering, MAD scaling (MAD=0→1), clipping ±1, Spearman distance and Ward.D2; two clusters. Constant PSI rows are explicitly excluded because their Spearman correlation is undefined. Large heatmaps use clipping ±2 and the reference ComplexHeatmap defaults (Euclidean/complete linkage), drawn with existing ggplot2; samples cluster within each condition. This preserves the two different clustering choices found in the source scripts.
- Sequence length and GC plots from VASTDB are explicitly separate from Matt. Missing reference IDs remain in splicing/GO results and are excluded from unavailable sequence annotations. IDs whose original coordinates differ from the downloaded reference are recorded and excluded from sequence/Matt/RBP inputs. Available PROT_IMPACT is v3; v2.3 is not available at the original URL. All annotation coverage and exclusions are exported.
- Matt uses the paper container, GENCODE vM10, GRCm38/mm10, species Mmus, 150 nt for branch-point analysis and `-notrbts`. Unselected background events are sampled to at most 10,000 with seed 42. No three-drug union is possible for Spire alone. Matt deduplicates genomic events and excludes events not found in the GTF; retained rows are checked and counted. All numeric features are tested with Mann–Whitney, `exact=FALSE`, ≥3 observations/group and ≥10 total. Paper-style stars use raw p; an additional BH column exposes sensitivity to multiple testing. Tests use full values, with clipping only in plotted coordinates.
- GC gradients reuse the repository's binning functions and sequences from Matt. The reference script accidentally extracts the last intron bases for its 5′ flank and the first bases for its 3′ flank; this is corrected and checked with an asymmetric sequence fixture. The ndiff pool is capped at 1,500 with seed 42; introns have fine/coarse layouts. Bands are descriptive mean ±1.96 SE across events, not independent biological replicate confidence intervals.
- GO uses the same `clusterProfiler::enrichGO`, mouse `org.Mm.eg.db`, BH correction, 0.05 p/adjusted-p cutoffs and default q cutoff 0.2, gene-set sizes 10–500. The primary `paper_OrgDb` universe matches the paper; `tested_N10` is a separately named sensitivity analysis. Gene mappings, full results and objects are exported. GO annotations/package versions are pinned here but the paper did not pin its OrgDb release, so exact historical GO replication is not claimed.
- Drug intersections/convergence require several contrasts and are not reproduced from a single Spire contrast. Gene expression, protein abundance/stability and microscopy panels are outside the requested splicing-only scope. No mechanistic or meiotic phenotype is inferred from these tables alone.

## RBP handoff: user intervention required

Use `OUTPUT/external_inputs/spire_TYPE_{up,down,ndiff}.txt` in rMAPS2 with mm10. Types: EX, IR, Alt5 and Alt3. These are the paper's **custom** coordinate tables, with VASTDB 1-based closed coordinates, not rMATS result files. On the minus strand the flanking exons are sorted by genomic coordinate, correcting the old export's strand ambiguity. Check the current custom-upload coordinate convention in the server interface before uploading; its help endpoint was inaccessible during preparation. Do not select a rMATS upload mode for these files.

Use 50 bp exon windows and 250 bp intron windows as stated in the paper. Background is all eligible unselected events per type (not Matt's sampled subset). Keep up/down defined relative to DKO. Use SE for EX and RI for IR; alternative-site coordinate tables are provided, but their mapping to the server's A5SS/A3SS custom format needs checking in the current interface before use. Export manifests identify every input event and counts/exclusions.

Before uploading, copy `OUTPUT/external_inputs/input_checksums.csv` into a new local result directory, e.g. `data/processed/spire_rmaps/`. Save downloaded tables as:

```
spire_rmaps/input_checksums.csv
spire_rmaps/EX/pVal.up.vs.bg.RNAmap.txt
spire_rmaps/EX/pVal.dn.vs.bg.RNAmap.txt
spire_rmaps/IR/pVal.up.vs.bg.RNAmap.txt
spire_rmaps/IR/pVal.dn.vs.bg.RNAmap.txt
```

EX expects RBP plus eight region columns; IR expects four; optional Alt5/Alt3 expect six (as in the reference analysis). Set `RMAPS_OUTPUT_DIR` to this directory and rerun `analysis/run.R` into a new output directory. The renderer checks input hashes, validates p-values, takes the minimum motif p per RBP/region and generates heatmaps and region counts. These are raw minimum p-values, not corrected motif significance. Zero p-values get a documented finite plotting floor. No external upload has been performed.

## Checks

```bash
scripts/run_R.sh tests/analysis_checks.R
scripts/run_R.sh tests/analysis_ingest.R
scripts/run_R.sh tests/cluster_smoke.R
BETAS_DATASET=fmn2 scripts/run_R.sh tests/cluster_smoke.R
```

The analysis checks cover exact threshold boundaries, invalid statistics, scaling dimensions, constant rows, both strands of coordinate export, RBP parsing with zero p-values, and the orientation of intron GC flanks. The ingest fixture also runs the new-result path and tests rejection of missing completion markers, reordered events and modified sample PSI. Generated runs validate alignment, sample identity, delta orientation, external input membership/hashes and output files.
