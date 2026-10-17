# Spire and FMN2 splicing figures

This is the splicing-only adaptation of the paper's code at commit `e231983fd82736971aef7fc938cf34bed5ce3121`. See [SOURCES.md](SOURCES.md), [panels.tsv](panels.tsv) and the preserved upstream [MIT notice](LICENSE.reference).

`scripts/` runs cluster simulations. `analysis/` consumes completed results and produces figures. `final_input_tables/` holds the unchanged inclusion tables and final simulation exports. `legacy/` preserves provisional results, historical tables, cached Spire Matt extraction and the incomplete Spire run. Its dated archive manifest records original paths and SHA-256 hashes. `final_input_tables/inclusion_manifest.csv` records the moved inclusion tables. Data, reference downloads, libraries, logs and figures stay out of Git. Cluster input paths, simulation configuration and the root lockfile are unchanged.

## Reproduce the final analyses

Run from the repository root with the matching filtered objects at `data/processed/spire_N10.rds` and `data/processed/fmn2/fmn2_N10.rds`. The local raw inputs are `final_input_tables/spiredko_INCLUSION_LEVELS_FULL-mm10.tab` and `final_input_tables/fmndko_INCLUSION_LEVELS_FULL-mm10.tab`.

```bash
# Extra GO/plot dependencies; isolated from the simulation library.
scripts/run_R.sh analysis/setup.R

# Download and verify VASTDB + GENCODE M10 resources; never change the input table.
python3 analysis/prepare_references.py \
  final_input_tables/spire_N10_nsim1000_run02/csv/spire_all_events.csv \
  final_input_tables/fmn2_N10_nsim1000_run01/csv/fmn2_all_events.csv

# Validate completed statistics and export sequence-tool inputs.
scripts/run_R.sh analysis/prepare_inputs.R \
  final_input_tables/spire_N10_nsim1000_run02/csv/spire_all_events.csv data/processed/spire_final_inputs
BETAS_DATASET=fmn2 scripts/run_R.sh analysis/prepare_inputs.R \
  final_input_tables/fmn2_N10_nsim1000_run01/csv/fmn2_all_events.csv data/processed/fmn2_final_inputs

# Same Matt container as the paper. Requires rootless Podman, or CONTAINER_RUNTIME=docker.
bash analysis/run_matt.sh data/processed/spire_final_inputs data/processed/spire_matt_final
bash analysis/run_matt.sh data/processed/fmn2_final_inputs data/processed/fmn2_matt_final

# Independent studies can run concurrently; each writes a separate directory.
scripts/run_R.sh analysis/run.R \
  final_input_tables/spire_N10_nsim1000_run02/csv/spire_all_events.csv results/spire_final &
BETAS_DATASET=fmn2 scripts/run_R.sh analysis/run.R \
  final_input_tables/fmn2_N10_nsim1000_run01/csv/fmn2_all_events.csv results/fmn2_final &
wait

# Requires both completed study analyses, including their Matt statistics.
scripts/run_R.sh analysis/compare.R results/spire_final results/fmn2_final results/comparison
```

Each output directory must be new. Open `OUTPUT/index.html` for the figure gallery; `figures/` holds PDFs and 300 dpi PNGs, and `tables/` holds plotted values, event lists, tests, exclusions and provenance. `COMPLETE` means the requested run finished; `tables/analysis_status.csv` records pending external steps. The analysis does not submit a simulation or recompute betAS statistics. For the present Spire run, `MATT_OUTPUT_DIR=legacy/matt/spire_legacy` reuses the archived extraction after verifying exact input hashes and retained event values; downstream tests and figures are regenerated.

The additional analysis library is under `analysis/renv/library/` and has its own `analysis/renv.lock`. It does not change the running betAS library. Reference downloads need roughly 5 GB on disk (compressed files plus the genome/GTF and subset annotations); the Matt container requires additional space. Matt's own reports are retained alongside the feature tables. Matt upstream filenames must be repository-relative without spaces.

Keep the matching filtered RDS locally. The loader validates every sample PSI and annotation against it, checks complete event order, and rejects unfinished runs. If serialized filtered-object hashes differ between machines, it also requires and validates the completed statistics checkpoint, its run metadata, sample sheet and supplied environment lockfile; both hashes and the validation method are exported. This confirms exported PSI/statistics, but does not compare the unavailable cluster-side Qual object or establish why the RDS bytes differ. Matt outputs are reused only if both exact exported input hashes match; regenerate them when event classifications change. Never reuse old RBP outputs with new inputs unless their hashes match. Set `BETAS_DATASET=fmn2` for FMN2 ingestion and rendering; otherwise Spire is selected.

## Compare the studies

`compare.R` matches EVENT IDs and requires matching gene, coordinate and event annotations. Event and gene overlap, selection rates, fitted delta plots and raw-sample PSI plots use only the common tested events. Alternative-site suffixes are also collapsed to loci in a separate overlap test. Exported candidates distinguish concordant events, opposing events and shared affected genes whose selected events differ.

The comparison runs BP/CC/MF GO for each study, the shared genes and the concordant-event genes on one common-tested Entrez background. Shared GO terms include their overlapping driver genes. Matt comparisons retain every tested feature and report both raw-p and BH significance, effect direction and median shifts scaled by the dataset's ndiff IQR. Full-dataset protein-impact proportions remain descriptive because the tested universes differ.

Overlap Fisher tests use a uniform event/locus/gene null. Events within genes and complementary alternative-site rows are dependent; gene length and event-count biases are not fully corrected. Direction agreement is reported separately, with an exploratory sign test. These independent experiments cannot estimate a genetic interaction, and a missing selection is not evidence of no effect.

## Methods and limits

- Spire: Cont_a/b/c versus DKO_a/b/c, exactly as named in the inclusion table. Independent animals, pooling, pairing and batches are not established by these column names. FMN2: two merged Het control pools versus two KO pools; these are four profiles, not twelve biological replicates. See the root README for the GEO mapping. Positive delta means higher inclusion/site use/retention in the knockout group.
- N≥10 and complete PSI in every sample (six for Spire, four for FMN2), as filtered by the existing workflow. Paper selection: upstream `FDR <= 0.05` and `abs(deltapsi) >= 0.10`; Pdiff is not an extra selection criterion for these figures. Delta is a fraction: 0.10 is 10 percentage points. The reference scripts use inclusive cutoffs, even where prose says “lower than”.
- betAS calls its empirical per-event false-positive rate `FDR`. It is not a BH-adjusted multiple-testing FDR. Exact zeros have finite simulation resolution, and tiny fitted changes at PSI boundaries can have FPR zero. Raw values are preserved; the volcano floor, compression and jitter affect display only. The effect-size filter matches the paper and does not calibrate the FPR.
- Legacy reconstruction requires identical statistical vectors across the three full `spire_pdiff_*.csv` files, and exact annotation/sample-PSI agreement at every identified row. ANN annotations come from the inclusion table at verified positions. The historical R/betAS revision, simulation settings and seed remain unverified. Historical “differential” and “combination” exports are not pooled into this analysis.
- Main EX counts include all `MmuEX` events, including ANN, matching the paper scripts. Counts by COMPLEX and the missing-annotation list are also exported. ANN is an annotation module, not proof that every event is constitutive.
- Events without gene labels remain in event counts and plots, but are excluded from gene lists, gene counts and GO input-gene counts. Their event counts are exported in `gene_annotation_coverage.csv`.
- Miniheatmaps use median centering, MAD scaling (MAD=0→1), clipping ±1, Spearman distance and Ward.D2; two clusters. Constant PSI rows are explicitly excluded because their Spearman correlation is undefined. Large heatmaps use clipping ±2 and the reference ComplexHeatmap defaults (Euclidean/complete linkage), drawn with existing ggplot2; samples cluster within each condition. This preserves the two different clustering choices found in the source scripts.
- Sequence length and GC plots from VASTDB are explicitly separate from Matt. Missing reference IDs remain in splicing/GO results and are excluded from unavailable sequence annotations. IDs whose original coordinates differ from the downloaded reference are recorded and excluded from sequence/Matt/RBP inputs. Available PROT_IMPACT is v3; v2.3 is not available at the original URL. All annotation coverage and exclusions are exported.
- Matt uses the paper container, GENCODE vM10, GRCm38/mm10, species Mmus, 150 nt for branch-point analysis and `-notrbts`. Unselected background events are sampled to at most 10,000 with seed 42. No three-drug union is possible for Spire alone. Matt deduplicates genomic events and excludes events not found in the GTF; retained rows are checked and counted. All numeric features are tested with Mann–Whitney, `exact=FALSE`, ≥3 observations/group and ≥10 total. Paper-style stars use raw p; an additional BH column exposes sensitivity to multiple testing. Tests use full values, with clipping only in plotted coordinates.
- GC gradients reuse the repository's binning functions and sequences from Matt. The reference script accidentally extracts the last intron bases for its 5′ flank and the first bases for its 3′ flank; this is corrected and checked with an asymmetric sequence fixture. The ndiff pool is capped at 1,500 with seed 42; introns have fine/coarse layouts. Bands are descriptive mean ±1.96 SE across events, not independent biological replicate confidence intervals.
- GO uses the same `clusterProfiler::enrichGO`, mouse `org.Mm.eg.db`, BH correction, 0.05 p/adjusted-p cutoffs and default q cutoff 0.2, gene-set sizes 10–500. The primary `paper_OrgDb` universe matches the paper; `tested_N10` is a separately named sensitivity analysis. Gene mappings, full results and objects are exported. GO annotations/package versions are pinned here but the paper did not pin its OrgDb release, so exact historical GO replication is not claimed.
- Paper drug intersections/convergence are adapted to a two-study comparison, with the different coverage and controls made explicit. Gene expression, protein abundance/stability and microscopy panels are outside the requested splicing-only scope. No mechanistic or meiotic phenotype is inferred from these tables alone.

## RBP handoff: user intervention required

Use `OUTPUT/external_inputs/DATASET_TYPE_{up,down,ndiff}.txt` in rMAPS2 with mm10 (`DATASET=spire` or `fmn2`). Types: EX, IR, Alt5 and Alt3. These are the paper's **custom** coordinate tables, with VASTDB 1-based closed coordinates, not rMATS result files. On the minus strand the flanking exons are sorted by genomic coordinate, correcting the old export's strand ambiguity. Check the current custom-upload coordinate convention in the server interface before uploading; its help endpoint was inaccessible during preparation. Do not select a rMATS upload mode for these files.

Use 50 bp exon windows and 250 bp intron windows as stated in the paper. Background is all eligible unselected events per type (not Matt's sampled subset). Keep up/down defined relative to the knockout. Use SE for EX and RI for IR; alternative-site coordinate tables are provided, but their mapping to the server's A5SS/A3SS custom format needs checking in the current interface before use. Export manifests identify every input event and counts/exclusions.

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
BETAS_DATASET=fmn2 scripts/run_R.sh tests/analysis_ingest.R
scripts/run_R.sh tests/comparison_checks.R
scripts/run_R.sh tests/cluster_smoke.R
BETAS_DATASET=fmn2 scripts/run_R.sh tests/cluster_smoke.R
```

The analysis checks cover exact threshold boundaries, invalid statistics, scaling dimensions, constant rows, both strands of coordinate export, RBP parsing with zero p-values, and the orientation of intron GC flanks. Both ingest layouts test missing completion markers, reordered events, modified PSI, differing serialized hashes and corrupt checkpoints. Comparison fixtures check the common universe, exact-ID joining, concordance, alternative-site loci and conflicting annotations. Generated runs validate alignment, sample identity, delta orientation, external input membership/hashes and output files.
