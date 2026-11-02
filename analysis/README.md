# Spire and FMN2 splicing figures

This is the splicing-only adaptation of the paper's code at commit `e231983fd82736971aef7fc938cf34bed5ce3121`. See [SOURCES.md](SOURCES.md), [panels.tsv](panels.tsv) and the preserved upstream [MIT notice](LICENSE.reference).

`scripts/` runs cluster simulations; `analysis/` consumes completed results and produces figures. Current SPIRE inputs are in `data/inputs/spire_fastq/`; `final_input_tables/` preserves historical SPIRE inputs and the unchanged FMN2 inputs. `legacy/` holds superseded analyses. The 2 November 2026 relocation manifest and historical output hashes are in `audit/spire_reanalysis_2026-11-02/`. Data, reference downloads, libraries, logs and figures stay out of Git. Cluster configuration and environment lockfiles are unchanged.


## Figure locations after cleanup

`results/index.html` is the curated figure gallery. Published PNG/PDF files are
under its five numbered folders; full reports, tables and individual rMAPS maps
are under `reports/`. Each report's `figures/` entries link to the published files,
so established plotting filenames remain usable without duplicate images.

Imported pipeline outputs are in `data/inputs/spire_fastq/`. CRG submission scripts
retain their remote `results/spire_fastq/` output locations; their environment
explicitly supplies the inclusion-table path. This local organization does not
change existing cluster runs. Original provenance paths and historical completion
commands are preserved as run; use `audit/results_cleanup_2026-11-02/relocations.csv`
for their new locations. Regenerate the viewing pages with:

```bash
python3 analysis/render_results.py
```

## Reproduce the current SPIRE analysis (2 November 2026)

The imported outputs are in `data/inputs/spire_fastq/`: RNA-seq, supplied DESeq2
results, the new VAST inclusion table and the completed N10/nsim1000 betAS run.
The original FASTQ condition mapping is `config/spire_fastq.tsv` (Cont 1/2/7,
DKO 4/5/6). FMN2 inputs and numerical results are unchanged.

Run from the repository root. Output destinations must be new. The local
wrapper uses a container, reuses installed R 4.5.2/packages and host runtime
libraries read-only, disables network access and automatic renv activation,
and grants writes only to analysis output directories. It installs nothing.
The existing CRG runner remains `scripts/run_R.sh`.

```bash
# Recreate N10 and validate every exported PSI/statistic against checkpoints.
bash scripts/run_local_R.sh scripts/01_filter.R
bash scripts/run_local_R.sh analysis/prepare_inputs.R \
  data/inputs/spire_fastq/betas_N10_nsim1000/csv/spire_fastq_all_events.csv \
  data/processed/spire_fastq/inputs

# New sequence features; source and references are mounted read-only.
bash analysis/run_matt.sh data/processed/spire_fastq/inputs data/processed/spire_fastq/matt
bash scripts/run_local_R.sh analysis/run.R \
  data/inputs/spire_fastq/betas_N10_nsim1000/csv/spire_fastq_all_events.csv reports/spire_final

# Current cross-study comparison and old/new SPIRE review, including supplied DGE.
bash scripts/run_local_R.sh analysis/compare.R reports/spire_final reports/fmn2_final reports/comparison
bash scripts/run_local_R.sh analysis/reanalyze_spire.R

```

Matt completed after the core figure run. The recorded completion command
`bash scripts/run_local_R.sh audit/spire_reanalysis_2026-11-02/complete_matt.R`
adds its validated figures/tables without repeating unchanged betAS/GO statistics.
A fresh run with Matt already complete uses the ordinary `analysis/run.R` command.
The six-panel Figure 2 was completed using the unchanged plotting modules; its
completion script, input hashes and temporary-output cleanup are recorded under
`audit/spire_reanalysis_2026-11-02/rmaps_completion/`.
Fresh SPIRE rMAPS scans use the corrected FASTQ groups. Reproduce them after the
core analysis, using a new output directory:

```bash
PYTHON_CPU_COUNT=6 pixi run --manifest-path ../rMAPS3/pixi.toml --locked python analysis/run_rmaps3.py \
  reports/spire_final spire_fastq reports/spire_final/rmaps3_outputs ../rMAPS3
pixi run --manifest-path ../rMAPS3/pixi.toml --locked python -B tests/rmaps3_output_checks.py spire_fastq
bash scripts/run_local_R.sh analysis/render_rmaps3.R
```

The renderer joins new SPIRE scans with the unchanged FMN2 scans and checks their
input hashes and dataset provenance. Historical SPIRE motif results stay archived.

Original inclusion tables, DESeq2 outputs and simulation checkpoints are retained.
The review does not rerun alignment, DESeq2 or betAS simulations.

Previous SPIRE reports, sample-identity diagnostics, comparisons and RBP figures
are archived in `legacy/pre_fastq_2026-11-02/results/`. The previous reproduction
instructions are preserved in `legacy/pre_fastq_2026-11-02/ANALYSIS_README.md`.
The three diagnostic scripts `check_spire_raw.py`, `check_spire_labels.R` and
`review_plots.R` explicitly target that historical archive. Current sample
correspondence and locus review are generated by `reanalyze_spire.R`.

`COMPLETE` marks a successful run; `tables/analysis_status.csv` records panel
coverage. Checkpoints, sample sheets, event order, PSI, annotations, contrast
orientation and environment-lock checksums are validated before figure generation.
If the local filtered RDS has different serialized bytes, exact PSI/annotation
agreement and the completed statistics/run/sample/lock checkpoints are required.
Old Matt/RBP outputs are never reused with new inputs unless their hashes match.

## Compare the studies

`compare.R` matches EVENT IDs and requires matching gene, coordinate and event annotations. Event and gene overlap, selection rates, fitted delta plots and raw-sample PSI plots use only the common tested events. Alternative-site suffixes are also collapsed to loci in a separate overlap test. Exported candidates distinguish concordant events, opposing events and shared affected genes whose selected events differ.

The comparison runs BP/CC/MF GO for each study, the shared genes and the concordant-event genes on one common-tested Entrez background. Shared GO terms include their overlapping driver genes. Matt comparisons retain every tested feature and report both raw-p and BH significance, effect direction and median shifts scaled by the dataset's ndiff IQR. Full-dataset protein-impact proportions remain descriptive because the tested universes differ.

Overlap Fisher tests use a uniform event/locus/gene null. Events within genes and complementary alternative-site rows are dependent; gene length and event-count biases are not fully corrected. Direction agreement is reported separately, with an exploratory sign test. These independent experiments cannot estimate a genetic interaction, and a missing selection is not evidence of no effect.

## Methods and limits

- Current Spire: Cont_1/2/7 versus DKO_4/5/6, exactly as named in the original FASTQ manifest. Historical Cont_a/b/c and DKO_a/b/c results are archived. The investigator confirmed three independent animals per condition on 2 November 2026; batch identities remain unverified. FMN2: two merged Het control pools versus two KO pools; these are four profiles, not twelve biological replicates. See [the cluster guide](../docs/CLUSTER.md) for the GEO mapping. Positive delta means higher inclusion/site use/retention in the knockout group.
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
- Paper drug intersections/convergence are adapted to a two-study comparison, with the different coverage and controls made explicit. These legacy modules cover splicing; the publication assembly below adds gene expression. Protein abundance/stability and microscopy measurements are not available here. No mechanistic or meiotic phenotype is inferred from these tables alone.

## RBP motif analysis

Current joint figures use the fresh SPIRE FASTQ scans and unchanged FMN2 maps.
Previous joint figures and raw SPIRE maps remain in
`legacy/pre_fastq_2026-11-02/results/` with their historical labels and selections.
See [RMAPS3.md](RMAPS3.md) for commands, methods and the maintained counter fork.

## Checks

Shared colours are defined once in `palettes.R`; see
[the palette guide](../docs/COLOURS.md). Current biological conclusions are in
[the evidence summary](../docs/CONCLUSIONS.md). Historical diagnostics remain
in `check_spire_raw.py`, `check_spire_labels.R` and `review_plots.R`; they read
and write the dated archive explicitly.

```bash
scripts/run_local_R.sh tests/analysis_checks.R
scripts/run_local_R.sh tests/palette_checks.R
scripts/run_local_R.sh tests/rbp_paper_checks.R
scripts/run_local_R.sh tests/analysis_ingest.R
BETAS_DATASET=fmn2 scripts/run_local_R.sh tests/analysis_ingest.R
scripts/run_local_R.sh tests/comparison_checks.R
scripts/run_local_R.sh tests/cluster_smoke.R
BETAS_DATASET=fmn2 scripts/run_local_R.sh tests/cluster_smoke.R
```

The analysis checks cover exact threshold boundaries, invalid statistics, scaling dimensions, constant rows, both strands of coordinate export, RBP parsing with zero p-values, and the orientation of intron GC flanks. Both ingest layouts test missing completion markers, reordered events, modified PSI, differing serialized hashes and corrupt checkpoints. Comparison fixtures check the common universe, exact-ID joining, concordance, alternative-site loci and conflicting annotations. Generated runs validate alignment, sample identity, delta orientation, external input membership/hashes and output files.

## Publication assembly

The two main and four supplementary publication layouts are generated by
`publication.R`, with the existing ggplot2/patchwork/OrgDb environment.
They add imported DESeq2 results, current-only PSI PCA, individual-animal
sensitivity, explicit division annotations and event-count-stratified
cross-study overlap sensitivity. Existing statistics and reports are preserved.

```bash
bash scripts/run_local_R.sh tests/publication_checks.R
bash scripts/run_local_R.sh analysis/publication.R results/publication
python3 analysis/package_publication.py results/publication
python3 tests/publication_output_checks.py results/publication
```

Use a new output directory for another render. Packaging uses Python's standard
library and optionally the installed Pandoc for Word/HTML manuscript export;
it installs nothing. See [the publication audit](../docs/PUBLICATION_AUDIT.md),
[methods](../docs/PUBLICATION_METHODS.md) and [legends](../docs/PUBLICATION_LEGENDS.md).
The completed package includes every plotted value, complete statistical tables,
file hashes, annotation metadata, session versions and an editable manuscript.
The local wrapper requires the already-installed rig R/packages and container
image; it is not a portable environment installer. Use the pinned lockfiles in
an isolated container to provision another machine; a clean-room rebuild was
not performed for this figure assembly. Do not run setup/renv bootstrap scripts
against host or existing analysis libraries without approval.
