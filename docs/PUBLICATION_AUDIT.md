# Publication audit — updated 3 November 2026

## Scope and evidence

This review covered imported SPIRE RNA-seq/DESeq2/VAST/betAS outputs, the
completed FMN2 analysis, the common-event comparison, GO, Matt and rMAPS3,
the plotting/report callers, configuration, tests and Git contents. The
worktree already contained an uncommitted FASTQ reanalysis and relocation
of figures/reports. Those changes were preserved; the initial diff/status
is recorded locally under `audit/publication_2026-11-02/`.

The new publication assembly consumes the current analyses and adds nine
composed figures, individual-animal effect sensitivity, explicit division
annotation/ranking, a gene event-count-matched overlap sensitivity analysis,
complete tables and manuscript text. It does not replace the upstream
statistical methods or manufacture evidence for the requested phenotype.

## Working components

| Component | Evidence and status |
|---|---|
| Current SPIRE groups | Original manifest gives controls 1/2/7 and DKO 4/5/6; the investigator confirmed independent animals |
| betAS ingestion | Event/sample/annotation/statistic comparisons against current exports; synthetic tests reject corrupt checkpoints, changed PSI and reordered events |
| Expression | Full original DESeq2 and annotated tables agree; saved object verifies `~0 + condition`, Wald test and parametric fit |
| Genomic references | Completed RNA-seq and DGE parameters specify GRCm38/Ensembl 88; cluster defaults now agree |
| Splicing figures | Current all-sample N10 universe, unchanged effect/FPR selection, all individual PSI retained |
| Comparison | Common tested events only; shared events, genes and directional agreement distinguished; corrected FASTQ labels retained |
| Sequence/motif analysis | Completed Matt and corrected rMAPS3 outputs retained, with exploratory limits visible |
| Publication graphics | Two main plus seven supplementary figures, 180 mm width, embedded-font PDF and 450 dpi PNG, fixed color meanings |
| Release artifacts | Full CSV tables, source data, field definitions, input/output hashes, session and annotation metadata, tested ZIP packages |

## Findings that constrain the paper

- Historical SPIRE condition labels were inconsistent with original FASTQ
  identities. Old candidate lists must not be used as current results.
- PCA shows substantial animal-to-animal heterogeneity. No animals were
  excluded and no batch was inferred solely from filename dates.
- betAS's `FDR` is a per-event empirical FPR. The selected list is not
  controlled at 5% multiple-testing FDR. Exact zero does not mean certainty.
- The DESeq2 HTML template's wording does not describe the downstream
  filtered table accurately. The stored filter parameters require BH P ≤0.05
  and a twofold effect: 50 genes. The manuscript's BH-only list has 111.
- Total Spire1 expression retains signal from some exons. RNA support at
  selected target exons and DNA genotype are different measurements.
- Division-associated genes exist among the candidates, but neither the
  tested-background GO analysis nor four focused division-category tests
  establishes pathway enrichment.
- Shared SPIRE/FMN2 genes exceed uniform and event-count-stratified null
  expectations, but directional concordance is weak and no shared concordant
  event separates every raw sample in both experiments. Different controls,
  coverage and sampling units limit phenocopy claims.
- Raw motif-window minima and sparse-count asymptotic p-values do not identify
  a common causal RBP. The main motif landscape remains exploratory.
- RNAmotifs is an optional historical handoff. Its revised container recipe
  has not been built/validated and contributes no publication result.

## Validation and remaining reproducibility work

Checks cover existing thresholds, scaling, annotation/export geometry,
palette, RBP summaries, comparison joining and overlap statistics, both
SPIRE/FMN2 ingest layouts, the betAS synthetic pipeline, raw-PSI omission
sensitivity, event-count-stratified permutations, output sample labels,
file hashes, archive integrity and report links. Figure PDFs were rendered
and visually inspected for labels, clipping and legend semantics; embedded
fonts were inspected with `pdffonts`. Logs are local under `logs/publication_*`.

The assembly imports the existing fitted statistics. Alignment, complete
betAS simulations, DESeq2 fitting, sequence/motif scans and a fresh-machine
environment build were not repeated. The local wrapper depends on the
installed rig R/packages and cached container image. A public data package
and SPIRE accession are still needed for a fresh clone to reproduce the
paper. Author completion of animal/library methods is listed in the methods.

## Git release boundary

Only code, tests, configuration, environment locks and documentation belong
in the commit. Data, generated figures/tables/manuscripts, BAMs/FASTQs,
libraries, logs, audits and editor settings remain ignored. The existing
reachable Git history contains source-sized blobs, not the local data
directories. Local Git storage also contains old temporary pack files;
these are not reachable objects and are not part of a push. They were not
deleted or used as a reason to rewrite history. No remote push is performed.

## Publication redesign, 3 November

The two main figures now cover (1) target RNA, expression, splicing burden
and ranked division-candidate profiles; (2) six reference sequence features,
full distributions for two shared associations and the regional EX RBP
heatmap. Seven supplements retain QC, sensitivity, the full FMN2 comparison,
all six feature distributions, GC profiles, IR/regional motif results and
functional interpretation. Methods, results, legends and panel map were
updated together. The previous package is preserved under
`audit/publication_redesign_2026-11-02/previous_publication/`.

Rendering checks additionally compare each feature population's group size
and median with the saved tests, validate complete RBP-region matrices and
GC bins, and check candidate rankings, source data, embedded fonts and ZIPs.
The requested Lena visual style is documented in COLOURS.md; no dependency
was installed and the reference repository was read only. Fontconfig emits
compatibility warnings from the installed runtime; PDF font inspection
verifies the explicitly chosen embedded Liberation Sans face.

## Study-separated revision after manuscript review

The latest package supersedes the nine-figure composition above: two SPIRE
main figures/five supplements, one FMN2 main figure/four supplements, and
one compact convergence figure. Sample terminology replaces animal labels
in the deliverable. All 119 selected events in division-annotated genes
replace the ranked gene examples. STAR mapping QC stays in Table S1.

Main feature selection now has an explicit editorial rationale rather
than inheriting six reference features as a discovery ranking. All six
remain in supplements. MaxEnt ECDFs retain all observations and report
values outside the labelled view. Visible dendrograms document PSI and
RBP row clustering. RBP ranking/clustering is performed independently
within each study.

Sequence-level reconciliation reproduced Matt intron GC within 1e−12.
SPIRE selected IR mean GC is 48.01/48.10% versus 45.27% background; the
spatial-profile background subset is 45.41%. Separate absolute and
background-relative metaplots expose the modest shifts without imputing
the SPIRE effect to FMN2. Per-event checks and subset membership are
exported in S7. The prior release is preserved under
`audit/publication_restructure_2026-11-03/previous_publication/`.
