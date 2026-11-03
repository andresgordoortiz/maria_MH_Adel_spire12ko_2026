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

## Feature displays, GO figures and RBP presentation

Compared box/point, violin/box/point and median/quantile displays against
the local oocyte project. Chose boxplots with selected measurements, group
sizes and raw medians for all six reference features, retaining all values
in summaries/tests and labelled MaxEnt limits. Proofs and the previous
publication package are preserved under
`audit/feature_go_revision_2026-11-03/`.

RBP trees are hidden while clustering and exported row orders are retained.
EX/IR regional distributions use the actual region names at 50 degrees.
The plot title describes motif enrichment; methods identify the exact
statistical limit: raw selected minima over overlapping windows and motif
variants have no final RBP-level multiplicity calibration.

Added SPIRE S6 and FMN2 S5 for GO: top-three-term unions in BP/CC/MF under
each background, plus nine separate manuscript-driven themes. No selection
requires significance. Saved gene-set sizes validate all “Not tested”
labels as exclusions by the original 10–500 size filter. Tested-background
GO remains nonsignificant; full-background results include terms relevant
to spindle, chromosome segregation and nuclear specks. These results are
shown together rather than suppressing either background.

## Inference and background audit after scientific review

The current revision supersedes earlier figure selections and raw motif
interpretations. Direct symbols plus unique-alias rescue raise SPIRE mapping
to 971/999 non-target selected symbols. The division annotation contains
127 events/79 names. All original GO tails were independently reproduced
before rebuilding backgrounds. Testable genes are primary; all GO
ontologies share a BH family. Term eligibility is fixed by that background
and reused across full-annotation and event-selection sensitivities. The
18 event-count-matched question tests are corrected together. Expression GO
uses the separate DESeq2 universe; no term passes BH.

The new standard-library motif counter unions catalogue variants and the
fixed native regional search, counting each event once. Explicit-window
fixtures verify geometry. Exact Fisher tests replace sparse asymptotic
minimum-window scores as primary evidence. BH covers 5278 tests/study;
PCBP2/EX and PPRC1/Alt5 pass in SPIRE, none in FMN2. Native scans are
retained. All complete families are checked independently in Python.

Main GO now shows observed/expected gene counts, and the main motif panel
shows the two supported regional associations. Full heatmaps, GO,
background diagnostics and all division events are supplemental. The
family ledger, complete tests, mapping/alias candidates, null draws and
FPR calibration diagnostics accompany the package. Previous package and
review proofs are preserved under `audit/inference_revision_2026-11-03/`.
No alignment, DESeq2 fit, betAS fit or Matt run was repeated.

## Compact publication revision after figure review

The current package supersedes the figure counts and presentation choices
above. It contains two SPIRE main figures, one FMN2 main figure, one compact
convergence figure and exactly three shared supplementary figures total.
SPIRE 1D restores all 127 division-annotated events in 79 genes. SPIRE 2A is a
standard GO dot plot: top three terms per ontology, gene ratio, count and q
from the full tested-background BH family. EX RBP matrices return to both
studies' main figures, with hidden dendrograms and exact-test BH dots.

S1 retains sample quality, target-exon support, clustered expression and both
studies' omission checks; S2 adds exon number and IR splice-site strength;
S3 retains IR RBP matrices and background-relative IR GC. GO background and
selection diagnostics, full features and omitted profiles remain in tables
and the inference review. No correction family or primary selection changes.
The previous package and render proofs are preserved under
`audit/compact_publication_2026-11-03/`. Source checks, complete-family BH,
plotted values, archives, local links, embedded fonts and seven final layouts
are checked before promotion. No upstream fit, alignment or remote push is
part of this presentation revision.

## rMAPS3 heatmap restoration

At the investigator's request, SPIRE 2D, FMN2 1E and S3A/B again show native
rMAPS3 regional raw P values. The original variant eligibility, minimum-P
aggregation and top-30 EX/top-15 IR ranking are restored. Reversed Cividis
encodes −log10 raw P, saturating at 8. Rows remain clustered with hidden
trees; exact-test BH dots and odds-ratio colors are removed from heatmaps.
Legends, results, methods, source matrices, gallery and output checks now
identify the statistic explicitly. Each displayed P/score and selected RBP
label is independently checked against the exported native rMAPS3 table.
The complete exact regional-presence analysis remains in Table S8 unchanged.
The seven-figure composition and three-supplement limit are unchanged. The
previous package is preserved in `audit/rmaps_heatmaps_2026-11-03/`.
