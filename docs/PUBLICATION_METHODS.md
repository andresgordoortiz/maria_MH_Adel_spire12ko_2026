# RNA expression and alternative-splicing analysis

## Samples and processing

SPIRE RNA-seq comprised three independent control biological replicates (Cont_1, Cont_2,
Cont_7) and three independent Spire1/2 double-knockout biological replicates (DKO_4, DKO_5,
DKO_6). Independence was confirmed by the investigator on 2 November 2026.
The original FASTQ manifest defined the groups. Historical tables contained
inconsistent condition labels and were excluded from the manuscript analyses.
All six samples were retained. Filename date labels were not treated as a
verified batch, and no batch or pairing term was fitted.

The completed single-end RNA-seq workflow used nf-core/rnaseq 3.27.0,
Nextflow 25.10.4, fastp trimming and STAR/Salmon against the GRCm38 primary
assembly and Ensembl release 88 annotation. VAST used mm10/Mm2 resources.
These are the references recorded in the completed run parameters; earlier
cluster defaults referring to GRCm39 were not the references used for these
results. The publication assembly imported completed processing outputs;
it did not realign reads.

The independent FMN2 comparison used GSE103718: two Fmn2+/− control pools
and two Fmn2−/− pools, each comprising 50 oocytes. The three technical
sequencing runs per pool were merged upstream and were not counted as
biological replicates. The accession-to-pool mapping is documented in
[CLUSTER.md](CLUSTER.md). These independent studies have different controls
and sampling units and do not constitute a factorial experiment.

## Gene expression

Differential expression was imported from nf-core/differentialabundance
2.0.0 (revision 30ed774), DESeq2 1.34.0 and R 4.1.3. The saved DESeqDataSet
directly verifies the design `~0 + condition`, a Wald test and a parametric
dispersion fit. The contrast is DKO minus Cont; this condition-only contrast
is equivalent to the corresponding contrast under an intercept model.
The pipeline received Salmon gene counts and effective gene lengths;
the saved parameters specify ratio size-factor estimation, independent
filtering, BH adjustment, a zero LFC test threshold and `alpha=0.1` for
DESeq2 independent-filter optimization. Exported shrunken log2 fold changes
and the original p-values were retained. Shrinkage was enabled; the pinned
pipeline documents ashr as its implemented shrinkage method, and ashr
2.2-54 appears in the run session information.

For the manuscript, selection used exported BH-adjusted P ≤0.05, without
an additional effect-size cutoff. The stricter subset with |shrunken
log2FC| ≥1 was reported separately. The pipeline's filtered export uses
this additional twofold threshold; it must not be confused with the full
DESeq2 table. Missing adjusted p-values were retained in supplementary
tables and excluded from selection. No DESeq2 model was refitted for
figure production. Normalized counts and VST values came from the supplied
outputs. Expression PCA used the 500 most variable VST genes after
excluding Spire1/2, with centered, unscaled variables.

## Alternative splicing and sensitivity

The completed betAS 1.2.1 analysis (revision
cd37737c75ab17a64652baf2fdb00d1b7c67c294) required corrected inclusion plus
exclusion support ≥10 and finite PSI in every sample, using six SPIRE
samples or four FMN2 pools. The run used 500 fitting points, 1,000
simulations, the package's fixed seed (`seed=TRUE`, internal seed 21122023)
and `CoverageWeight=FALSE`. All event IDs, sample PSI and exported
statistics were checked against the completed SPIRE export; ingest tests
also verify checkpoint, sample-order and annotation integrity.

Events were selected at empirical per-event FPR ≤0.05 and absolute fitted
ΔPSI ≥0.10 (10 percentage points). The upstream column named `FDR` is
an empirical **per-event false-positive rate**, not a BH-adjusted FDR
across events. It is not converted to a multiple-testing p-value here.
Zero estimates have finite simulation resolution. Consequently these
are selected candidate events rather than a list controlled at 5%
multiple-testing FDR. Pdiff was retained but was not an additional filter.
Positive effects mean increased inclusion, alternative-site use or intron
retention in knockout. EX includes the ANN annotation module; IR, Alt5 and
Alt3 are kept distinct. Alternative-site rows at one locus can be dependent.

For every selected event, raw mean PSI differences and all leave-one-sample
raw mean differences were calculated. A direction-stable event retained
the fitted effect's sign after each omission; a stricter descriptive
criterion required ≥10 percentage points in that direction after every
omission. Complete separation required every knockout PSI to lie strictly
beyond every control PSI in the fitted direction. These checks do not
refit betAS, provide new significance tests or remove any sample. Threshold
sensitivity used FPR cutoffs 0.01/0.025/0.05 and fitted-effect cutoffs
10/15/20 percentage points. Splicing PCA used only the current SPIRE PSI:
variance ≥25 percentage-points squared, one highest-variance event per
named gene (event-ID tie break), excluding Spire1/2, then centered,
unscaled PCA. It does not depend on historical sample matching.

## Functional annotation, GO and correction families

The primary functional analysis uses all testable genes, not the full
annotation catalogue. Original GO objects and raw tests were independently
reproduced before revising the inference. Input gene symbols were resolved
in the pinned org.Mm.eg.db by direct SYMBOL mapping, then by ALIAS only when
the alias resolved to exactly one Entrez ID. Ambiguous mappings were excluded;
candidate aliases, unresolved names and mapping routes are exported. Entrez
IDs were deduplicated. Engineered Spire1/2 loci were excluded from functional
queries and backgrounds. The event-level annotation retains the original
input gene name linked to its resolved Entrez ID and inherited GOALL evidence.

For each study and ontology, the primary universe is the set of resolved,
N10-testable Entrez genes annotated in that ontology. Term eligibility is
10–500 genes in this primary universe. Exactly the same term list is reused
for the full-annotation background and for both stricter selections: events
retaining a raw change ≥10 pp after every sample omission, and events with
strict raw-PSI group separation. The original event selection remains primary.
One-sided hypergeometric/Fisher tails test overrepresentation. BH adjustment
pools BP, CC and MF within each study/background/selection analysis, including
zero-overlap terms. No correction is restricted to the terms chosen for a
figure. Full-annotation comparisons are reference diagnostics.

The main GO dot plot shows the three top-ranked terms in each ontology,
ordered by complete-family BH q, nominal P and GO ID. The x axis gives
Count/query genes annotated in that ontology; point size gives gene count
and color gives q from the full tested-background BP/CC/MF family.
This post-analysis display does not require significance or a connection
to the manuscript themes; related ontology terms are retained. All background
and selection sensitivities and nine email-motivated terms remain in complete
tables and the inference audit. The main division-event panel retains every
selected event whose host gene maps to one of four division GO categories;
annotation is distinguished from pathway overrepresentation. Original
four-term-only adjusted P values remain legacy reference results.

A second sensitivity test preserves the selected-gene count within strata
floor(log2(number of testable events per gene)). Ten thousand gene allocations
(seed 42) supply expected counts and upper-tail Monte Carlo P=(1+exceedances)/
10,001 for all nine terms in each study. BH corrects these 18 tests together.
This controls event-count opportunity, with remaining coverage, expression,
gene-length and annotation biases possible. Mapping, event opportunities,
observed/expected counts, null draws and all selection sensitivities are exported.

Expression GO is a distinct analysis: BH-selected DESeq2 genes versus genes
with finite DESeq2 padj, with the same mapping/exclusion rules, 10–500 term
sizes and one BH family across BP/CC/MF. No expression term is promoted from
its nominal P value. Descriptive PCA, heatmaps, counts and sequence profiles
are not treated as new hypothesis tests. Table S1 records every inferential
family and size. Full details and rationale accompany INFERENCE_REVIEW.md.

## Cross-study comparison

Only identically annotated common tested events are compared. Shared genes,
shared events and concordant directions are distinguished. The seven existing
one-sided Fisher overlap tests, a two-sided binomial direction test and an
event-count-stratified gene-overlap permutation test form one nine-test BH
family. The matched overlap uses floor(log2(common tested events per gene)),
10,000 permutations and seed 42. Global similarity is described by Spearman
correlation of raw mean changes. These event/locus tests are limited by
within-gene dependence and the independent studies' unequal controls,
replication and detection opportunities. No factorial interaction is estimated.

## Sequence features and exact regional motif presence

Matt 1.3.2, GENCODE M10 and GRCm38 inputs are unchanged. Unselected Matt
backgrounds are capped at 10,000 events with seed 42. Existing two-sided
Mann–Whitney tests retain BH across every exported EX/IR feature and all
three contrasts within study; displayed feature choices do not define a smaller
correction family. All six reference features remain in the exported measurements. VASTDB transcript consequences remain predictions.

The primary RBP analysis now tests regional motif presence rather than a
selected minimum window P. The original rMAPS3 asymptotic window tests and
maps are reference outputs: sparse examples can have strongly optimistic
P values, which BH alone cannot repair. The new standard-library Python
counter reads the validated native FASTA and catalogue: 114 motif variants
under 91 RBP labels. Variants of an RBP are unioned; each event contributes
one positive/negative observation for presence in any window of the fixed
regional search. No P value is used to select a variant or window.

Search geometry matches the existing rMAPS grids: 50-nt windows, 1-nt steps,
50 exon or 250 intron positions; IR body positions are scaled as in the
validated native counter. A fully contained motif occurrence, including overlapping occurrences,
must occur in at least one searched window. Non-scaled grid unions span 99 or 299 nt, subject to source
sequence length; scaled intron windows can have gaps in long introns.
Synthetic tests compare membership with every explicit window, including
both edges, scaled gaps and motif unions. Every FASTA/input is hashed, and
foreground/background event counts are checked against the prepared manifests.

One-sided exact Fisher/hypergeometric tests compare motif-positive events
in the selected direction versus its unselected background. BH is applied
once across all 91 RBPs, 29 regions spanning EX/IR/Alt5/Alt3, and both
directions: 5,278 tests per study. Main EX heatmaps and supplemental IR
heatmaps show effects and corrected dots; complete Alt5/Alt3 results remain
in Table S8, including the passing PPRC1/Alt5 association.
Event correlation, GC/length composition and catalogue specificity remain
limitations of this sequence-level association analysis. It measures neither
RBP expression nor binding or causal activity.

## Figure construction

R/ggplot2/patchwork layouts use Liberation Sans, embedded in vector PDFs,
with 450-dpi PNGs at 180-mm width. Direction is blue (#2166AC), red (#B2182B)
and background grey (#BDBDBD); study identity uses purple/green. Restrained
borders and grids follow the Lena reference. Both studies have distinct
main figure sets and a compact convergence figure, with three shared
supplementary figures in total: QC/stability, additional features and IR
sequence context.

The EX RBP heatmaps are main panels and the IR heatmaps are supplemental.
Thirty EX and fifteen IR labels per
study are ordered by their smallest regional q, then largest absolute effect
and label; all tests were corrected before this display selection. Cell
color shows log2 odds ratio with a 0.5 continuity correction for display,
saturated at ±2. Dots mark exact regional tests with BH q ≤0.05. Unsaturated
row effects are standardized for Euclidean/complete-linkage clustering;
constant rows are centered without division by zero. Trees are hidden.
Matt boxplots show quartiles/1.5-IQR whiskers and individual selected-event
measurements with horizontal jitter. Background measurements are summarized
by boxes; separate outlier symbols are suppressed, with all values retained
in statistics. White points mark raw medians, dotted lines the background
median. Length/number axes are logarithmic, GC is percent. Group sizes and
medians match saved tests. MaxEnt uses a labelled −2 to 12 view with all
outside-view counts exported; negative scores remain valid observations.

Whole-intron GC recalculated from each sequence agrees with Matt within
1e−12. Positional profiles retain reference bins (EX 68, IR 88) and a
1,500-event background. Difference profiles subtract the background bin
mean; descriptive SE=sqrt(SE_selected²+SE_background²). Ribbons ±1.96 SE
are across events, not biological-replicate or multiplicity-adjusted
uncertainty. PSI heatmaps retain every selected event, center each row on
its median, and use Euclidean/complete-linkage clustering with visible
trees and saturation at ±50 pp. Expression heatmaps are mean-centered VST,
clustered the same way and saturated at ±3. Sample columns retain manifest
order; no unreported exclusions or sample clustering are applied.

## References and author completion

- Love MI, Huber W, Anders S. DESeq2. *Genome Biology* 2014;15:550.
  [Method and reference](https://doi.org/10.1186/s13059-014-0550-8).
- [DESeq2 documentation](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
  distinguishes model tests, independent filtering and LFC shrinkage.
- [nf-core/differentialabundance 2.0.0](https://registry.nextflow.io/pipelines/nf-core/differentialabundance%402.0.0)
  documents the pinned pipeline and shrinkage behavior.
- [Original FMN2 experiment](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103718).
- Upstream splicing code, betAS/Matt/rMAPS attribution, exact resource
  hashes and additional method references: [SOURCES.md](../analysis/SOURCES.md)
  and [RMAPS3.md](../analysis/RMAPS3.md).

- [GO reference-list guidance](https://geneontology.org/docs/go-enrichment-analysis/).
- [Young et al. 2010: RNA-seq selection bias](https://pmc.ncbi.nlm.nih.gov/articles/PMC2872874/).
- [betAS FPR definition](https://pmc.ncbi.nlm.nih.gov/articles/PMC10946425/).
- [R multiple-testing documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html).

Before submission, authors must supply the mouse allele/genotyping details,
donor ages, oocyte collection/stage, RNA isolation/library preparation,
sequencer details and the SPIRE data-accession statement. These details
cannot be reconstructed from the analysis tables. No microscopy or
protein measurements were available for this computational figure set.
