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

## Functional interpretation and study comparison

Division relevance was defined by membership in GO:0051301 (cell division),
GO:0051321 (meiotic cell cycle), GO:0007051 (spindle organization) or
GO:0007059 (chromosome segregation), including inherited GOALL membership
in the pinned mouse OrgDb. Evidence codes and all matching genes were
exported. All selected events in genes annotated to these categories are
shown, excluding the engineered Spire1/2 loci. No top-gene ranking or
extreme-profile selection is used. The units are splicing events; one gene
can contain multiple events and one event can belong to multiple categories.
Four targeted one-sided Fisher
tests used all tested named SPIRE genes as background, with BH adjustment
across these four tests. GO membership is not evidence that an isoform
change impairs division.

The existing clusterProfiler enrichGO results were retained, using mouse
org.Mm.eg.db, BH adjustment, p/adjusted-p cutoffs 0.05, q cutoff 0.2 and
gene-set sizes 10–500. Both the full OrgDb and tested-N10 backgrounds were
reported; comparison analyses use a common tested Entrez background.
All tested-term results, including nonsignificant terms, accompany the
selected results. Package versions and annotation metadata are exported.

Cross-study analysis required matching event IDs and gene/coordinate/event
annotations and used only the common tested universe. Shared gene
targeting, shared events and concordant fitted effects were distinguished.
Uniform event/locus/gene overlap used one-sided Fisher tests; the existing
seven reported overlap tests have a separate BH correction. Direction
agreement used an exploratory two-sided binomial test against 0.5.
Global similarity was described by Spearman correlation of raw mean PSI
differences across all common events. An additional gene-level sensitivity
analysis fixed SPIRE membership and permuted FMN2 selection within strata
defined by floor(log2(common tested events per gene)), preserving the
FMN2 selected count in each stratum. With 10,000 permutations and seed 42,
the upper-tail Monte Carlo P was (1 + exceedances)/(10,001). This controls
coarsely for event-count opportunity, not all gene-length, expression or
coverage biases. Dependent events limit the binomial and event-level tests.

## Exploratory annotations and presentation

Sequence features came from the pinned Matt 1.3.2 container, GENCODE M10
and GRCm38; unselected backgrounds were capped at 10,000 with seed 42.
Existing two-sided Mann–Whitney tests and their BH correction across the
exported feature-test family were retained. Length ratios compare genomic
features of selected and unselected events, not changes in intron length.
VASTDB PROT_IMPACT v3 supplies predicted transcript consequences; it does
not establish protein depletion. rMAPS3 scans used the corrected window
counter, 114 motifs under 91 RBP labels, and directional Mann–Whitney
window tests. Shared labels required raw P ≤0.05 in the same region and
direction in each study. Searching motifs, variants and windows without
family-wide calibration makes this exploratory; sparse-count asymptotic
p-values can be optimistic. No causal regulator is inferred.

Figures were composed in R/ggplot2/patchwork using the requested hybrid
of Lena QuantSeq styling and the reference RBP heatmaps. Liberation Sans provides the installed Helvetica/Arial-compatible face;
PDFs embed the resolved fonts (explicit family selection avoids the R
runtime’s generic Helvetica alias to DejaVu). Direction uses blue (#2166AC) and red
(#B2182B); study identity uses purple (#6A3D9A) and green (#33A02C).
Fine borders, restrained horizontal grids and explicit facet labels follow
the Lena visual reference. The RBP matrices retain Cividis and regional usage annotation strips.
Thirty EX and fifteen IR labels are selected independently within each
study by maximum raw score across eligible motifs, directions and regions.
Motif eligibility requires raw P ≤0.05 somewhere within study/direction.
Each label/region/direction cell retains the maximum −log10 raw P across
eligible variants; absent combinations are displayed as zero. Unsaturated
row profiles are standardized to mean zero and unit SD for Euclidean,
complete-linkage clustering; constant rows are centered without division
by zero. Dendrograms show this ordering, while cell colors show unscaled
scores saturated at 8. This separates regional profile similarity from
absolute score magnitude. No RBP expression or binding is inferred.

The six reference Matt features are all retained in study-specific
supplements. Main sequence panels focus on directly interpretable
associations: SPIRE EX upstream-intron length and IR whole-intron GC;
FMN2 EX upstream/downstream-intron lengths. This editorial selection is
post-analysis, not a new hypothesis-testing family. Annotation-dependent
co-occurrence counts and all other feature tests remain in Table S7.
Empirical cumulative distribution functions (ECDFs) show overlapping
populations without smoothing choices. Group sizes and medians are checked
against the saved tests. Feature-family BH-adjusted P values are labelled
q. Length/number axes are logarithmic; GC is shown in percent. MaxEnt
ECDFs use a −2 to 12 display window, with counts outside that window
reported on the plot and in source data. All values contribute to the
ECDF and tests; negative MaxEnt scores are retained.

Whole-intron GC was recalculated as (G+C)/(A+C+G+T) from each Matt
SEQ_INTRON and checked against INTRON_GCC at tolerance 1e−12. The audit
exports each event, its GC values and membership in the spatial-profile
subset. GC profiles retain reference bins (EX 68, IR fine 88) and the
1,500-event unselected background. Absolute plots and difference-from-
background plots are shown separately for each study. For each bin,
ΔGC is the selected-group mean minus background mean, and descriptive SE
is sqrt(SE_selected² + SE_background²). Ribbons are ±1.96 SE across
events, not biological-replicate uncertainty or multiplicity-adjusted
positional tests. Whole-intron GC and scaled local bins are different
summaries. Regional motif boxes are descriptive only.

All figures are 180 mm wide, with vector PDF layouts and 450 dpi PNGs.
Every selected PSI row is retained and centered on its sample median,
then clustered using Euclidean distance and complete linkage. The
splicing heatmaps show dendrograms and saturate at ±50 pp. The expression
supplement is gene-mean-centered VST, clustered with the same distance
and linkage, and saturated at ±3. Sample columns retain manifest order;
no sample clustering or unreported exclusion is applied. Row orders and
individual measurements are exported. Main figure sets are separate for
SPIRE and FMN2, with a compact common-universe convergence figure.

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

Before submission, authors must supply the mouse allele/genotyping details,
donor ages, oocyte collection/stage, RNA isolation/library preparation,
sequencer details and the SPIRE data-accession statement. These details
cannot be reconstructed from the analysis tables. No microscopy or
protein measurements were available for this computational figure set.
