# Figure legends

## SPIRE Figure 1. Expression and splicing changes in Spire1/2 DKO samples

**A**, Mean normalized counts versus shrunken DESeq2 log2 fold change.
Red/blue mark increased/decreased genes at BH-adjusted P ≤0.05; grey
includes unselected genes and genes without finite adjusted P. All
positive-mean genes with finite effects are shown; Spire1/2 and Mt1 are
labelled. **B**, Selected EX, IR, Alt5 and Alt3 events by fitted direction.
EX includes ANN; Alt5/Alt3 denote alternative splice-site usage. Selection
requires corrected support ≥10 and finite PSI in every sample, empirical
betAS per-event FPR ≤0.05 and |fitted ΔPSI| ≥10 pp. This is not BH FDR
control. **C**, Every selected splicing event, with row-median-centered
PSI, Euclidean/complete-linkage row clustering and a visible dendrogram.
Colors saturate at ±50 pp. Columns retain manifest order; replicates 1–3
are controls (original IDs 1/2/7), 4–6 DKO (IDs 4/5/6). Each is an independent
biological replicate. **D**, All 127 selected events in 79 named genes
linked to four division GO categories by direct/unique-alias Entrez mapping,
excluding engineered loci. Each point is a splicing event, colored by event
class; a gene can contribute multiple events and an event can appear in
several categories. There are 207 memberships: 77 events/47 genes for cell
division, 42/26 for meiotic cell cycle, 29/17 for spindle organization and
59/36 for chromosome segregation. Jitter changes only vertical position.
These annotations do not establish pathway enrichment or impaired division.

## SPIRE Figure 2. GO context, sequence architecture and regional RBP motifs

**A**, Conventional GO dot plot of the three top-ranked terms in each of
BP, CC and MF, ordered by complete-family BH q, then nominal P and GO ID.
Gene ratio is the number of selected Entrez genes with the term divided by
all selected genes annotated in that ontology; point size gives gene count,
and color gives BH q. The background comprises all N10-testable,
ontology-annotated genes. Eligible terms have 10–500 background genes;
BH pools every eligible BP/CC/MF test before selecting terms to display.
No primary term passes BH q ≤0.05. Terms are ranked without requiring a
connection to division or the other manuscript themes. **B**, EX upstream
intron length. **C**, IR whole-intron GC. Boxes show quartiles/1.5-IQR
whiskers; points are individual selected events with horizontal jitter.
Background events are summarized by boxes; separate outlier symbols are
suppressed. White points give raw medians; the dotted line gives the
background median. All observations enter tests and summaries. Length uses
a logarithmic axis. q is BH across all Matt feature tests and contrasts
within study. These are genomic contexts of affected events, not DNA
sequence changes; measurements are events, not biological replicates.
**D**, rMAPS3 motif P-value heatmap for 30 EX RBP labels. Colors show
−log10 of the minimum raw rMAPS3 P value across windows and eligible motif
variants within each region/direction; colors saturate at 8. A motif variant
is eligible within study/direction when its raw P ≤0.05 in at least one region,
following the original rMAPS3 panel selection. Each RBP is ranked by its
largest score across regions/directions, with label tie breaking; the top
30 are retained. Rows are clustered on standardized unsaturated scores with
Euclidean distance and complete linkage; trees are hidden. Scores use a
1e−300 floor for numerical zeros; missing eligible profiles are shown as
P=1/score=0. The displayed P values are unadjusted selected minima, not
BH-adjusted RBP-level tests. Direction refers to event use, not protein
abundance. Native P values and the separate exact regional-presence analysis
with its complete BH correction are both retained in Table S8. Motif patterns
do not establish physical binding or causal RBP activity.

## FMN2 Figure 1. Independent FMN2 splicing and sequence-context analysis

**A**, Selected events by class and direction under the same criteria used
for SPIRE. **B**, All 275 selected events, row-centered PSI and visible
Euclidean/complete-linkage dendrogram. Replicates 1–2 are Fmn2 heterozygous
controls and 3–4 knockout; each biological replicate is a pool of 50 oocytes,
with technical sequencing runs merged upstream. **C,D**, EX upstream and
downstream intron lengths, with study-specific background and feature-family
q values; display conventions follow SPIRE Fig. 2B,C. **E**, Regional EX
rMAPS3 P values for 30 RBP labels, using the same raw-P display, ranking
and hidden-tree clustering as SPIRE Fig. 2D. The separate exact regional
presence tests in Table S8 yield no BH-significant FMN2 association.
A corresponding FMN2 expression analysis was not supplied here. The studies
have different controls, replication and tested universes.

## Convergence Figure 1. Shared targeting and directional concordance

All comparisons use 48,323 identically annotated common tested events in
4,229 genes. **A**, Study-specific and shared selected events and genes;
a shared gene may contain different selected events in each study.
**B**, Binned density of raw mean ΔPSI for all common events, with log-scaled
counts per bin. Symbols mark the 11 shared selected events by concordance
of fitted effects. **C**, Fitted effects of those 11 events; circles denote
SPIRE and triangles FMN2. Lines connect estimates, not uncertainty bounds.
Seven events are concordant and four opposing (direction-test q≈0.706
in the nine-test overlap family); global raw-change Spearman correlation is
0.0065. Table S5 includes individual sample PSI, complete overlaps and
event-count-stratified gene-overlap sensitivity. Shared targeting does not
establish an identical splicing program or a global molecular phenocopy.

## Figure S1. Sample quality, target-exon support and splicing stability

**A**, SPIRE expression PCA of the 500 most variable VST genes, excluding
Spire1/2. **B**, SPIRE raw-PSI PCA, one highest-variance event per named gene
with variance ≥25 pp², excluding Spire1/2. Variables are centered, not
variance scaled; axes give explained variance. **C**, VAST corrected
inclusion-plus-exclusion junction support for target exon events before
N10 filtering, displayed as log10(support+1). This is junction support,
not gene abundance or PSI. IDs 1/2/7 are controls and 4/5/6 DKO. RNA
measurements do not replace DNA genotyping. **D**, All 111 selected SPIRE
expression genes, mean-centered VST, with Euclidean/complete-linkage row
clustering, hidden row labels and saturation at ±3. Columns retain sample
order. **E,F**, Fractions of selected SPIRE and FMN2 events retaining the
fitted direction, or ≥10 pp in that direction, after every one-sample omission.
These descriptive checks do not refit betAS; FMN2 has only two pools per
condition. All samples were retained. STAR mapping QC, normalized target
gene counts, fitted/raw effects and threshold sensitivities remain in tables.

## Figure S2. Additional exon architecture and intron splice-site context

**A,B**, EX median exon number in SPIRE and FMN2. **C,D**, IR 5′ splice-site
MaxEnt scores in SPIRE and FMN2. Box/point conventions follow SPIRE Fig. 2B,C,
with each study's own background and complete Matt BH family. Group sizes
are events; nonsignificant results are retained. Exon-number axes are
logarithmic. MaxEnt scores are displayed from −2 to 12; captions/source data
report observations outside that view. All values enter summaries and tests;
negative scores are valid. Small FMN2 IR groups limit interpretation.
All six reference features and every Matt test remain in Table S7.

## Figure S3. Intron rMAPS3 motif P values and background-relative GC profiles

**A,B**, IR rMAPS3 P-value heatmaps for SPIRE and FMN2, respectively,
with 15 RBP labels per study. Raw-P selection, −log10 color scale, ranking
and hidden-tree row clustering follow SPIRE Fig. 2D. These colors do not
encode odds ratios or adjusted P values. **C,D**, SPIRE and FMN2 IR-centered GC
profiles after subtracting the corresponding background-bin mean. The 88
bins scale different physical regions; grey shading marks exonic regions.
The unselected background is capped at 1,500 events. Ribbons are ±1.96
sqrt(SE_selected²+SE_background²) across events. These are descriptive
pointwise ribbons, not biological-replicate uncertainty or adjusted
positional tests. Background subtraction is within study and bin.
Whole-intron GC recalculated from sequence agrees with Matt to 1e−12.
Whole-intron distributions and positional bins summarize different things;
Table S7 includes the sequence audit, subset membership, absolute profiles
and complementary EX profiles. The FMN2 profile does not establish the
SPIRE whole-intron GC association in FMN2.
