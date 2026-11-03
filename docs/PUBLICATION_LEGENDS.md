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
Colors saturate at ±50 pp. Columns retain manifest order; replicate 1–3
are controls (original IDs 1/2/7), 4–6 DKO (IDs 4/5/6). Each is an independent
biological replicate. **D**, All selected events in genes annotated to four
named division GO categories, excluding engineered Spire loci. Each point
is a splicing event, not a gene-expression measurement. An event may
appear in multiple categories; event/gene counts are per category. Small
vertical jitter separates points without changing fitted effects. The
union comprises 119 events in 75 genes. Gene annotation does not establish
pathway enrichment or an isoform's function. Tables S2–S4 retain every
measurement, event ID and event-to-gene-to-GO assignment.

## SPIRE Figure 2. Sequence context and exploratory regional RBP motifs

**A**, Upstream intron length for selected EX events and the unselected
Matt background. **B**, Whole-intron GC for IR events. ECDFs show the
cumulative fraction of events, retaining every finite value; length uses
a logarithmic axis. Down/Up denote decreased/increased exon inclusion or
intron retention. Labels give existing feature-family BH-adjusted P (q)
against background; events are not biological replicates. These are genomic
features of affected events, not changes to DNA sequence after knockout.
**C**, Top 30 EX RBP labels within SPIRE, ranked by maximum raw motif score.
Columns proceed through the annotated regions, increased then decreased
usage within each region. Cells are the maximum −log10 raw window P over
eligible motif variants; eligibility requires P ≤0.05 somewhere in that
study/direction. Absent combinations display zero. Dendrograms use
Euclidean/complete linkage on row-standardized unsaturated score profiles;
cell colors retain raw scores and saturate at 8. Regions are schematic,
not nucleotide-scaled. Searching variants/windows and sparse-hit asymptotic
tests make these scores exploratory. They do not measure RBP expression,
binding or causality. Tables S7–S8 retain complete results.

## SPIRE Figure S1. RNA quality control and targeted-locus support

**A**, PCA of the 500 most variable VST expression genes, excluding Spire1/2.
**B**, PCA of raw PSI, one highest-variance event per named gene with
variance ≥25 pp², excluding Spire1/2. Variables are centered, not variance
scaled; axes give explained variance. **C**, Normalized Spire1/2 gene
counts plus one on a logarithmic axis; each point is one sample. **D**,
VAST corrected inclusion-plus-exclusion junction support for target exon
events before N10 filtering, displayed as log10(support+1). This is junction
support, not gene abundance or PSI. IDs 1/2/7 are controls and 4/5/6 DKO.
All samples were retained. RNA measurements do not replace DNA genotyping.
STAR mapping percentages are provided separately as alignment QC in Table S1.

## SPIRE Figure S2. Sensitivity and clustered expression profiles

**A**, Fitted versus raw mean PSI effects; diagonal is identity. **B**,
Fraction of selected events retaining the fitted direction, or ≥10 pp in
that direction, after every one-sample omission. This does not refit betAS.
**C**, Counts under the indicated empirical-FPR and fitted-effect cutoffs.
**D**, All 111 selected expression genes, mean-centered VST, with
Euclidean/complete-linkage row clustering and color saturation at ±3.
Columns retain sample order. Raw measurements and row order are exported.

## SPIRE Figure S3. All six reference sequence features

**A–F**, ECDFs of EX upstream/downstream intron lengths and median exon
number; IR length, GC and 5′ splice-site MaxEnt score. The reference feature
set includes nonsignificant results. Group populations match the saved
Matt tests. Length/number axes are logarithmic; GC is percent. MaxEnt scores
are displayed from −2 to 12; the caption and source data report values
outside that view. All values contribute to ECDFs and tests. Negative
MaxEnt scores are valid and are not truncated in the data. q denotes
existing BH adjustment across the study's full feature-test family.

## SPIRE Figure S4. Absolute and background-relative sequence GC

**A,B**, Absolute EX-centered mean GC and its difference from the background.
**C,D**, Corresponding IR-centered profiles. Grey shading marks exonic
regions. EX uses 68 bins and IR 88; local regions differ in physical length.
The unselected background is capped at 1,500 events. Ribbons are mean ±1.96
SE across events; for differences, SE is sqrt(SE_selected²+SE_background²).
These are descriptive pointwise ribbons, not biological-replicate uncertainty
or multiplicity-adjusted positional tests. Background subtraction is within
study and bin. Whole-intron GC was independently recalculated from each
sequence and agreed with Matt to tolerance 1e−12. Whole-intron distributions
and scaled positional bins are different summaries; Table S7 records the
sequence audit and profile-subset membership.

## SPIRE Figure S5. Intron RBP landscape and regional motif distributions

**A**, Top 15 IR RBP labels using the within-study selection and clustering
rules of SPIRE Fig. 2C. **B,C**, EX and IR motif-score distributions across
regions R1–R8 and R1–R5 in heatmap order. Boxes show median/quartiles and
1.5-IQR whiskers, with outliers retained as points. The dotted line marks
raw P=0.05. Motif variants and regions are correlated and are not biological
replicates; no tests between the preselected score distributions are inferred.

## FMN2 Figure 1. Independent FMN2 splicing and sequence-context analysis

**A**, Selected events by class and direction under the same criteria used
for SPIRE. **B**, All 275 selected events, row-centered PSI and visible
Euclidean/complete-linkage dendrogram. Replicates 1–2 are Fmn2 heterozygous
controls and 3–4 knockout; each biological replicate is a pool of 50 oocytes,
with technical sequencing runs merged upstream. **C,D**, ECDFs of EX
upstream/downstream intron lengths, with study-specific background and
feature-family q values. **E**, Top 30 EX RBP labels selected and clustered
within FMN2, following SPIRE Fig. 2C. There is no gene-expression panel
because a corresponding FMN2 expression analysis was not supplied here.
The experiments have different controls, replication and tested universes.

## FMN2 Figure S1. Sensitivity across biological replicate pools

**A**, Fitted versus raw mean effects for all selected events. **B**,
Fractions retaining direction or ≥10 pp after every one-pool omission.
These descriptive checks use two pools per condition and do not refit
betAS or provide calibrated significance tests.

## FMN2 Figure S2. Reference sequence-feature distributions

**A–F**, The same six features and ECDF conventions as SPIRE Fig. S3, using
FMN2 populations and its own background and BH family. All observations
remain in the calculations, including MaxEnt values outside the labelled
view. Small IR groups limit the interpretation of overlapping distributions.

## FMN2 Figure S3. Absolute and background-relative GC profiles

**A–D**, The same layouts and uncertainty conventions as SPIRE Fig. S4,
computed independently within FMN2. Wider event-level ribbons reflect its
smaller selected sets. These plots do not establish the SPIRE GC association
in FMN2. Sequence-level checks and subset summaries accompany Table S7.

## FMN2 Figure S4. Intron motifs and regional score distributions

**A**, Top 15 IR labels ranked and clustered within FMN2. **B,C**, EX and IR
regional motif-score distributions. Selection, dendrogram, saturation and
exploratory interpretation follow SPIRE Fig. S5.

## Convergence Figure 1. Shared targeting and directional concordance

All comparisons use 48,323 identically annotated common tested events in
4,229 genes. **A**, Study-specific and shared selected events and genes;
a shared gene may contain different selected events in each study.
**B**, Binned density of raw mean ΔPSI for all common events, with log-scaled
counts per bin. Symbols mark the 11 shared selected events by concordance
of fitted effects. **C**, Fitted effects of those 11 events; circles denote
SPIRE and triangles FMN2. Lines connect estimates, not uncertainty bounds.
Seven events are concordant and four opposing (exploratory sign-test
P=0.549); global raw-change Spearman correlation is 0.0065. Table S5 includes
individual sample PSI, complete overlaps and event-count-stratified gene
overlap sensitivity. Shared targeting does not establish an identical
splicing program or a global molecular phenocopy.
