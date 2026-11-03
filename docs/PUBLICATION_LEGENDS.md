# Figure legends

## Figure 1. SPIRE loss alters expression and splicing, including division-associated genes

**A**, Normalized Spire1/2 counts plus one, on a logarithmic axis; each point
is one independent animal (three per condition). Residual Spire1 RNA is
interpreted with exon support in Fig. S1D. **B**, Mean normalized expression
versus shrunken log2 fold change. Red/blue mark increased/decreased genes
at BH-adjusted P ≤0.05; grey includes unselected genes and genes without
finite adjusted P. All positive-mean genes with finite effects are shown;
Spire1, Spire2 and Mt1 are labelled. **C**, Selected exon (EX), intron (IR),
alternative 5′ and alternative 3′ events by direction. Selection requires
all-sample corrected support ≥10, empirical betAS per-event FPR ≤0.05 and
|fitted ΔPSI| ≥10 percentage points. This is not multiple-testing FDR control.
**D**, Raw PSI for 15 division-associated genes, one event per gene, ranked
by minimum signed leave-one-animal raw effect, then absolute raw effect and
event ID. Division membership uses four explicit GO categories; engineered
Spire loci are excluded. Columns show all six animals in manifest order;
no sample clustering is applied. **E**, Raw PSI for the three leading
candidates, with group means as black horizontal bars. Strip labels give
event class and fitted ΔPSI. Cit: MmuALTA0004163-1/3; Gen1: MmuEX0020291;
Nek4: MmuEX0031456. These candidates do not establish pathway enrichment
or division impairment. Tables S2–S4 retain all tests and candidates.

## Figure 2. Sequence architecture and regional RBP motifs of affected events

**A**, Six features fixed by the reference analysis, including nonsignificant
results. Color shows (selected-group median − unselected-group median) /
unselected-group IQR, separately for each study and direction. Saturation
is ±1; black dots denote BH-adjusted P ≤0.05 from the existing two-sided
Mann–Whitney feature-test family within each study. A positive value means
a larger feature, not a stronger splicing effect. **B,C**, Distributions of
upstream intron length for EX events and 5′ splice-site MaxEnt score for IR
events. Bg denotes the study-specific unselected background; Down/Up denote
decreased/increased inclusion or retention. Violins have equal maximum width;
white boxes show quartiles and median, whiskers extend to 1.5 IQR. Outlier
points are suppressed but retained in the distributions and tests. Length
is plotted on a logarithmic axis. Labels give event counts and BH-adjusted
P (q) against background; events are not independent animals. **D**, Top
30 EX RBP labels ranked by maximum raw motif score across both studies,
directions and regions, with Euclidean/complete-linkage row clustering.
Within each region, columns are SPIRE increased/decreased then FMN2
increased/decreased. Annotation strips identify study and usage direction.
Cells show −log10 of the minimum raw window P across eligible motif variants
for that label/region. Motif eligibility requires P ≤0.05 somewhere in the
corresponding study/direction; absent combinations are displayed as zero.
Color saturates at 8; full values are exported. Regions proceed 5′ to 3′
along the exon neighborhood and are schematic, not drawn to nucleotide scale.
Scores are exploratory, without correction for searching motifs and windows;
sparse-count asymptotic tests can be optimistic. They do not measure RBP
expression, binding or causal activity. Tables S7–S8 contain all tests,
feature measurements and motif-region results.

## Figure S1. Sample quality and support at targeted exons

**A**, PCA of the 500 most variable VST expression genes, excluding Spire1/2.
**B**, PCA of current SPIRE raw PSI, one highest-variance event per named gene
with variance ≥25 pp², excluding Spire1/2. Both center variables without
variance scaling. Axis percentages indicate explained variance; point labels
retain original FASTQ identifiers. **C**, STAR uniquely aligned read percentage.
**D**, Corrected VAST junction support for observed Spire1/2 exon events before
all-sample N10 filtering, displayed as log10(support+1); columns 1/2/7 are controls and 4/5/6 are DKO. All animals were
retained. RNA support does not replace DNA genotyping, and PCA does not
identify a biological batch by itself. Table S1 provides the underlying QC.

## Figure S2. Splicing profiles and sensitivity

**A**, Fitted versus raw mean changes for all selected SPIRE events; diagonal
is identity. **B**, Percentages retaining the fitted sign, or a ≥10 pp shift
in that direction, after every one-animal omission. These are descriptive
checks without refitting betAS. **C**, Counts under different empirical-FPR
and fitted-effect thresholds. **D**, All 1,547 selected events, row-median
centered and ordered by event class, fitted effect and event ID. Samples
retain manifest order; colors saturate at ±50 pp. Table S3 retains every
sample value, selection threshold and sensitivity result.

## Figure S3. Shared SPIRE and FMN2 targets have heterogeneous effects

Comparison is restricted to 48,323 identically annotated events in 4,229
genes. SPIRE has three animals per group; FMN2 has two pools per group,
with heterozygous controls. **A**, Shared and study-specific selected events
and genes within the common universe. A shared gene can contain different
selected events. **B**, Binned density of raw mean ΔPSI for all common events;
grey intensity is log-scaled event count. Symbols mark the 11 shared selected
events by concordance of fitted effects. **C**, Fitted effects of all shared
events; circles denote SPIRE and triangles FMN2. Lines connect estimates,
not uncertainty bounds. Seven events are concordant and four opposing
(two-sided sign-test P=0.549). **D**, Raw PSI for those events in every
animal/pool. None of the seven concordant events strictly separates both
groups in both studies. Table S5 includes overlap tests, tested universes
and event-count-stratified permutation sensitivity. Shared targeting does
not establish global molecular phenocopy or genetic interaction.

## Figure S4. Distributions of the six reference sequence features

**A–F**, Upstream and downstream intron lengths and median exon number for
EX events; intron length, GC content and 5′ splice-site strength for IR events.
Plots include every finite measurement from the Matt populations used by
the saved tests. Each study uses its own background, capped upstream at
10,000 unselected events. Violins, boxes, n and q follow Fig. 2B,C. Length
and exon-number axes are logarithmic. Associations describe genomic
features of selected events, not knockout-induced changes to DNA sequence.

## Figure S5. GC composition across exon and intron neighborhoods

**A**, Exon-centered profiles in SPIRE and FMN2; 68 sequence bins across
five regions. **B**, Intron-centered profiles; 88 bins spanning flanking
exons, splice-site flanks and scaled intron body. Lines are event-level
mean GC%; ribbons are mean ±1.96 SE across events and are descriptive,
not animal-level uncertainty. Grey regions mark exonic sequence. The
unselected background is capped at 1,500 events. Bins follow the retained
reference implementation; sequence regions have unequal physical lengths.
Per-bin event counts, means and SEs accompany the source data.

## Figure S6. Intron RBP landscape and regional motif distributions

**A**, Top 15 IR RBP labels, using the selection, clustering and display
rules in Fig. 2D. **B,C**, EX and IR motif-score distributions for eligible
motif variants, by study, direction and region (R1–R8/R1–R5 in heatmap order).
Points represent motif variants; boxes show median and quartiles, with
1.5-IQR whiskers. The dotted line is raw P=0.05. Variants, windows and RBP
labels are correlated and are not biological replicates. No significance
stars or tests between these selected score distributions are inferred.
Different variants may supply a label's minimum in each study.

## Figure S7. Expression profiles and functional interpretation

**A**, Expression volcano; only genes with finite adjusted P and effects.
**B**, VST for all 111 selected genes, centered on each gene's six-sample
mean, ordered by shrunken log2FC, with saturation at ±3 VST units.
**C**, Expression-only/shared/splicing-only selected genes among named genes
evaluable for both assays. **D**, VASTDB v3 predicted transcript consequences
of selected SPIRE events. The annotation and observed direction must be
considered together; neither establishes protein loss. **E**, Enriched GO
BP/CC/MF term counts using full OrgDb versus tested-N10 backgrounds, at
BH P ≤0.05 and q ≤0.2. **F**, Four focused one-sided Fisher tests for division
categories, using tested named genes and BH correction across four tests.
The dashed line is adjusted P=0.05. Tables S4, S6 and S7 preserve complete
membership, evidence codes, tests and predicted consequences.
