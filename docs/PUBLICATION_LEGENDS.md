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
biological replicate. **D**, Compact GO context for seven processes raised
by the manuscript. Each point describes genes hosting selected splicing
events, after direct/unique-alias mapping and exclusion of engineered loci.
Filled points show observed gene counts; open points show tested-background
expectations. Connecting lines show their differences. Labels give genes and q from
BH across the full tested-background BP/CC/MF family. q tests overrepresentation. No displayed category passes
BH. All 127 division-annotated events in 79 named genes appear in S7D.

## SPIRE Figure 2. Sequence contexts and two regional motif associations

**A**, EX upstream intron length. **B**, IR whole-intron GC. Boxes show
quartiles/1.5-IQR whiskers; points are individual selected events with
horizontal jitter. Background events are summarized by boxes; separate
outlier symbols are suppressed. White points give raw medians; the dotted
line gives the background median. All observations enter tests and summaries.
Length uses a logarithmic axis. q is BH across all Matt feature tests and
contrasts in that study. Measurements are events, not biological replicates.
These are genomic contexts of affected events, not changes in DNA sequence.
**C**, The two regional motif-presence tests passing BH in SPIRE. Bars show
the fraction of events positive for any catalogue variant of the RBP in the
fixed regional search; labels give positive/total events. PCBP2 uses EX
upstream-exon 3′ regions; PPRC1 uses the Alt5 intron after the long exon.
The one-sided exact presence tests are corrected across 5,278 RBP/region/
direction tests covering all four event types. Neither a window nor a
variant is selected by P. Regional sequence association does not establish
RBP abundance, physical binding or causal activity. Full counts, geometry,
tests and clustered effect profiles accompany Table S8 and Fig. S5.

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

**A–F**, Boxplots and individual selected measurements for EX upstream/downstream intron lengths and median exon
number; IR length, GC and 5′ splice-site MaxEnt score. The reference feature
set includes nonsignificant results. Group populations match the saved
Matt tests. Length/number axes are logarithmic; GC is percent. MaxEnt scores
are displayed from −2 to 12; the caption and source data report values
outside that view. All values contribute to box summaries and tests; plotting conventions follow SPIRE Fig. 2A,B. Negative
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

## SPIRE Figure S5. Regional RBP effects and complete-family correction

**A,B**, EX and IR effect heatmaps for 30 and 15 RBP labels, respectively.
Rows are chosen by their smallest complete-family q, then largest absolute
effect and label. Colors show log2 motif-presence odds ratios with a 0.5
continuity correction for display, saturated at ±2. Dots mark exact regional
presence tests with BH q ≤0.05. Clustering uses row-standardized unsaturated
effects and Euclidean/complete linkage; trees are hidden. **C,D**, Effect
summaries across all 91 RBP labels in the named EX and IR regions. Boxes
show median/quartiles, 1.5-IQR whiskers and outliers; the dotted line is zero.
These descriptive distributions receive no additional between-group tests.
Axes use the anatomical region names at 50°. Up/down refers to event use,
not protein abundance. Complete tests and positive/total counts are in S8.

## SPIRE Figure S6. GO ontologies under two reference backgrounds

**A**, Union of the top three raw-P terms per ontology/background; all
terms were tested before display selection. **B**, Nine explicit terms
motivated by division, actin, nuclear envelope, chromatin, speckles and
splicing. Columns compare the primary tested-gene and full-annotation
backgrounds using exactly the same eligible term list, fixed by 10–500
primary-background genes. BP, CC and MF are faceted. The x axis is observed/
expected gene enrichment, size gives contributing genes, and color gives
BH-adjusted P pooled across all three ontologies. Related parent/child
terms are retained. Engineered loci are excluded; aliases are rescued only
when unambiguous. No primary tested-background term passes BH. Full-annotation
results are sensitivity/reference evidence, not a substitute primary result.

## SPIRE Figure S7. Why annotation counts and GO enrichment differ

**A**, Proportion of non-engineered named genes with a selected event,
grouped by number of testable events. The 128+ group combines higher bins
for display only; null matching retains log2 strata. **B**, Observed selected
gene counts versus three expectations: all ontology-annotated genes, tested
genes, or tested genes matched for event-count opportunity. Matching uses
10,000 allocations with fixed selected counts per stratum; all 18 question-
driven study/term tests receive BH correction together. **C**, Observed/
expected ratios, selected gene counts and full-family q across the original,
omission ≥10 pp and complete raw-PSI-separation selections. The same GO
family is retained; the primary selection does not change. Colors saturate
at 0–2. **D**, All 127 selected events in 79 named genes linked to four
division GO categories by resolved Entrez mapping. Each point is an event;
a gene may have several events, and one event may appear in several rows.
There are 207 category memberships. Counts are 77/47 for cell division,
42/26 for meiotic cell cycle, 29/17 for spindle organization and 59/36 for
chromosome segregation (events/genes). Jitter changes only vertical position.
Neither these annotations nor GO tests measure an impaired cellular process.

## FMN2 Figure 1. Independent FMN2 splicing and sequence-context analysis

**A**, Selected events by class and direction under the same criteria used
for SPIRE. **B**, All 275 selected events, row-centered PSI and visible
Euclidean/complete-linkage dendrogram. Replicates 1–2 are Fmn2 heterozygous
controls and 3–4 knockout; each biological replicate is a pool of 50 oocytes,
with technical sequencing runs merged upstream. **C,D**, Boxplots and individual selected measurements of EX
upstream/downstream intron lengths, with study-specific background and
feature-family q values. There is no gene-expression panel
because a corresponding FMN2 expression analysis was not supplied here.
The experiments have different controls, replication and tested universes.

## FMN2 Figure S1. Sensitivity across biological replicate pools

**A**, Fitted versus raw mean effects for all selected events. **B**,
Fractions retaining direction or ≥10 pp after every one-pool omission.
These descriptive checks use two pools per condition and do not refit
betAS or provide calibrated significance tests.

## FMN2 Figure S2. Reference sequence-feature distributions

**A–F**, The same six features and boxplot/point conventions as SPIRE Fig. S3, using
FMN2 populations and its own background and BH family. All observations
remain in the calculations, including MaxEnt values outside the labelled
view. Small IR groups limit the interpretation of overlapping distributions.

## FMN2 Figure S3. Absolute and background-relative GC profiles

**A–D**, The same layouts and uncertainty conventions as SPIRE Fig. S4,
computed independently within FMN2. Wider event-level ribbons reflect its
smaller selected sets. These plots do not establish the SPIRE GC association
in FMN2. Sequence-level checks and subset summaries accompany Table S7.

## FMN2 Figure S4. Regional RBP effects and correction

**A–D**, EX/IR regional effect heatmaps and regional effect distributions,
using the same complete-family exact presence tests and display rules as
SPIRE Fig. S5. No FMN2 test passes BH across all 5,278 comparisons. Trees
remain hidden; named regions are angled in C,D. Original minimum-window
asymptotic maps are reference outputs rather than this figure's test.

## FMN2 Figure S5. GO ontologies and manuscript processes

**A,B**, Top-ranked GO terms and the same nine question-driven terms,
using the rules and visual encodings of SPIRE Fig. S6. Each result is
computed within FMN2's own tested or full-annotation background. No
FMN2 primary tested-background term has BH-adjusted P ≤0.05. Identical
primary-background eligibility is used in both columns and all sensitivities.

## Convergence Figure 1. Shared targeting and directional concordance

All comparisons use 48,323 identically annotated common tested events in
4,229 genes. **A**, Study-specific and shared selected events and genes;
a shared gene may contain different selected events in each study.
**B**, Binned density of raw mean ΔPSI for all common events, with log-scaled
counts per bin. Symbols mark the 11 shared selected events by concordance
of fitted effects. **C**, Fitted effects of those 11 events; circles denote
SPIRE and triangles FMN2. Lines connect estimates, not uncertainty bounds.
Seven events are concordant and four opposing (direction-test q≈0.706
in the nine-test overlap family); global raw-change Spearman correlation is 0.0065. Table S5 includes
individual sample PSI, complete overlaps and event-count-stratified gene
overlap sensitivity. Shared targeting does not establish an identical
splicing program or a global molecular phenocopy.
