# Figure legends

## Figure 1. Selective expression and splicing changes in Spire1/2 DKO oocytes

**A**, Normalized Spire1 and Spire2 gene counts, plotted after addition of
one on a logarithmic axis. Each point is one independent animal (three per
condition). Total Spire1 counts can include retained exons and should be
interpreted alongside exon support in Fig. S1D. **B**, Imported DESeq2
expression results: shrunken log2 fold change versus −log10 BH-adjusted P.
Blue and orange identify decreased and increased genes at adjusted
P ≤0.05; grey indicates unselected genes. Genes without finite adjusted
P are absent from the plot and retained in Table S2. Labels identify
Spire1/2 and the three lowest-adjusted-P genes. **C**, Selected EX, IR,
Alt5 and Alt3 events by direction. Selection requires all-sample N10,
empirical betAS per-event FPR ≤0.05 and |fitted ΔPSI| ≥10 percentage
points; this is not multiple-testing FDR control. **D**, All selected
SPIRE events, row-median centered and ordered by event class, fitted
effect and event ID. Samples retain their manifest order. Colors saturate
at ±50 percentage points; no rows are removed. **E**, Individual raw PSI
and group means (horizontal bars) for three division-associated examples.
Candidates were defined using four explicit GO categories and ranked by
minimum signed leave-one-animal raw effect, with one event per gene and
the engineered Spire loci excluded. Cit is an alternative 3′-site-use
event; Gen1 and Nek4 are EX events. These examples do not establish pathway
enrichment or functional impairment. Tables S2–S4 contain all tests and
candidates, not only the examples shown.

## Figure 2. Limited shared splicing targets in SPIRE and FMN2 knockouts

The comparison is restricted to 48,323 events with matching annotations,
covering 4,229 named genes. SPIRE compares three independent animals per
group; FMN2 compares two biological pools per group, with heterozygous
controls. **A**, Study-specific and shared selected events and genes
within this common universe. A gene can be selected in both experiments
through different events. **B**, Binned density of raw mean ΔPSI in both
experiments for every common tested event; grey intensity is the number
of events per bin on a log scale. Symbols identify the 11 shared selected
events, classified by concordance of their *fitted* effects. Raw and fitted
effect directions can differ. **C**, Fitted effects of all shared events;
circles denote SPIRE and triangles FMN2. Lines connect the two estimates
for the same event and are not uncertainty intervals. The top four rows
have opposing fitted directions and the remaining seven are concordant.
The exploratory two-sided binomial sign test gives P=0.549. **D**, Raw
PSI for the same 11 events and every animal/pool, on a common 0–100% scale.
None of the seven concordant events has strict complete group separation
in both experiments. Shared targeting does not establish an identical
splicing program, genetic interaction or phenotypic equivalence. Table S5
includes the common universe, all candidates, overlap tests and a
permutation sensitivity analysis stratified by tested event count per gene.

## Figure S1. Animal-level quality control and targeted-locus support

**A**, PCA of the 500 most variable VST expression genes, excluding
Spire1/2. **B**, PCA of current SPIRE raw PSI, retaining one highest-variance
event per named gene with variance ≥25 pp² and excluding Spire1/2. Both
PCAs center variables without variance scaling; axis percentages indicate
explained variance. Point labels retain the original FASTQ identifiers.
**C**, Percentage of uniquely aligned reads from STAR, with a zero-based
axis. **D**, Corrected VAST inclusion-plus-exclusion junction support for
observed Spire1/2 exon events before all-sample N10 filtering, displayed
as log10(support+1). Events are not independent gene-expression counts.
Residual Spire1 support and loss of other target exons can coexist. These
RNA measurements do not replace DNA genotyping. All animals were retained;
neither PCA nor read-alignment rate alone identifies a biological batch.

## Figure S2. Sensitivity of selected splicing events and tested universes

**A**, Fitted versus raw mean PSI changes for every selected SPIRE event;
the diagonal indicates identity. **B**, Within each event class, the
percentage retaining the fitted sign after every one-animal omission and
the percentage retaining a ≥10 pp raw change in that direction after
every omission. These are descriptive checks without refitting betAS.
**C**, Selected event counts across the indicated empirical-FPR and
absolute fitted-effect thresholds. A lower empirical-FPR threshold does
not convert this measure into multiple-testing FDR control. **D**, All
tested and common-tested event counts, showing the unequal detection
opportunities before restricting the comparison universe. Raw data,
selection thresholds and sensitivity summaries are in Tables S3 and S5.

## Figure S3. Expression profiles and functional annotation

**A**, VST expression of all 111 genes selected at BH-adjusted P ≤0.05,
centered on each gene's six-sample mean and ordered by shrunken log2FC.
Color saturation is ±3 VST units. **B**, Expression-only, shared and
splicing-only selected genes among named genes with both evaluable
splicing and finite expression adjusted P. Gene symbols are deduplicated;
this is a descriptive intersection. **C**, Numbers of enriched GO BP,
CC and MF terms for all selected SPIRE genes, comparing the full OrgDb
background with the tested-N10 background. Existing enrichGO cutoffs
include BH P ≤0.05 and q ≤0.2. **D**, Four focused one-sided Fisher tests
of division-category overrepresentation, using all tested named SPIRE
genes as background and BH correction across the four tests. The dashed
line denotes adjusted P=0.05. Category membership includes inherited
GOALL annotations. Candidate division genes are not evidence of an
enriched division pathway. Tables S4 and S6 provide membership, evidence
codes and complete enrichment results.

## Figure S4. Exploratory sequence features, transcript impact and RBP motifs

**A**, Ratios of median genomic flanking-intron length for selected EX
events relative to each study's own unselected background. Open/filled
symbols distinguish BH P >0.05/≤0.05 in the existing two-sided
Mann–Whitney feature-test family. These are descriptive feature ratios,
not fold changes in intron length after knockout; no confidence intervals
are implied. Table S7 gives group sizes, medians and every feature test.
**B**, VASTDB v3 predicted consequence categories for selected SPIRE
events, split by effect direction. Inclusion/exclusion annotations must
be combined with the observed direction; the category alone is not a
prediction of protein loss in DKO. **C**, RBP labels with raw P ≤0.05 in
the same broad region and event direction in both studies. Bar labels
are counts out of 91 tested RBP labels; the axis shows percentages.
Different motif variants or windows can supply a label's evidence in
each study. Raw p-values do not correct the search across motifs and
windows, and sparse-count approximation can overstate evidence.
The widespread overlap does not identify a specific causal regulator.
Complete annotations and exploratory motif summaries are in Tables S7–S8.
