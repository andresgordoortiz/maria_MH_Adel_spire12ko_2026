# Inference, backgrounds and the manuscript decision

This audit changes the statistical presentation as well as the figures.
The main conclusions use a tested-gene background and complete, disclosed
BH families. Whole-annotation results are reference comparisons. The
publication is conditional on the existing betAS event selection; it does
not claim calibrated event-level splicing FDR.

## Why many annotated genes do not imply GO enrichment

The former figure displayed 196 category memberships for 119 unique events
in 75 genes. Categories overlapped and individual genes contributed multiple
events. Those counts describe annotation; enrichment asks whether the
fraction exceeds the fraction already available in the assay.

The mapping audit recovered 41 selected SPIRE symbols through unique
aliases in the pinned OrgDb, improving mapping from 930 to 971 of 999
non-engineered selected symbols. Ambiguous aliases remain unresolved.
For example H1foo maps to H1f8 and March2 to Marchf2. This was applied to
the background as well as selected genes, before interpreting enrichment.
The corrected division union is 127 events in 79 named genes. Its larger
annotation coverage does not establish enrichment. Pooling the four
categories correctly gives 79 unique selected genes, versus 76.36 expected
among the 780 testable BP-annotated genes in their union. Repeated category
memberships must not be summed as independent genes.

| Term | Selected genes observed | Expected: tested background | Expected: full annotation | Expected: event-count matched |
|---|---:|---:|---:|---:|
| Cell division | 47 | 47.38 | 22.67 | 52.79 |
| Chromosome segregation | 36 | 34.95 | 15.32 | 40.21 |
| Nuclear speck | 32 | 32.19 | 15.77 | 34.59 |

For cell division specifically, the BP query contains 900 annotated
selected Entrez genes and the tested background contains 9,193 annotated
genes, of which 484 have the term. Expected count is
900 × 484 / 9,193 = 47.38. The observed fraction is 47/900 = 5.22%,
versus 5.26% in the tested background. Its unadjusted overrepresentation
P is about 0.55: removing BH would not make this a significant enrichment.
The full annotation has only 648/25,728 = 2.52%, giving 22.67 expected.
These denominators exclude unresolved/ontology-unannotated genes and
engineered loci; they are not the original 1,000 named-gene count.

The three rows explain the apparent contradiction: affected genes
are present, but their numbers are close to or below the tested-background
expectations. The whole-genome comparison predicts far fewer because it
includes genes that could never enter the splicing selection. That can
recapitulate properties of the measured oocyte transcriptome and event
catalogue rather than a selective response to SPIRE loss.

Selection opportunity also matters. Approximately 2.7% of genes with one
testable event are selected, compared with 15.3% of genes with 64–127
events. The sensitivity test samples genes within log2 event-count strata,
preserving the number selected in each stratum. All 18 targeted tests
(nine terms in each study) remain nonsignificant after their complete BH
correction. Matching controls event-count opportunity, with residual
coverage, expression, gene-length and annotation biases still possible.

## Decisions for GO

- Resolve symbols directly, then rescue only unique aliases; export candidates,
  unresolved names and mapping routes. Deduplicate Entrez IDs.
- Exclude engineered Spire loci from functional interpretation in both studies.
- Use all testable, ontology-annotated genes as the primary background.
- Fix term eligibility at 10–500 primary-background genes. Reuse exactly those
  terms in the full-annotation comparison and all selection sensitivities.
  This removes the previous background-dependent size-filter exclusions.
- Apply BH once across BP + CC + MF for each study/background/selection
  analysis. Main question-driven terms retain q from the entire family;
  no favorable correction restricted to the displayed terms is used.
- Keep the original event selection as primary. Report leave-one-sample
  ≥10 pp and strict separation as sensitivities; do not promote the threshold
  producing the most attractive GO result.

No primary tested-background GO term passes BH in SPIRE or FMN2. Some
synapse/localization terms pass only in stricter SPIRE subsets; they do
not establish the requested division/speckle phenotype. A separate GO
analysis of BH-selected expression genes also has no BH-significant term.
It uses the finite-padj DESeq2 universe rather than a splicing universe.

The compact main GO panel therefore shows selected processes, gene counts,
tested-background expectations, with full-family q values. It is
functional context, not a claim of enriched pathways. Full ontologies,
background comparisons, alternative selections, matched nulls and the
individual division events are supplemental.

## One correction policy, with explicit test families

BH applies to nominal inferential tests across the complete family before
selecting rows to display. Families are separate for distinct outcomes:
DESeq2 genes; expression GO; splicing GO; Matt features; regional motif
presence; and cross-study overlap. Every family and size is exported in
Table S1. Raw P is retained for reproducibility; significance uses q ≤0.05.
Descriptive heatmaps, effect sizes, sample PCA, event counts and GC profiles
do not become hypothesis tests and do not receive decorative P values.

For GO, all three ontologies are pooled. Matt retains its complete feature
family, including all three contrasts, rather than correcting just the
six displayed features. Cross-study overlap now corrects all nine reported
tests together, including the direction and matched-gene sensitivities.

The betAS column named FDR is a simulation-tail estimate, as verified in
the installed source and original method. In SPIRE, 209,923 events have
constant raw PSI and zero reported FPR; none survive the effect-size
selection. Thus its numerical zeros must not be treated as generic,
calibrated P values and passed through BH to manufacture an event FDR
claim. This does not show that all selected changes are false. It means
that the current event list is an effect/FPR selection, supported by raw
PSI and omission checks. Calibration or another appropriate model is needed
before calling it an event-level FDR discovery set. Exact two-sided label
permutation has minimum P=0.10 with three versus three samples and 1/3 with
two versus two; low replication limits assumption-free confirmation.

## Why the motif analysis changed

The old regional scores selected the smallest asymptotic Mann–Whitney
P across windows and motif variants. One FMN2 EX example reported about
1.4×10⁻²⁵ from one motif-positive foreground event versus two background
events; its exact binary allocation P is about 0.0091. Correcting the
uncalibrated minimum would not repair the sparse-count approximation.

The revised primary motif test counts presence anywhere in the fixed native
regional search. Catalogue variants for one RBP are unioned, and each event
contributes once. The test is an exact one-sided Fisher/hypergeometric test
of motif-positive event frequencies versus the existing unselected sequence
background. BH includes all 91 RBPs, 29 regions across four event types,
and both directions: 5,278 tests per study. There is no P-based window or
variant selection. Coordinate geometry is checked against explicit windows.

Two SPIRE associations pass: PCBP2 in the EX upstream-exon 3′ region
(48/280 versus 7,654/89,623; q≈0.0165), and PPRC1 in the Alt5 intron after
the long exon (16/83 versus 1,530/27,506; q≈0.0337). None pass in FMN2.
These are regional sequence associations under the stated event-level
model. Correlated events from a gene, GC/length composition and catalogue
specificity remain limitations; BH does not fix them. The results do not
measure protein abundance, binding or causality. The two supported
associations enter the main figure. Clustered effect-size heatmaps and
regional distributions move to supplements, with corrected significance
indicated separately. Original window scans remain reference data.

## Sources and interpretation

The [Gene Ontology guidance](https://geneontology.org/docs/go-enrichment-analysis/)
recommends a reference containing genes that could enter the selected list.
[Young et al., 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2872874/)
demonstrate RNA-seq selection bias; event-count matching here is a
study-specific sensitivity analysis, not a claim to reproduce goseq.
The [R p.adjust documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html)
describes BH's application to defined test families and dependence assumptions.
The [betAS paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10946425/)
defines FPR as a simulated null-tail proportion, distinct from
multiple-testing FDR. Input hashes, original-test reproduction, complete
new tests and exact motif counts accompany this audit.
