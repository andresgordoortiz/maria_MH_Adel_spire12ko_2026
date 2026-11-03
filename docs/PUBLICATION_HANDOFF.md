# SPIRE RNA analysis: final figure handoff

Prepared for Marie-Hélène and Maria, 3 November 2026.

Use **SPIRE Figure 1 and SPIRE Figure 2** as the two RNA-analysis figures.
Their numbering is local to this contribution; Maria can renumber them when
combining them with the imaging figures. No additional GO significance is
needed to describe the expression and splicing results. The conclusions
below are conditional on the existing analyses and their stated limitations.

## The two figures

- **SPIRE Figure 1: expression and splicing.** DESeq2 expression changes,
  counts by splicing class, all selected event-level PSI profiles, and every
  selected event in the four division annotation categories. This directly
  answers the request for expression and division-related splicing results.
- **SPIRE Figure 2: functional and sequence context.** The tested-background
  GO result, exon-flanking intron lengths, retained-intron GC, and native
  rMAPS3 motif profiles. GO is explicitly nonsignificant. Motif heatmap
  colors are raw selected P values, not corrected discoveries or evidence
  of an RBP mechanism. Exact regional-presence tests are separate in S8.

S1–S3 provide quality/stability checks and additional sequence context.
FMN2 and convergence figures are separate supporting comparisons; the
SPIRE conclusion does not depend on a global molecular phenocopy claim.

## What the data support

Among three independent biological replicates per condition, DESeq2
identified **111 genes at BH-adjusted P ≤0.05** (40 increased, 71 decreased).
The splicing analysis selected **1,547 events in 1,000 named genes** using
empirical per-event FPR ≤0.05 and |fitted ΔPSI| ≥10 percentage points.
These events include **127 events in 79 genes annotated to division
processes**, excluding the engineered loci. The betAS selection does not
provide 5% multiple-testing FDR control; call these selected candidate
events, and retain the raw-PSI and sample-omission checks.

Among the 127 division-annotated events, 116 retain the fitted direction
after every one-sample omission, 51 retain a raw change of at least 10
percentage points after every omission, and 43 strictly separate every
control and DKO PSI value. These are descriptive stability checks, not
independent validation or a new significance threshold (S4).

Suggested manuscript conclusion:

> Spire1/2 DKO oocytes exhibited gene-expression differences and selected
> alternative-splicing changes, including events in genes annotated to
> cell-division processes. These processes were not preferentially enriched
> relative to the genes testable in this experiment. The RNA results
> identify candidates for integration with the cellular phenotype, without
> establishing impaired division or a global FMN2-like molecular response.

## Why the GO result is not a failed figure

The two questions are different: whether division-annotated genes contain
selected events, and whether those genes occur more often than expected.
For cell division, **47 genes are observed versus 47.38 expected** with the
tested-gene background. For chromosome segregation, the counts are
36 versus 34.95; for nuclear speck, 32 versus 32.19. No primary GO term
passes BH correction. The counts support annotation, not enrichment.

Keep the tested-gene background and complete correction family. A
whole-annotation background includes genes that could not enter the
selection and produces a different expectation. The
[Gene Ontology reference-list guidance](https://geneontology.org/docs/go-enrichment-analysis/)
recommends using the genes from which the selected list could have arisen.
Full results, background comparisons and nonsignificant tests remain in S6.

## Files to share

- `SPIRE_handoff.zip`: both SPIRE figures as embedded-font PDFs and
  450-dpi PNGs, manuscript text/Word/HTML when Pandoc is available, this
  handoff, and the table guide. The manuscript also documents the separate
  supporting FMN2 analysis.
- `Figures_Supplementary.zip`: the three shared supporting figures.
- `Supplementary_tables.zip`: complete S1–S8 CSV tables, plotted source
  data, field definitions, provenance and the R session. This archive is
  about 37 MB; share it by a file link if email attachment limits require it.
- `index.html`: the local gallery, including the independent FMN2 and
  convergence figures mentioned in the full analytical manuscript.

Results, methods and legends are in `Manuscript.docx` and `text/`; table
definitions are in `text/TABLES.md`. S2 contains the expression results,
S3 the complete SPIRE splicing analysis, and S4 every division-annotated
candidate. Keep event IDs as text when opening CSV files in a spreadsheet.

## Suggested email

Dear Marie-Hélène,

I have assembled the current Control versus Spire1/2 DKO RNA analyses into
two figures, with figure legends, methods, results and complete supplementary
tables. Figure 1 summarizes expression and splicing changes, including all
selected events in genes annotated to division processes. Figure 2 provides
the functional and sequence context, with supporting quality and sensitivity
checks in the supplementary material.

The expression analysis identifies 111 genes at BH-adjusted P ≤0.05. The
splicing analysis selects 1,547 candidate events, including 127 events in
79 division-annotated genes. GO does not show significant overrepresentation
against the testable-gene background, so the text distinguishes affected
division-related genes from pathway enrichment. The splicing selection and
its statistical limitations are stated in the methods and legends.

These figures can be integrated with Maria's cellular phenotype results.
The independent FMN2 comparison is available separately; it does not
establish a global molecular phenocopy.

Best wishes,
Andrès

## Remaining author information

The analytical figures are assembled; submission still needs the allele
and genotyping details, donor ages, oocyte collection/stage, RNA isolation
and library preparation, sequencer details and SPIRE data accession. No
imaging, protein or functional validation was supplied for this analysis.
Alignment, DESeq2 fitting and betAS fitting were not repeated for this handoff.
