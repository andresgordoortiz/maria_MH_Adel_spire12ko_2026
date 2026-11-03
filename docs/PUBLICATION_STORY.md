# Publication story and panel decisions

The primary question is the RNA phenotype of SPIRE loss. Two SPIRE main
figures show expression/splicing, division-annotated events, GO context,
sequence features and regional RBP motif patterns. FMN2 remains an independent
main figure; a small convergence figure assesses shared targets and effect
agreement. There are exactly three shared supplementary figures in total.

| Question | Main evidence | Decision |
|---|---|---|
| What changes? | SPIRE 1A–C: DESeq2 MA, splicing classes, all selected PSI profiles | Direct expression and event-level results with sample heterogeneity visible |
| Are division-related genes affected? | SPIRE 1D: all 127 splicing events in 79 division-annotated genes | Retain every matching event; no favored-gene ranking or claim of pathway enrichment |
| What does GO rank highest? | SPIRE 2A: standard gene-ratio/count/BH-q dot plot, three top terms per ontology | Transparent ranking without forcing division terms; none passes primary BH |
| What characterizes affected sequences? | SPIRE 2B,C: EX upstream intron length and IR GC | Distributions and individual measurements, with complete-family Matt q |
| What regional motif patterns appear? | SPIRE 2D: clustered EX RBP effects | Effects and exact-test BH dots shown separately; hidden dendrogram |
| What does FMN2 show independently? | FMN2 1: counts, PSI, EX flanks and RBP heatmap | Independent controls, replication and detection; no significant motif association promoted |
| What converges? | Convergence 1: common-universe overlap and effect agreement | Limited shared targeting with weak global directional agreement |

The restored division-event panel answers the supervisor's question directly:
it displays splicing events in annotated host genes. GO answers a separate
question about overrepresentation. Cell division has 47 selected genes versus
47.38 expected from testable genes. The union of four division categories has
79 versus 76.36 expected. Presence and preferential targeting are different
claims. No primary GO term passes BH in either study. The main GO plot ranks
all terms without requiring a connection to the manuscript; related parent/
child terms remain visible. Its q values retain the full BP+CC+MF family.

The motif heatmaps restore the regional overview, with row clustering and
hidden trees. Colors encode log2 presence odds ratios, and dots indicate exact
regional presence tests corrected across all 5,278 tests per study. Catalogue
variants and fixed regional windows are unioned before counting each event
once. PCBP2/EX and PPRC1/Alt5 pass in SPIRE, none in FMN2; the Alt5 result
remains in the complete table rather than receiving a separate main panel.
These associations measure sequence context, not RBP expression, binding
or causal activity. Original asymptotic minimum-window scans remain reference
outputs; they do not define publication significance.

The three supplements are selected for distinct supporting questions:

1. **S1, quality and stability:** SPIRE expression/PSI PCA, clearly labelled
   VAST target-exon support, clustered expression, and both studies'
   leave-one-sample checks. This addresses heterogeneity and reproducibility.
2. **S2, complementary features:** exon number and IR splice-site strength
   in each study. These add information beyond the main intron-length/GC
   panels and retain nonsignificant comparisons and FMN2's small groups.
3. **S3, intron sequence context:** both IR RBP matrices and background-relative
   IR GC profiles. This complements the main EX matrices and tests the visual
   interpretation of the whole-intron GC association without duplicating it.

Additional fitted/raw plots, threshold grids, full GO background diagnostics,
all features, absolute/EX GC profiles and uninformative motif-effect distribution
boxes do not receive extra supplementary figures. Their underlying results
remain in Tables S1–S8, source data and the inference audit. This keeps the
figure set compact without deleting nonsignificant results or narrowing any
correction family to displayed panels.

Matt boxplots with individual selected-event points were chosen after comparing
box/point, violin and quantile displays with the local oocyte project. They
show medians, spread and unequal group sizes without density smoothing. MaxEnt
views report all outside-range values; sequence-level GC checks reproduce Matt.
All feature choices are post-analysis editorial decisions, not a discovery
ranking or a smaller testing family.

BH applies to all nominal inferential families. Counts, PCA and heatmaps are
descriptive. betAS retains its empirical per-event FPR/effect selection,
without a calibrated event-level FDR claim. The background/inference audit
explains alias rescue, testable genes and event-count matching. Imaging or
functional experiments are still required for actin, nuclear positioning,
chromatin, speckles or impaired division; RNA alone does not establish these
cellular phenotypes or a global FMN2 molecular phenocopy.
