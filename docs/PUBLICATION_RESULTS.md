# Results: SPIRE expression and RNA processing

## SPIRE loss affects expression and several classes of splicing events

RNA sequencing of three independent biological replicates per condition
identified 111 differentially expressed genes at BH-adjusted P ≤0.05
(40 increased and 71 decreased; SPIRE Fig. 1A; Table S2). Fifty also had
an absolute shrunken log2 fold change ≥1. Spire2 RNA was depleted, whereas
some Spire1 exons retained RNA support (SPIRE Fig. S1C,D). The saved
DESeq2 model included condition only. All samples were retained; PCA
revealed substantial between-sample heterogeneity (SPIRE Fig. S1A,B).

Among 287,393 events with complete PSI and corrected support ≥10 in every
sample, 1,547 events in 1,000 named genes met empirical per-event FPR ≤0.05
and absolute fitted ΔPSI ≥10 percentage points. These comprised 551 EX,
271 IR, 330 alternative 5′ and 395 alternative 3′ events (SPIRE Fig. 1B).
Changes were bidirectional: 784 increases and 763 decreases. In particular,
139 IR events increased and 132 decreased; there was no general
predominance of increased intron retention. These selections do not have
5% multiple-testing FDR control. Clustered event-level PSI profiles show
the sample heterogeneity within this selected set (SPIRE Fig. 1C).

## Affected genes include relevant processes without preferential GO targeting

The corrected mapping links 127 selected events to 79 named genes in four
division categories, excluding the engineered loci (SPIRE Fig. S7D; Table S4).
This adds eight events/four names to the direct-symbol-only annotation through
unique alias rescue. In the full non-target selection, 971/999 symbols map,
including 41 rescued aliases; unresolved or ambiguous names remain explicit.
A gene can contribute several events and belong to several GO categories.
The union has 79 selected genes versus 76.36 expected from its 780
testable annotated genes, so pooling the categories does not establish enrichment.

SPIRE Fig. 1D places seven manuscript-relevant GO categories in context,
showing observed/expected gene counts and complete-family q. Cell division has
47 selected genes versus 47.38 expected from the tested background;
chromosome segregation has 36 versus 34.95; nuclear speck has 32 versus
32.19. Thus these annotations occur among the changes, but their numbers
are not greater than expected for the assay's accessible genes.

The full-annotation background predicts only 22.67, 15.32 and 15.77 genes,
respectively, giving apparently strong enrichment. The difference reflects
background composition, not removal of the observed events. Using identical
term eligibility and BH across BP/CC/MF, no primary tested-background term
passes in either study (SPIRE Fig. S6; FMN2 Fig. S5). Matching selected genes
by testable-event count increases the SPIRE expectations to 52.79, 40.21 and
34.59; none of the 18 targeted matched tests passes BH (SPIRE Fig. S7).
Stricter selections yield a few synapse/localization terms in SPIRE, with
no preferential enrichment of the nine manuscript themes. These are
sensitivity results, not grounds to change the primary event threshold.
Expression GO using finite-padj DESeq2 genes as background also yields no
BH-significant term (Table S6). Functional annotation supports integration
with the cellular phenotype; these analyses do not establish impaired division.

## SPIRE events have distinctive sequence contexts and regional motif patterns

Selected EX events had longer upstream introns than the unselected
background: medians were 2,118 nt for increased-inclusion events and
2,108 nt for decreased-inclusion events, versus 1,444 nt in the background
(SPIRE Fig. 2A; Table S7; exact statistics in the source tables).
The main sequence panels emphasize interpretable genomic measurements,
not the smallest P values across all Matt features. Annotation-dependent
transcript co-occurrence metrics, correlated features and all nonsignificant
comparisons remain in the complete feature table.

Selected IR events had higher whole-intron GC content. Mean GC was 48.01%
for increased-retention events and 48.10% for decreased-retention events,
versus 45.27% in the Matt background; medians were 47.22%, 46.52% and
44.47%, respectively (SPIRE Fig. 2B). Existing feature-family BH-adjusted
P values were 0.0020 and 0.0047 for increased and decreased retention.
Recalculation from each intron sequence reproduced Matt GC to numerical
precision. The spatial-profile background is a 1,500-event subset, whose
mean whole-intron GC was 45.41%; this modest sampling difference does not
explain away the GC association. Whole-intron distributions and positional
profiles answer different questions. Study-specific difference-from-background
profiles make the approximately 2–3 percentage-point SPIRE differences
visible without claiming uniform shifts at every sequence position
(SPIRE Fig. S4; GC audit in Table S7).

An exact regional presence analysis identifies two SPIRE motif associations
after BH across all 5,278 tests (SPIRE Fig. 2C; Table S8). PCBP2-labelled
motifs occur in 48/280 increased-use EX events in the upstream-exon 3′
search region, versus 7,654/89,623 background events (17.14% versus 8.54%;
q=0.0165). PPRC1-labelled motifs occur in 16/83 increased-use Alt5 events
in the intron after the long exon, versus 1,530/27,506 background events
(19.28% versus 5.56%; q=0.0337). Each event is counted once, with catalogue
variants combined and the regional search fixed before testing.

These are regional sequence associations, with residual event dependence
and sequence-composition confounding possible. They do not establish RBP
binding, expression or causal activity. Clustered effects for all displayed
EX/IR labels and full-region summaries are in SPIRE Fig. S5. The original
minimum-window asymptotic scores remain reference data because sparse
counts can give optimistic P values; they are not the primary motif test.

## FMN2 is presented as an independent comparison dataset

FMN2 has two biological replicate pools per condition and heterozygous
controls. Of 50,415 tested events, 275 events in 208 genes met the same
selection criteria. Its event counts and clustered PSI profiles are shown
separately (FMN2 Fig. 1A,B). Increased-inclusion EX events had longer
upstream and downstream introns than the FMN2 background (medians 3,311
versus 1,509 nt and 2,321 versus 1,416 nt; BH-adjusted P=0.00102 and
0.00219; FMN2 Fig. 1C,D). The regional RBP effects
are shown separately in FMN2 Fig. S4; no exact regional motif association
passes the complete BH family.

The FMN2 IR sample sizes are smaller (29 decreased and 14 increased events
with Matt measurements). Intron GC was not significantly shifted. Predicted
5′ splice-site strength was lower in decreased-retention events (median
7.64 versus 8.84; BH-adjusted P=0.0242), but overlapping distributions and
the small selected set argue against making this a principal mechanistic
claim. It is retained alongside all six reference features in FMN2 Fig.
S2; the MaxEnt box/point panel shows −2 to 12 with explicit counts outside the view.
Negative MaxEnt scores are valid model scores. All observations remain
in calculations and source tables. Across the six features, boxplots and
individual selected-event points make the population spread and small
FMN2 IR sets visible without density smoothing.

FMN2 Fig. S5 presents its GO results with the same mapping, engineered-locus
exclusion and fixed-family rules. No primary BP/CC/MF term passes pooled
BH with its tested-gene background. Whole-annotation and sensitivity
results remain visible as reference comparisons.

## Convergence is limited to selected targets and sequence contexts

Within 48,323 common tested events in 4,229 genes, SPIRE selected 98 events
and FMN2 selected 263. Eleven events and 22 genes were shared (Convergence
Fig. 1A). Gene overlap exceeded a uniform expectation (3.44; one-sided
Fisher q≈5.60×10⁻¹²) and an event-count-stratified permutation expectation
(4.77; Monte Carlo q≈0.00015; Table S5; BH across all nine overlap/direction tests). Seven shared events had concordant
fitted effects and four had opposing effects. Across all common events,
raw changes had Spearman correlation 0.0065 (Convergence Fig. 1B,C).

Thus, shared exon-flanking architecture and selected target overlap coexist
with different overall responses. These results do not establish a global
molecular phenocopy. Actin, nuclear positioning, nuclear-envelope forces,
chromatin, speckles and division outcomes require independent imaging or
functional evidence; they were not measured by this RNA analysis.
