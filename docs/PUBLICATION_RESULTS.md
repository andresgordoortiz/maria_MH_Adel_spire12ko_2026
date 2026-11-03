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

## Events in division-annotated genes address the manuscript's functional question

The manuscript's division question was operationalized using four named GO
categories: cell division, meiotic cell cycle, spindle organization and
chromosome segregation. Their membership identified 119 selected splicing
events in 75 genes after excluding the engineered Spire loci (SPIRE Fig.
1D; Table S4). Each point represents an exon, retained-intron or alternative
splice-site event; GO annotation belongs to the host gene. An event can
appear in more than one category. All matching events are shown, with
no ranking of favored genes or selection of illustrative extreme profiles.

This is functional annotation, not evidence of pathway enrichment or
impaired division. The four focused overrepresentation tests were
nonsignificant after BH correction (all adjusted P≈0.597), and broad GO
enrichment did not persist with the tested-gene background (Tables S4,S6).
SPIRE Fig. S6 now shows ranked BP, CC and MF terms under both backgrounds,
and a separate set of nine terms motivated by the manuscript. With the
full-annotation background, chromosome segregation (BH-adjusted P=0.00160),
mitotic spindle (2.63×10⁻⁷) and nuclear speck (0.00978) were among the
significant terms. Their tested-background adjusted P values were approximately
1.00, 0.476 and 1.00. Thus the relevant annotations are present, while
robust overrepresentation relative to genes testable in this assay is not
established. The strongest tested-background BP terms were axon development
and axonogenesis (both adjusted P=0.0643); these overlapping annotations do
not demonstrate neuronal phenotypes in oocytes.
The event list provides directly traceable hypotheses for integration with
cellular measurements, without establishing isoform function.

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

EX RBP motif profiles provide exploratory regional hypotheses (SPIRE
Fig. 2C). Rows are selected within SPIRE and clustered by standardized
regional score patterns; the trees are hidden. The IR
landscape and motif-score distributions are supplemental (SPIRE Fig. S5).
Scores are minima over motif variants and windows, not calibrated RBP-level
p-values, expression measurements or binding assays. Sparse-hit asymptotic
scores and correlated motifs limit regulator specificity. No causal RBP
is established.

## FMN2 is presented as an independent comparison dataset

FMN2 has two biological replicate pools per condition and heterozygous
controls. Of 50,415 tested events, 275 events in 208 genes met the same
selection criteria. Its event counts and clustered PSI profiles are shown
separately (FMN2 Fig. 1A,B). Increased-inclusion EX events had longer
upstream and downstream introns than the FMN2 background (medians 3,311
versus 1,509 nt and 2,321 versus 1,416 nt; BH-adjusted P=0.00102 and
0.00219; FMN2 Fig. 1C,D). The EX motif landscape is ranked and clustered
within FMN2 (FMN2 Fig. 1E).

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

FMN2 Fig. S5 presents its GO results separately. Full-background terms
included nuclear speck (BH-adjusted P=0.0171) and actin binding (0.00845);
tested-background values were 0.999 and 0.403. No BP, CC or MF term retained
BH-adjusted P ≤0.05 with the FMN2 tested-gene background.

## Convergence is limited to selected targets and sequence contexts

Within 48,323 common tested events in 4,229 genes, SPIRE selected 98 events
and FMN2 selected 263. Eleven events and 22 genes were shared (Convergence
Fig. 1A). Gene overlap exceeded a uniform expectation (3.44; one-sided
Fisher P=6.22×10⁻¹³) and an event-count-stratified permutation expectation
(4.77; Monte Carlo P≈0.00010; Table S5). Seven shared events had concordant
fitted effects and four had opposing effects. Across all common events,
raw changes had Spearman correlation 0.0065 (Convergence Fig. 1B,C).

Thus, shared exon-flanking architecture and selected target overlap coexist
with different overall responses. These results do not establish a global
molecular phenocopy. Actin, nuclear positioning, nuclear-envelope forces,
chromatin, speckles and division outcomes require independent imaging or
functional evidence; they were not measured by this RNA analysis.
