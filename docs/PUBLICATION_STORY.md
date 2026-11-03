# Publication story and panel decisions

The primary question is the RNA phenotype of SPIRE loss. Two main SPIRE
figures cover expression/splicing, functional context, sequence features
and the motif associations supported by the revised tests. FMN2 remains
an independent comparison set. The compact convergence figure assesses
shared targets without implying a globally concordant molecular phenocopy.

| Question | Main evidence | Decision |
|---|---|---|
| What changes? | SPIRE 1A–C: DESeq2 MA, splicing counts, all selected PSI profiles | Direct expression and event-level results, with sample heterogeneity visible |
| Which processes do the affected genes include? | SPIRE 1D: seven question-driven GO categories, observed versus expected genes and full-family q | Functional context; no preferential pathway enrichment is claimed |
| What characterizes affected sequences? | SPIRE 2A,B: EX upstream-intron length and IR GC | Interpretable distributions, all selected measurements and complete-family Matt q |
| Which motif associations survive a calibrated sparse-count test? | SPIRE 2C: PCBP2/EX and PPRC1/Alt5 regional presence | Exact tests with BH across all 5,278 RBP/region/direction comparisons |
| What does FMN2 show independently? | FMN2 1: counts, PSI and EX flanks | Independent controls, replication and detection; no BH-significant regional motif is promoted |
| What converges? | Convergence 1: common-universe overlap and effect agreement | Limited shared targeting, weak overall directional agreement |

GO is now in the main figure, with a clear primary background: testable
ontology-annotated genes. The seven terms represent the supervisor's themes;
they are an editorial choice, not a preregistered set or the seven strongest
P values. Their q values use the full BP+CC+MF family. Top-ranked terms,
both backgrounds and all nine question-driven terms are in the GO
supplements. Full annotation is a reference sensitivity, not an equally
valid replacement chosen for significance.

The GO audit resolves the apparent paradox: cell division has 47 affected
genes versus 47.4 expected from the tested background. More category
memberships do not create overrepresentation. Direct-symbol mapping also
missed old names: unique alias rescue adds eight division events/four
names, giving 127 events in 79 genes. Individual event annotations move
to SPIRE S7D; all matching events remain visible without favored-gene ranking.
S7 also exposes gene/event opportunity and selection sensitivities. The
full inference/background review is in PUBLICATION_INFERENCE_REVIEW.md.

The former minimum-window RBP scores had sparse-count calibration problems.
The main motif evidence therefore uses event-level regional presence,
unioning catalogue variants and the fixed native search before exact Fisher
testing. Two SPIRE associations survive the complete BH family; none do
in FMN2. Broader clustered EX/IR effect heatmaps and regional summaries are
supplemental. They retain row clustering with hidden trees and use actual,
angled region names. Cell color encodes effect; corrected significance is
shown separately. Original scans remain reference data. No RBP expression,
binding or causal activity is inferred.

Matt layouts were compared with the local oocyte project: boxes with
individual measurements, violins with boxes/points, and quantile intervals.
Boxplots with selected-event points were chosen for clear medians, spread,
unequal population sizes and FMN2's small IR sets. No density smoothing or
value clamping is used. All six reference features and all other Matt tests
remain available; main feature choice is post-analysis and never defines a
smaller correction family. MaxEnt is supplemental with labelled display
limits and outside-view counts. GC agrees independently with sequence.

Target-locus RNA measurements, PCA and sensitivity remain supplemental;
STAR mapping QC is table-only. Expression GO is also a table: no term passes
its complete BH family. BH is used for all nominal inferential families;
counts, PCA and heatmaps are descriptive. The betAS simulation-tail field
retains its original effect/FPR selection interpretation, with no calibrated
event-level FDR claim. The inference ledger makes these distinctions explicit.

SPIRE has seven supplements: QC, splicing/expression sensitivity, all six
reference features, GC profiles, regional motif effects, full GO and the
background/event-opportunity audit. FMN2 has five: sensitivity, features,
GC, regional motifs and GO. Tables S1–S8 contain all results, including
nonsignificant ones, mapping candidates, null draws and original references.
No microscopy, protein measurements or functional division assays were
available to establish the cellular phenotype in this RNA analysis.
