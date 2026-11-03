# Figure selection and narrative

The supervisor requested expression and splicing results for SPIRE,
methods/results and supplemental tables, including relevance to division.
The primary manuscript set therefore has two SPIRE figures. FMN2 is an
independent comparison with different controls, sampling and detection
opportunities, and has a separate figure set. No FMN2 expression analysis
is invented from the available splicing-only comparison inputs.

| Question | Main evidence | Why it belongs here |
|---|---|---|
| What changes in SPIRE? | SPIRE 1A–C: expression MA, splicing class counts, all selected PSI profiles | Direct expression and splicing outcomes; effect sizes and sample heterogeneity remain visible |
| Does this include division-related genes? | SPIRE 1D: all 119 selected events mapped to four named division GO categories | The email explicitly asks this; the annotation is visible and each point is an event, not a gene |
| What characterizes affected SPIRE events? | SPIRE 2A,B: EX flanking length and IR GC boxplots/measurements | Interpretable sequence contexts with supported population shifts; no claim that the knockout changes DNA sequence |
| What regulator hypotheses follow? | SPIRE 2C: within-SPIRE regional EX motif heatmap | Motif enrichment with clustered row order; final RBP-level multiplicity remains uncalibrated |
| What does FMN2 show independently? | FMN2 1: counts, PSI, EX flanks and motif landscape | Avoids implying the studies are interchangeable or that all features converge |
| What actually converges? | Convergence 1: common-universe overlap, global effects and all shared-event effects | Distinguishes shared targeting from directional agreement; the latter is weak |

Target-gene expression moved to SPIRE S1: it checks the engineered loci
and is redundant as a second main expression panel. STAR uniquely mapped
read percentages are alignment QC, neither gene expression nor VAST PSI;
they remain in Table S1 rather than taking figure space. Target-exon
junction support is explicitly labelled VAST and retained in S1.

The former ranked Cit/Gen1/Nek4 profiles and top-15 heatmap were removed.
The email did not request that ranking, and it overemphasized selected
examples. The replacement displays all matching events with named GO
categories and event classes. The full event-to-gene-to-GO mapping, event
IDs, PSI and evidence codes are provided. GO membership does not mean
that division is enriched or experimentally impaired.

The six former feature tiles were inherited from a reference analysis,
not the six most important discoveries. They are now a complete reference
feature supplement. Main-panel choice uses direct interpretability,
observed effect and distribution overlap, not smallest P alone. Matt's
strongest co-occurrence/annotation counts are retained in the all-feature
table rather than presented as biological mechanisms. Main features are
post-analysis selections; all tests and nonsignificant features remain available.

The 5′ splice-site result is supplemental. A valid negative MaxEnt score
had stretched the full-range violin, and the small FMN2 IR groups and
large distribution overlap did not justify a main mechanistic panel.
The box/point panel uses a labelled −2 to 12 view with outside-range
counts, while every measurement remains in summaries and tests.

PSI matrices retain visible row dendrograms. RBP matrices retain clustered
row order with the trees hidden, as requested. RBP labels are selected
within each study, and clustering uses row-standardized score profiles.
The supplemental motif boxes use actual exon/intron region names, angled
to fit the panels.

Matt displays were compared using three native R layouts, informed by the
oocyte project's `FiguresS2_S3_features_matt_plots.R`: boxes with event
points, violins with boxes/points, and median/50%/90% intervals. Boxplots
with points were chosen because they show the selected measurements and
quartiles directly, work with FMN2's small IR sets and require no density
smoothing. The background is summarized by its box, and group sizes and
raw medians remain explicit. Full measurements and summary tables are
retained. Local proofs are preserved in `audit/feature_go_revision_2026-11-03/feature_display_review.pdf`.

“Exploratory” has a specific statistical meaning for the RBP panels:
the strongest raw window P is selected across variants, and this selected
minimum has no final RBP-level multiple-testing calibration. The plots
therefore say “motif enrichment” with a raw-P scale, and the methods explain
the selection. This analysis supports regional motif associations;
expression, physical binding and causal RBP activity need other measurements.

GC was checked against sequence, not adjusted to match an expected result.
Matt and sequence GC agree. In SPIRE the whole-intron mean difference is
approximately 2.7–2.8 percentage points. The former plot mixed this modest
shift with FMN2's broad uncertainty and displayed absolute values on a
wide range. Separate absolute and difference-from-background profiles
show what is present. FMN2 does not show the same GC association. The
spatial background subset and whole-sequence summaries are documented.

The six SPIRE supplements cover QC, sensitivity/expression, six reference
features, GC profiles, additional motifs and GO. The five FMN2 supplements
cover sensitivity, reference features, GC, additional motifs and GO.
GO deserves a figure rather than being buried in tables: each study's
first GO panel shows the union of the top three terms per ontology and
background; its second shows nine processes motivated by the email.
These are distinct, disclosed selection rules. The tested and full
backgrounds appear side by side, with gene enrichment, counts and adjusted
P. Relevant terms such as spindle, chromosome segregation and nuclear
speck are visible. Broad terms excluded by the original 10–500 gene-size
filter are labelled “Not tested”. All tested-background terms are
nonsignificant after BH adjustment, so the complete ontology comparison
belongs in the supplement while the main division panel presents direct
event annotations. The full hierarchy can produce related parent/child
terms; no unreported semantic filtering is applied. Protein-impact
annotations and complete nonsignificant results remain in tables. Original outputs are
preserved locally. The limitations in the results text must accompany
any claim about phenocopy, pathway enrichment or causal RBPs.
