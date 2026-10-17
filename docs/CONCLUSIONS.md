# Spire and FMN2: what the splicing data support

Reviewed 17 October 2026 using the final betAS tables and corrected rMAPS3 runs.
The strongest common finding is **overlap in affected targets and sequence
architecture**, with limited agreement in direction. These data do not establish
a common causal RBP or direct spliceosome regulation. Spire genotype labels
remain unresolved; its conclusions, including overlaps with FMN2, are provisional.
No samples have been relabelled.

## What changes in each dataset

Selection is N≥10 and complete PSI in every profile, empirical betAS FPR≤0.05
and absolute fitted ΔPSI≥10 percentage points (pp). This FPR is not a
multiple-testing-adjusted FDR. Pdiff and raw sample separation provide additional
descriptive support without changing the paper's selection rule.

| Readout | Spire1/2 DKO versus labelled Cont | Fmn2 KO versus Het |
|---|---:|---:|
| Profiles per group | 3 | 2 merged pools |
| Tested events | 289,583 | 50,415 |
| Selected events | 1,436 (0.50%) | 275 (0.55%) |
| Affected named genes | 945 | 208 |
| EX / IR / Alt5 / Alt3 | 537 / 241 / 293 / 365 | 145 / 46 / 25 / 59 |
| Increased / decreased event use | 711 / 725 | 141 / 134 |
| Increased / decreased intron retention | 127 / 114 | 15 / 31 |
| Selected in the 48,908-event common universe | 94 (0.19%) | 264 (0.54%) |

**Spire:** the labelled contrast contains bidirectional changes across event
types, without generalized increased intron retention. Its larger absolute count
accompanies a much larger tested universe and cannot establish a stronger
phenotype than FMN2. Even matched-universe rates remain descriptive because
controls, coverage and replication differ.

Spire1's retained-intron event MmuINT0151614 is a target-locus/QC candidate,
not evidence of autoregulation: fitted ΔPSI is +70.1 pp, arithmetic mean ΔPSI
is +23.4 pp, and samples are heterogeneous. Exon support at both Spire loci
separates Cont_a/b and DKO_c from Cont_c and DKO_a/b. Independent PCA and held-out
marker checks do not validate the original identities or a c-label swap.
Spire2 has no event passing all-profile N10; its differential splicing is
untestable. See the [identity report](../results/spire_identity_check/REPORT.md)
and [target-locus review](../results/biological_review/REPORT.md).

**FMN2:** changes are also selective and bidirectional. Selected introns more
often decrease retention (31 versus 15), arguing against generalized
intron-retention failure. Camk2g exon MmuEX0008978 (−53.7 pp) and Fmn1 exon
MmuEX0019373 (−48.6 pp) have complete raw-PSI separation of the two pools per
group; these are candidates, without demonstrated protein consequences.
Lower Fmn2 exon support in both KO pools is consistent with the manipulation,
but no Fmn2 event survives all-profile N10. Controls are heterozygous, not wild type.

## What is common

The comparison uses **48,908 events in 4,229 named genes tested in both studies**.

- **25 affected genes overlap**, versus 3.26 expected under a uniform-gene null
  (one-sided Fisher p≈1.14×10⁻¹⁶). Shared targeting exceeds that null, which does
  not correct gene length or event-count bias.
- **18 exact event IDs overlap at 15 loci:** 11 concordant, seven opposing.
  Direction agreement is not convincingly above 50% (exploratory sign-test
  p=0.481); raw mean ΔPSI across all common events is almost uncorrelated
  (Spearman 0.0096). The knockouts are not interchangeable splicing phenocopies.
- **Rev3l alternative 3′-site choice** is the clearest shared computational
  candidate after inspecting complete raw separation and Pdiff≥0.95 in both.
  Its two complementary rows represent one locus. **Zfp957** exon exclusion
  also separates both groups (−11.2/−15.3 pp; FMN2 Pdiff=0.927). Ktn1 and Asap1
  are concordant secondary candidates with overlapping FMN2 pools. Nexn and
  Stag1 retained introns change in opposing directions.
- **Longer flanking introns accompany increased exon inclusion in both.**
  Upstream medians: 1,867 versus 1,402 nt in Spire; 3,311 versus 1,509 nt in FMN2,
  comparing up events with each study's ndiff background. Downstream medians:
  2,284 versus 1,339 nt and 2,321 versus 1,416 nt. All four comparisons pass the
  exported BH threshold. These are pre-existing genomic features of the affected
  exons, not introns becoming longer in KO. Ten feature comparisons have matching
  median-shift direction and BH significance; several features are correlated.
- **No common enriched GO pathway is established.** No term passes in both
  separate study sets on the common-tested background; shared-gene GO has no
  significant BP/CC/MF terms. Full-dataset BP themes with the paper's OrgDb
  background (112 terms in Spire, seven in FMN2) disappear with each study's
  tested-gene background. Those background-sensitive themes should not headline
  the interpretation.

See the [comparison report and candidate tables](../results/comparison/REPORT.md).

## What the RBP figures add

Corrected scans restore motif counting and reproduce the requested paper panels,
but do **not** identify a defensible common driver. Same-label/region/direction
raw-p≤0.05 selects 85/91 labels for IR-up and 86/91 for IR-down in both studies:
the criterion is too broad to establish specificity.

Sparse counts can produce extreme asymptotic Mann–Whitney p-values. One FMN2
exon window has one positive event among 61 versus two among 19,992 background
events: native p≈1.4×10⁻²⁵, versus exact conditional p≈0.0091 before multiple
searching. Shared SRSF1 scores use different motif variants/windows: 21/127
positive events in the Spire winner versus 1/15 in FMN2. Motif presence is not
binding, expression, protein abundance or causal activity. See the
[RBP report](../results/rbp/REPORT.md) and support-count audits.

## Biological interpretation and next priorities

A connection through actin-dependent cell organization is plausible: Spire and
FMN2 cooperate in mouse-oocyte actin assembly
([Pfender et al., 2011](https://pubmed.ncbi.nlm.nih.gov/21620703/)). Independent
oocyte work connects cytoplasmic forces with nuclear-condensate remodeling and
RNA processing ([Al Jord et al., 2022](https://www.nature.com/articles/s41467-022-32675-5)).
**Our interpretation is a hypothesis:** perturbing this system could affect
overlapping splicing targets, with different directions and additional targets
in each dataset. These tables cannot distinguish direct regulation from altered
cell state, mechanics or technical effects.

First recover Spire genotype/sample provenance or allele-aware alignment evidence.
Prioritize Rev3l and Zfp957 for shared splice-junction validation, and Camk2g/Fmn1
for FMN2 follow-up. RBP prioritization needs calibrated sparse-count statistics
and multiple-search control. No simulation, expression analysis or label
correction was performed for this review. The shared [palette](COLOURS.md)
changes appearance only; biological cutoffs and statistical results are unchanged.
