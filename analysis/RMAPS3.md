# rMAPS3 and RBP paper figures

Current figures: `results/rbp/index.html`. Counts are already complete for both
datasets; regenerating figures does not rerun motif scanning.

```bash
bash scripts/run_R.sh --vanilla analysis/render_rmaps3.R
```

The renderer reads the checked `results/{spire,fmn2}_final/rmaps3_outputs/`
directories. It writes all RBP figures and tables into `results/rbp/`, including
a paper-panel manifest and a report with the exact selections and limitations.

| Paper panel | File prefix |
|---|---|
| 2H: combined exon heatmap | `Fig2H_EX_RBP_heatmap` |
| S6A: exon regional tracks | `FigS6A_EX_RBP_tracks` |
| S6B: exon direction boxplots | `FigS6B_EX_RBP_boxplots` |
| S6C: intron regional tracks | `FigS6C_IR_RBP_tracks` |
| S6D: combined intron heatmap | `FigS6D_IR_RBP_heatmap` |
| S6E: intron direction boxplots | `FigS6E_IR_RBP_boxplots` |

Alternative-site panels and regional/total RBP counts are repository extensions.
RBP abundance-versus-enrichment needs matching expression measurements and is
outside the agreed splicing-only scope.

## Counter ownership and manual publication

`../rMAPS3` now maintains the correction and its tests, on branch
`fix/motif-window-counting`, commit `1076240`. Publication is left to the user;
see `../rMAPS3/docs/FORK_SETUP.md`. Its README has a brief change note and
`docs/COUNTER_CORRECTION.md` records the reproducible evidence and exact scope.
The previously separate patch installer and duplicate tests are archived under
`legacy/rmaps3/development/`.

The completed analyses used exactly the same counting function as this fork.
The archived source under `legacy/rmaps3/source_used_2026-10-17/` preserves the
as-run files. It includes the former local patch installer in the module;
removing that unused installer from the maintained fork does not change counting.

The native counter misses a known motif in unchanged splice-site sequence when
900 bases are added at the opposite end, discards valid internal RI hits before
normalization, and can round a positive window width to zero. Those are
reproducible counting contradictions. An author's intentions cannot be proved
from code alone. Half-open edge conventions and suppression of unused temporary
matrix files are separately documented choices, not described as proven bugs.

## Fresh motif scans

Use the existing isolated Pixi environment. No packages are installed into the
betAS or R analysis environments. From this repository root, with prepared final
analysis inputs and **new, empty output destinations**:

```bash
pixi run --manifest-path ../rMAPS3/pixi.toml --locked python analysis/run_rmaps3.py \
  results/spire_final spire results/spire_final/rmaps3_outputs ../rMAPS3
pixi run --manifest-path ../rMAPS3/pixi.toml --locked python analysis/run_rmaps3.py \
  results/fmn2_final fmn2 results/fmn2_final/rmaps3_outputs ../rMAPS3
```

The current destinations already exist and are protected against overwrite.
Archive an earlier run explicitly before starting a replacement. Each dataset
uses three motif workers by default and retains all prepared backgrounds.

The driver checks input checksums and event manifests, converts VASTDB 1-based
closed coordinates to 0-based half-open rMAPS coordinates, and invokes `se`,
`ri`, `a5ss`, and `a3ss`. EX and IR use target/flanking geometry; alternative
sites combine the variable segment with the appropriate adjacent exon using
strand and event type. Multi-site events retain their individual VAST site-usage
direction rather than becoming binary long/short probabilities.

Settings are mm10, exon/intron extents 50/250, window 50, step 1 and one-sided
`mannwhitney_greater`. Existing foreground groups are empirical betAS FPR≤0.05
and absolute fitted ΔPSI≥0.1 after the all-sample N10 filter. The background is
the prepared complement. Explicit-coordinate mode does not reclassify events
using the CLI significance switches. Up/down describe KO event usage, not RBP
expression. The catalogue contains 114 motifs under 91 labels and omits the
paper's separately added SF3B1 motif.

For sequence length L and grid j=0,…,extent−1, start windows are `[j,j+50)`
and end windows `[L−extent−50+j,L−extent+j)`. The last end endpoint is L−1,
exclusive. RI body starts are j×L/250; very long introns are sparsely sampled.
These explicit conventions are not verified historical server coordinates.
A hit must fit fully in a tested window. Matching remains non-overlapping;
one instance can contribute to multiple overlapping windows.

## Figure methods and interpretation

The paper scripts `Figure2H_S4A_B.R` and `FigureS4C_E.R` supply the selections,
region order, clustering, tracks and score-distribution comparisons.
Spire and FMN2 replace the three drug contrasts. Existing ggplot2/patchwork tools
provide the layout without adding packages. Heatmaps select 30 exon and 15
intron labels across datasets, with Euclidean/complete row clustering and
columns ordered region → dataset → direction. The shared reversed-Cividis scale saturates at 8;
exported scores retain their full values. Colours follow the common project
palette in [../docs/COLOURS.md](../docs/COLOURS.md).

Tracks follow the deposited code's top 15 labels per dataset/direction across
regions (ties retained), with three distinct labels per region. The manuscript
caption instead says top 15 per position; this difference is explicit. Boxplots
retain motifs with raw p≤0.05 anywhere in their dataset/direction, then plot all
of their regions. Raw Wilcoxon stars compare motif-score distributions, not
biological replicates. Zero p-values remain in tables and use a display floor.
The deposited LOESS span 0.3 is singular with only five/eight unique region
positions. The adaptation uses span 0.75; all 16 dataset/direction fits and
prediction standard errors were checked finite without warnings. Curves are
descriptive, and region positions are schematic rather than binding coordinates.

Sparse, tied and imbalanced counts can make asymptotic Mann–Whitney p-values
extremely optimistic. The FMN2 EX spot check has 1/61 positive foreground events
versus 2/19,992 background events: approximate p≈1.4×10⁻²⁵ versus exact binary
allocation p≈0.0091, before multiple searching. See the per-dataset
`real_window_audit.json` and the report. [SciPy documents the approximation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html).
Raw regional/motif minima are not FDR. Shared labels can have different winning
motifs/windows; aliases, dependent events and sequence composition also matter.
Spire identities remain unresolved; FMN2 has Het controls and two pools per group.
No common causal RBP is established and no whole-catalogue calibration is claimed.

## Checks and archive

```bash
python3 -B tests/rmaps3_checks.py
pixi run --manifest-path ../rMAPS3/pixi.toml --locked python -B ../rMAPS3/tests/test_motif_counter.py
bash scripts/run_R.sh --vanilla tests/rbp_paper_checks.R
bash scripts/run_R.sh --vanilla tests/analysis_checks.R
pixi run --manifest-path ../rMAPS3/pixi.toml --locked python -B tests/rmaps3_output_checks.py spire
pixi run --manifest-path ../rMAPS3/pixi.toml --locked python -B tests/rmaps3_output_checks.py fmn2
```

Output checks cover all 741,000 window p-values per dataset, regional minima,
FASTA identities/order, map files and independently recounted real windows.
The exact binary diagnostic is a spot check, not a replacement of all tests.

Earlier Fisher baselines, incomplete attempts, figures and the as-run source are
under `legacy/rmaps3/`. Preflights and the illustrated explanation are under
`audit/rmaps3/`. The cleanup relocation CSV records SHA-256 hashes and maps old
paths to current files. Historical provenance is retained as originally written.
Comparisons with Fisher baselines change both the counter and statistical test.
