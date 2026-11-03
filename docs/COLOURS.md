# Figure colours

`analysis/palettes.R` is the single palette definition for current report plots.
There is no universally best palette: these choices separate experimental
identity from signed effects and keep numeric magnitudes ordered in lightness.

| Meaning | Colour |
|---|---|
| Spire dataset | Purple `#CC79A7` |
| FMN2 dataset | Green `#009E73` |
| Increased inclusion/site use/retention, positive effect | Vermillion `#D55E00` |
| Decreased inclusion/site use/retention, negative effect | Blue `#0072B2` |
| Unselected events | Grey `#B3B3B3` |
| Control / heterozygous control | Grey `#777777` |
| Knockout | Charcoal `#222222` |
| Shared concordant / opposing targets | Charcoal / amber |

Categorical colours follow Okabe–Ito. Direction never changes between datasets
or event types. Centred PSI, z-scores and signed feature shifts use the same
blue–off-white–vermillion diverging scale. Numeric magnitudes use reversed
Cividis: pale yellow for low values, dark blue for high values. GO adjusted-p
legends reverse this numeric order so stronger evidence is also darker.
Protein-impact categories have fixed named colours, including categories absent
from an individual panel. Labels, position, facets and legends accompany colour.

Cividis was designed for perceptual uniformity and colour-vision accessibility:
[Nuñez et al., 2018](https://doi.org/10.1371/journal.pone.0199239).
The category palette comes from [Okabe and Ito's guidance](https://jfly.uni-koeln.de/color/).
These choices do not guarantee every reproduction or vision condition; greyscale
still needs the explicit labels for categorical colours.

This changes figure appearance only. Thresholds, event selection, clustering,
scores and tests are unchanged. Raw rMAPS3 maps are preserved as tool outputs
with their original rendering; their derived paper panels use this palette.
Previous report figures and checksums are archived under
`legacy/figures_before_palette_2026-10-17/` and `audit/palette_2026-10-17/`.

The publication assembly has a separate, user-requested Lena/RBP style.
It uses the installed Helvetica/Arial-compatible Liberation Sans face, 8-point base
text, 7-point ticks/legends, 9-point titles and 12-point panel letters.
Direction is blue `#2166AC` / red `#B2182B`, background grey `#BDBDBD`;
study is purple `#6A3D9A` / green `#33A02C`. Controls are grey and SPIRE DKO
purple. The candidate PSI heatmap uses sequential blue; regional RBP scores
use reversed Cividis. Labels, facet order and annotation tracks redundantly
identify categories. Borders and horizontal grids follow the Lena reference.

Style sources read locally: `lena_quantseq/R/go_analysis_publication_plot.R`
and the q7/q8 plotting scripts at commit
`7a84859bae1a1399606b5776a2dbed5a81fdf2ec`, plus this repository's
`analysis/rbp_paper_figures.R`. The publication code has no runtime dependency
on the Lena repository. Theme and compositions live in
`analysis/publication_helpers.R` and `analysis/publication_panels.R`.
Earlier exploratory reports retain their existing palette and layouts.

The study-separated revision clusters all selected PSI rows with visible
dendrograms. RBP labels are selected within study and clustered using
row-standardized regional profiles; cells retain raw scores. The shared
plotting code emits distinct SPIRE/FMN2 sets and a compact convergence
figure. Full rationale is in PUBLICATION_STORY.md.

The subsequent feature/GO revision uses boxplots with individual selected
event measurements in the same palette. RBP trees are hidden while row
clustering is retained; PSI dendrograms remain visible. GO facets name the
three ontologies and both backgrounds; grey-to-purple encodes the original
BH-adjusted P and point size gives contributing genes. The oocyte project's
feature scripts were read as visual references only; no source was executed
or modified there and no runtime dependency was added.

The current main EX and supplemental IR RBP heatmaps restore native rMAPS3
raw P values. Reversed Cividis encodes −log10 minimum raw P on 0–8, with
values ≥8 saturated and explicitly labelled. Direction strips retain blue/red;
row clustering remains and trees are hidden. Exact regional-presence tests
and BH results stay in Table S8, separately from heatmap values.
The main GO dot plot uses gene ratio on x, point area for gene count and
purple-to-grey for complete-family BH q on a log scale from 0.05 to 1.
Its subtitle explicitly states that none of the primary terms passes BH.
Three named ontology facets organize the top-three-per-ontology selection.
