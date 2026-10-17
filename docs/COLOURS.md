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
