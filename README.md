# Spire and FMN2 oocyte splicing

Splicing-only analyses of Spire1/2 DKO versus labelled controls, and Fmn2 KO versus
heterozygous controls. Input data and results are local and are excluded from Git.

Start with **`results/index.html`**. RBP paper figures are together in
**`results/rbp/index.html`**.

[Biological conclusions](docs/CONCLUSIONS.md) compare the knockouts and shared
targets. [Figure colours](docs/COLOURS.md) explains the consistent palettes.

| Location | Contents |
|---|---|
| `final_input_tables/` | Final inclusion and betAS input tables |
| `results/spire_final/`, `results/fmn2_final/` | Current splicing results and raw motif maps |
| `results/comparison/` | Cross-study splicing and GO comparisons |
| `results/rbp/` | Figure 2H, S6A–E, RBP extensions, tables and methods |
| `results/spire_identity_check/`, `results/biological_review/` | Sample-identity checks, PCA and target-locus review |
| `analysis/` | Reproducible analysis and plotting code |
| `scripts/`, `config/` | CRG setup, whole-table betAS jobs and dataset configuration |
| `audit/` | Validation evidence and relocation manifests |
| `legacy/` | Previous inputs, runs, figures and superseded attempts |
| `data/`, `reference/`, `metadata/` | Local processed inputs, references and manifests |

[Cluster installation and submission](docs/CLUSTER.md) ·
[Analysis workflow](analysis/README.md) ·
[RBP workflow and fork](analysis/RMAPS3.md) ·
[Sources and attribution](analysis/SOURCES.md)

Recreate RBP figures from the completed motif scans:

```bash
bash scripts/run_R.sh --vanilla analysis/render_rmaps3.R
```

Diagnostic code is kept in `analysis/`, outside generated results:

```bash
python3 analysis/check_spire_raw.py
bash scripts/run_R.sh --vanilla analysis/check_spire_labels.R
bash scripts/run_R.sh --vanilla analysis/review_plots.R
```

These diagnostics require the existing local input/review tables and do not
assign genotypes or change sample labels.

The corrected rMAPS3 source is in the sibling `../rMAPS3` directory, branch
`fix/motif-window-counting`, commit `1076240`. Its `docs/FORK_SETUP.md` contains
the manual GitHub fork/push commands. This repository does not carry another
maintained copy of that counter.

Spire sample identities remain unresolved. FMN2 controls are heterozygous and
have two pools per group. Motif p-values are exploratory: sparse-count
approximations and searches across many windows can exaggerate significance.
