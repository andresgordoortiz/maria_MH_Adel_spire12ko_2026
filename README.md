# Spire1/2 knockout oocytes: RNA expression and splicing

Start with the **[publication package](results/publication/index.html)**:
two SPIRE main figures and six supplements, a separate FMN2 main figure
and five supplements, and a compact convergence figure. The package includes
methods, results, legends, editable manuscript text and supplementary tables
S1–S8. Outputs are local and intentionally excluded from Git.

The [figure-selection rationale](docs/PUBLICATION_STORY.md) explains the
main/supplement boundary. SPIRE expression/splicing is the primary story;
all 119 selected events in division-annotated genes are shown without a
favored-gene ranking. Sequence context and exploratory motifs provide a
second layer. FMN2 is presented independently because study design and
responses differ. Shared targets do not establish global molecular phenocopy.

SPIRE uses the original FASTQ groups: **Cont_1/2/7 versus DKO_4/5/6**, with
three independent biological replicates per condition (investigator-confirmed). Historical
condition labels are superseded. FMN2 has two biological pools per condition
and heterozygous controls. betAS's `FDR` column is an empirical per-event
false-positive rate, not a multiple-testing-adjusted FDR.

## Reproduce the publication assembly

With the completed input analyses and the existing isolated R environment:

```bash
bash scripts/run_local_R.sh tests/publication_checks.R
bash scripts/run_local_R.sh analysis/publication.R results/publication
python3 analysis/package_publication.py results/publication
python3 tests/publication_output_checks.py results/publication
```

Choose a new output directory for a new render. The wrapper runs R in a
network-disabled container with source/references/libraries mounted read-only;
it installs nothing. It reuses the local rig R installation and is not a
fresh-machine installer. The two renv lockfiles pin the upstream analysis
environments. Setup/bootstrap scripts must be run only in an approved isolated
environment. See [full reproduction instructions](analysis/README.md).

## Repository map

| Location | Purpose |
|---|---|
| `analysis/`, `scripts/`, `config/`, `tests/` | Source code, execution, settings and checks |
| [Publication audit](docs/PUBLICATION_AUDIT.md) | What works, limitations and release checks |
| [Methods](docs/PUBLICATION_METHODS.md), [legends](docs/PUBLICATION_LEGENDS.md) | Manuscript source text |
| `results/publication/` | Final composed figures, tables, source data and manuscript exports |
| `results/index.html` | Extended analysis figure gallery |
| `reports/` | Completed upstream analyses and exploratory reports |
| `data/inputs/spire_fastq/` | Imported RNA-seq, DESeq2, VAST and betAS outputs |
| `data/processed/`, `data/reference/` | Intermediates and reference resources |
| `legacy/`, `audit/` | Preserved earlier work and audit records |

Git contains source, configuration, tests, lockfiles and documentation only.
Inputs, figures, tables, environments and raw reads are ignored. A fresh clone
therefore needs the analysis data package; a public SPIRE data accession has
not yet been supplied. Original analyses and investigator changes are preserved.

[FASTQ workflow](docs/FASTQ_CLUSTER.md) · [Cluster/betAS](docs/CLUSTER.md) ·
[Source attribution](analysis/SOURCES.md) · [Palette](docs/COLOURS.md) ·
[Detailed prior conclusions](docs/CONCLUSIONS.md)
