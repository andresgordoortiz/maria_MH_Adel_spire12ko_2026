# Spire1/2 knockout oocytes: RNA expression and splicing

Start with the **[publication package](results/publication/index.html)**:
two main figures, four supplementary figures, methods, results, legends,
editable manuscript text and complete supplementary tables S1–S8.
These outputs are generated locally and intentionally excluded from Git.

The computational result is selective, bidirectional RNA change, including
division-associated splicing candidates and limited shared FMN2 targets.
The data do not establish a global molecular phenocopy, a division-pathway
enrichment or a causal RBP mechanism. Microscopy and functional measurements
are needed for those parts of the manuscript.

SPIRE uses the original FASTQ groups: **Cont_1/2/7 versus DKO_4/5/6**, with
three independent animals per condition (investigator-confirmed). Historical
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
