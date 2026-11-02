# SPIRE FASTQ to VAST and DGE on CRG

Local presentation: downloaded pipeline outputs are kept in `data/inputs/spire_fastq/`;
`results/` contains the current figures, and `reports/` contains full analysis reports.
The CRG output paths in the commands below are unchanged.

Run from `/users/aaljord/agordo/git/maria_MH_Adel_spire12ko_2026` after copying
these repository changes to CRG. The FASTQs stay in `data/raw/renater_files/`.

The tracked manifest is [config/spire_fastq.tsv](../config/spire_fastq.tsv):
controls **1, 2, 7** and double KO **4, 5, 6**, using the supplied filenames.
`filename_date_label` preserves `12mars` / `14mars`; it is not a verified batch
or collection date. The investigator confirmed one independent animal per file on 2 November 2026.
Historical sample correspondence is now documented in `docs/CONCLUSIONS.md`;
the historical condition labels are superseded. Samples are not merged or paired.

**Input layout: single-end**, with strandedness inferred by nf-core/rnaseq
(`auto`). Complete scans of controls 1 and 2 found only R1-labelled reads and
no interleaved mates (see the checks below). The other four files are assumed
to use the same layout. This workflow assumes ordinary full-length RNA-seq
without UMIs, not a 3′ tag-counting assay. If additional mate files are found,
adapt the sample-sheet preparation before using them.

## References and submission

[config/fastq_cluster.env](../config/fastq_cluster.env) reuses the existing
CRG VAST workflow at `/users/aaljord/agordo/git/vast-tools_nextflow`, VASTDB at
`/users/aaljord/agordo/VASTDB`, and scratch under
`/nfs/scratch01/aaljord/agordo/maria_MH_Adel_spire12ko_2026`.

VAST uses **mm10 / Mm2**. RNA-seq and DGE now default to the completed
run's **GRCm38 / Ensembl 88** references under
`$PROJECT_DIR/reference/mm10_ensembl88/`. Earlier defaults used GRCm39;
the imported run parameters verify that the actual analysis used GRCm38.
Both RNA-seq and DGE receive the same `GENES_GTF`. Reference paths can
still be overridden; the submission script checks their existence.

```bash
cd /users/aaljord/agordo/git/maria_MH_Adel_spire12ko_2026

# Optional: select another existing, matched FASTA/GTF pair.
# export GENOME_FASTA=/absolute/path/genome.fa
# export GENES_GTF=/absolute/path/genes.gtf

bash scripts/submit_fastq.sh
squeue --me
```

This checks all six raw files, the reference paths and `VASTDB/Mm2`, generates
`data/processed/spire_fastq/samples/rnaseq.csv`, and queues three driver jobs:

1. **rnaseq**: nf-core/rnaseq **3.27.0**, STAR + Salmon, fastp (minimum length
   35), QC and saved trimmed reads.
2. **vast**, dependent on RNA-seq success: the existing VAST Nextflow workflow,
   using the saved `fastp/*.fastp.fastq.gz` files. It skips repeated trimming
   and FASTQ QC, runs align/combine/compare, and skips MATT and betAS. The
   existing betAS workflow remains separate.
3. **dge**, also dependent on RNA-seq success: nf-core/differentialabundance
   **2.0.0**, `rnaseq` profile (DESeq2), with Salmon gene counts **and gene
   length offsets**. The manifest supplies sample groups; the contrast is
   `DKO_vs_Cont`, so positive log2 fold change means higher expression in DKO.

The saved DESeqDataSet verifies `~0 + condition` with a DKO-minus-Cont contrast,
without a batch term. Animal independence is investigator-confirmed; filename
date labels alone do not establish batches. Check the QC/PCA before interpretation.
DESeq2's adjusted p-values use BH correction; the pipeline's report defaults
select adjusted p-value ≤0.05 and fold change ≥2 (in either direction).
Complete differential tables are also retained.

Each driver uses Nextflow **25.10.4** and launches Slurm tasks on `genoa64`,
with the existing CRG time/QOS mapping and Singularity containers. It loads
the existing Nextflow, Java and Python modules; it does not install tools in
the project R environments or host Python libraries. Source/raw/reference
binds are read-only; task outputs go to scratch. Scratch must be outside the
project's read-only mount. Nextflow downloads pipeline code/plugins and
containers into its caches as needed, so network access is required on the
first run. VAST retains the shared workflow's `latest` container tag; preserve
the cached image and shared workflow revision for reproducibility.

If `TOWER_ACCESS_TOKEN` is already exported, drivers add `-with-tower`.
No token is required for ordinary Slurm execution.

## Outputs and resuming

All new results are isolated under `results/spire_fastq/`:

- `rnaseq/`: QC, alignments, saved `fastp/` reads and `star_salmon/` matrices.
- `vast/inclusion_tables/spire_fastq_INCLUSION_LEVELS_FULL-mm10.tab`:
  the new inclusion table, using the manifest's `SPIRE_*` identifiers.
- `vast/compare_results/`: VAST group comparisons.
- `dge/`: DESeq2 tables, exploratory plots and the differential-abundance report.

Driver logs are `logs/spire_{rnaseq,vast,dge}.JOBID.log`. Submitted job IDs
are appended to `runs/spire_fastq/slurm_jobs.tsv`. Nextflow sessions are kept
separately in `runs/spire_fastq/{rnaseq,vast,dge}` with work in scratch.

After fixing a failed stage, submit only that stage from the repository root:

```bash
bash scripts/submit_fastq.sh rnaseq
# Once RNA-seq has completed successfully:
bash scripts/submit_fastq.sh vast
bash scripts/submit_fastq.sh dge
```

Drivers use `-resume`; keep their work and run directories. Do not run two
copies of the same stage simultaneously. `afterok` dependents remain pending
if RNA-seq fails: cancel those old pending IDs before submitting replacements.
The existing processed tables, historical results, `metadata/samples.csv`
and betAS configuration are unchanged. The new inclusion table needs an
explicit betAS sample mapping before use; do not substitute it for the old
table with the old `Cont_a/b/c` mapping.

## Local checks

```bash
python3 -I tests/fastq_checks.py
bash -n scripts/submit_fastq.sh scripts/fastq_pipeline.sbatch config/fastq_cluster.env
```

The check uses synthetic FASTQs and mocked Slurm/Nextflow executables to
verify file validation, groups, single-end sheets, `afterok` dependencies,
paths containing spaces, stage commands and length-offset hand-off. It does
not execute alignment, DESeq2 or live Slurm jobs.

The three configurations were also resolved with Nextflow 25.10.4 inside an
isolated container, checking the Slurm executor, queue, read-only binds and
DESeq2 profile.

On 2026-10-29, both downloaded control FASTQs were scanned completely using
Python's standard library, without installing software. Every record had a
FASTQ header/separator and matching, nonempty sequence and quality lengths.
Both complete gzip streams passed their CRC checks.

| FASTQ | Reads | Read length | Mate labels | Adjacent identical read IDs |
|---|---:|---:|---|---:|
| `reads_1-controle-12mars.fastq.gz` | 94,973,298 | 118 nt throughout | R1 only | 0 |
| `reads_2-controle-12mars.fastq.gz` | 73,098,440 | 118 nt throughout | R1 only | 0 |

Neither file contains interleaved R1/R2 mates. In the first 10,000 records
of each, the sample indexes were respectively `TCCGCCAATT+CAGCACGGAG` and
`AGGACAGGCC+AGCTGGAATG`, with no shared read IDs in those subsets. These
observations support treating the files as separate, single-end inputs;
they do not establish biological independence, genotype or whether R2 files
were originally produced elsewhere. The other four FASTQs were not inspected.

Counts, compressed-file SHA-256 fingerprints and deletion confirmations are
kept locally in `audit/spire_fastq_checks/control_{1,2}.json`, with the header
subset check in `header_subset.json`. The two downloaded local FASTQs were
deleted after validation at the user's request, freeing 10,899,774,045 bytes
(10.90 GB). Audit artifacts remain excluded from Git.

Upstream documentation: [rnaseq 3.27.0](https://nf-co.re/rnaseq/3.27.0/docs/usage),
[differentialabundance](https://nf-co.re/differentialabundance/latest/docs/usage/),
and the [2.0.0 RNA-seq profile](https://github.com/nf-core/differentialabundance/blob/2.0.0/conf/profile/rnaseq/rnaseq.config).
