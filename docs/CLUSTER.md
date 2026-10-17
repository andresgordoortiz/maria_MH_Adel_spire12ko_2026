# Spire1/2 DKO and FMN2 KO: CRG betAS workflow

The repository separates the CRG simulation workflow (`scripts/`) from the Spire and FMN2 splicing figures and comparison (`analysis/`). See [the analysis guide](../analysis/README.md) for reproducible commands and [panel mapping](../analysis/panels.tsv). Locally, unchanged inclusion tables and completed simulation exports live in `final_input_tables/`; provisional results and historical tables are archived under `legacy/`. New galleries are generated under `results/spire_final/`, `results/fmn2_final/` and `results/comparison/`. Input data, legacy exports, audit results, generated sample sheets, filtered tables, logs and results are not committed.

## Input and environment

The default input is the existing CRG file:

```
/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/notebooks/inclusion_tables/spiredko_INCLUSION_LEVELS_FULL-mm10.tab
```

No input needs to be copied into this repository. To use another location, export `SPIRE_INCLUSION_TABLE=/absolute/path/to/table.tab` before preparation and keep the same setting when submitting jobs.

R 4.5.2 is managed by rig in user mode. renv restores the pinned dependencies, including betAS 1.2.1 at commit `cd37737c75ab17a64652baf2fdb00d1b7c67c294`, matching the [reference repository](https://github.com/andresgordoortiz/oocyte_splicing_ssa_pladb_tub_2025). betAS and renv are the direct dependencies; betAS requires its plotting/Shiny dependency closure. Bioconductor 3.22 metadata is also pinned. Restore the lockfile on CRG; do not copy compiled libraries from another machine.

## Prepare on CRG

From this repository's root, run:

```bash
export PATH="$HOME/.local/bin:$PATH"
srun --job-name=spire-setup --qos=short --time=02:00:00 \
  --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=8G \
  bash scripts/prepare_cluster.sh
```

If rig is missing, setup first installs the official rig 0.10.0 Linux binary into `~/.local/bin`, verifying its pinned SHA-256. No sudo is needed. An existing rig on PATH is reused. This follows the [official user-install method](https://rig.r-lib.org/install.html).

Setup then installs/restores R and packages, validates the configured sample/quality column pairs, generates `metadata/samples.csv` from the group mapping in `config/analysis.R`, and filters the existing inclusion table. It creates all required output directories. The allocated node needs internet access for installation and shared access to the project and home directory. The input is read without modification.

The filter is `betAS::filterEvents(..., N=10)`: no missing PSI and corrected inclusion+exclusion coverage >=10 in **every sample** (six for Spire, four for FMN2). An independent coverage check must agree exactly. ANN annotation-module exons are retained and exported separately. The filtered FULL table, paired PSI/Qual RDS are generated in `data/processed/`; QC reports are generated in `metadata/`.

The configured contrast is Cont_a/b/c (A) versus DKO_a/b/c (B). Positive delta means higher inclusion/site use/retention in DKO. The labels do not establish independent animals, pooling, pairing or batches.

If using an existing CRG R 4.5.2 instead of rig, activate that installation and set `RSCRIPT` to its executable. Run `"$RSCRIPT" scripts/setup.R` followed by `scripts/run_R.sh scripts/01_filter.R`. For a non-portable R build, `RENV_CONFIG_REPOS_OVERRIDE=https://cloud.r-project.org` selects source packages during restore; system build dependencies may be required.

## Select the FMN2 dataset

The default dataset remains Spire. Set `BETAS_DATASET=fmn2` for the second inclusion table. Its default input is:

```
/users/aaljord/aaljord/Al Jord Lab - CRG/Lab Members/Andres Gordo Ortiz/Data/Oocyte/FMN2KO and SPIRE/FMN2KO SPIRE sent to MHV 14 02 2025/oocyte_fmn2_spire_report_andresgordo/raw_data/fmndko_INCLUSION_LEVELS_FULL-mm10.tab
```

The spaces are handled inside R; no shell quoting of the default path is necessary. For another location, export `FMN2_INCLUSION_TABLE='/absolute/path with spaces/table.tab'` and keep the same override during submission.

```bash
# The R environment is already installed: filter the FMN2 table in an allocation.
BETAS_DATASET=fmn2 srun --job-name=fmn2-filter --qos=short --time=02:00:00 \
  --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=8G \
  bash scripts/run_R.sh scripts/01_filter.R

# After filtering succeeds, submit one 48-hour simulation/export job.
BETAS_DATASET=fmn2 scripts/submit_crg.sh results/fmn2_N10_nsim1000_run01
```

For setup plus filtering on a fresh installation, run `BETAS_DATASET=fmn2` with the earlier srun command for `scripts/prepare_cluster.sh` instead. Submission explicitly exports the selected dataset to Slurm. To resume FMN2, use `BETAS_DATASET=fmn2 sbatch --job-name=fmn2-betAS scripts/betas.sbatch results/fmn2_N10_nsim1000_run01`.

This input has **four merged profiles**, corresponding to **two Fmn2+/− control pools (A, Het) and two Fmn2−/− pools (B, KO)**. Positive delta means greater inclusion/site use/retention in KO. [GEO GSE103718](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103718) describes 50 oocytes per pool and three technical sequencing replicates per biological sample. The merged-run names and the [reference run manifest](https://github.com/andresgordoortiz/oocyte_splicing_ssa_pladb_tub_2025/blob/main/notebooks/metadata/PRJNA406820_fmn2dko_metadata.tsv) give this mapping, cross-checked against all 12 GEO sample records:

| Merged runs | Original pool | Group |
|---|---|---|
| SRR6026682–SRR6026684 | Fmn 7 | Het, Fmn2+/− |
| SRR6026685–SRR6026687 | Fmn 8 | Het, Fmn2+/− |
| SRR6026688–SRR6026690 | Fmn 9 | KO, Fmn2−/− |
| SRR6026691–SRR6026693 | Fmn 10 | KO, Fmn2−/− |

The 12 technical runs are not treated as 12 biological replicates. The names establish the expected grouping; verifying the original FASTQ merging commands remains outside this audit.

All 661,162 input IDs are unique mouse events, with correctly paired PSI/count columns and nonmissing PSI in range 0–100. N10 in **all four samples**, plus complete PSI, retains **50,415** events: 15,625 S/C1/C2/C3/MIC exons, 16,816 Alt3/Alt5 events, 13,457 IR events and 4,517 ANN exons. Controls individually have 62,908–68,872 events with coverage >=10 versus 106,330–111,183 in KO; this imbalance needs consideration alongside the small biological sample size.

FMN2 outputs are isolated in `data/processed/fmn2/`, `metadata/fmn2/` and the chosen run directory. CSV filenames use `fmn2_` / `differential_fmn2_`; job logs use `logs/fmn2-betAS.JOB_ID.out` and `.err`. Default Spire paths, groups and filenames remain compatible with existing filtered objects and runs. Data and audit artifacts remain excluded from Git.

## Submit one long simulation job

After preparation succeeds, from the login/submission node:

```bash
scripts/submit_crg.sh results/spire_N10_nsim1000_run02
```

This invokes **sbatch once**, requesting one CPU, 64 GB and `vlong`/48 hours. The same job runs the full-table betAS calculations and exports the CSVs. It prints a job ID and monitoring command. Existing run directories are refused; use a new name for each analysis. Before replacing a previously submitted array run, cancel its array and collection job IDs recorded in that run's `slurm_jobs.txt` using `scancel JOB_ID...`.

The command inside the allocation is:

```bash
scripts/run_R.sh scripts/02_simulate.R results/spire_N10_nsim1000_run02
scripts/run_R.sh scripts/03_collect.R results/spire_N10_nsim1000_run02
```

The simulation script calls these upstream functions, each on the **entire** paired, N10-filtered PSI/Qual tables:

```r
fdr <- do.call(betAS::prepareTableVolcanoFDR,
               c(common, list(nsim = config$nsim)))
pdiff <- do.call(betAS::prepareTableVolcano, common)
fstat <- do.call(betAS::prepareTableVolcanoFstat, common)
```

`common` supplies `psitable=psi`, `qualtable=qual`, `npoints=500`, the three Cont columns as A and the three DKO columns as B, `maxDevTable=maxDevSimulationN100`, `seed=TRUE`, and `CoverageWeight=FALSE` (plus labels/colors). `nsim=1000`. The fixed betAS internal seed is 21122023. All coverage-passing events are simulated. The old `chunk_size` configuration field is retained for compatibility with already filtered objects; submission ignores it and creates no array. Older `chunks.csv` files may remain and are unused. Filtering does not need to be repeated for this scheduling change.

The full-table calls retain simulated distributions across events, which is why memory is increased from 8 GB per old task to 64 GB for this job. The job requests 48 hours; full-dataset completion within that limit has not been benchmarked on CRG. The CPU count is one because these upstream wrappers run sequentially. The runtime loads the same GLPK and FlexiBLAS modules used for the working CRG environment whenever the module command is available.

Completed full-table stages are saved atomically as `FDR.rds`, `Pdiff.rds` and `Fstat.rds`; the joined table is `statistics.rds`. If a job times out, rerun the same initialized run without calling the initialization helper again:

```bash
sbatch scripts/betas.sbatch results/spire_N10_nsim1000_run02
```

The script reuses completed stages only after checking run identity and event order. An interrupted stage restarts from its beginning. Collection verifies the full event set, configuration and input-object/lockfile hashes. A `COMPLETE` marker appears only after successful CSV export. Do not alter configuration or processed inputs during a run. Jobs use the installed packages and never install dependencies.

## Results and statistical interpretation

Under `results/RUN/csv/`:

- `spire_all_events.csv`: all passing events, individual PSI and new statistics.
- `spire_statistics_{exons,introns,alt,ann}.csv`: complete family-specific results.
- `differential_spire_*.csv`: upstream FDR <=0.05 and |delta PSI| >=0.10.
- `spire_*_combination.csv`: the same limits plus Pdiff >=0.95.

These are explicit new selection defaults, not a reconstruction of the old unknown thresholds. Full results retain all unselected events. The simulation export stops at event statistics; `analysis/run.R` creates heatmaps, protein-impact summaries, GO and the other splicing figures.

Pdiff, Fstat and upstream FDR are joined by EVENT after checking that their delta PSI agrees. Delta PSI is on a 0–1 scale; 0.10 means 10 percentage points. The betAS variable named `FDR` is a per-event simulated null tail fraction, **not a BH-adjusted value**. Its source computes `sum(abs(simDeltaPsi) >= abs(foundPsi)) / nsim` separately for each event, without using the total number of rows. The formula therefore has no global multiple-testing denominator. However, this revision calls `simulate_reads()` before resetting its RNG inside `estimateFDR()`, so splitting a table can change simulated values through RNG state; exact split/full equality is not guaranteed. The workflow passes the full table; its empirical resolution with 1,000 simulations is 0.001. See the [pinned source](https://github.com/DiseaseTranscriptomicsLab/betAS/blob/cd37737c75ab17a64652baf2fdb00d1b7c67c294/R/jointBetas_nofitting.R).

Validation covers both the Spire and FMN2 sample layouts with a small synthetic fixture, full-table calculation/export, agreement with a direct full-table FDR/group-estimate call, input paths with spaces, checkpoint reuse, rejection of incomplete or corrupted results, and submission of exactly one sbatch job. Local data and test artifacts are excluded from commits. Live CRG execution remains untested.

## Checks

```bash
scripts/run_R.sh tests/cluster_smoke.R
BETAS_DATASET=fmn2 scripts/run_R.sh tests/cluster_smoke.R
scripts/run_R.sh tests/analysis_checks.R
scripts/run_R.sh tests/analysis_ingest.R
BETAS_DATASET=fmn2 scripts/run_R.sh tests/analysis_ingest.R
scripts/run_R.sh tests/comparison_checks.R
```
