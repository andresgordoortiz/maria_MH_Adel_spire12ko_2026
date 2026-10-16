# Spire1/2 DKO: CRG betAS workflow

This repository contains the environment, configuration and execution scripts only. Input data, legacy exports, audit results, generated sample sheets, filtered tables, logs and simulation results are not committed.

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

Setup then installs/restores R and packages, validates the six expected sample/quality column pairs, generates `metadata/samples.csv` from the group mapping in `config/analysis.R`, and filters the existing inclusion table. It creates all required output directories. The allocated node needs internet access for installation and shared access to the project and home directory. The input is read without modification.

The filter is `betAS::filterEvents(..., N=10)`: no missing PSI and corrected inclusion+exclusion coverage >=10 in **each of the six samples**. An independent coverage check must agree exactly. ANN annotation-module exons are retained and exported separately. The filtered FULL table, paired PSI/Qual RDS and chunk manifest are generated in `data/processed/`; QC reports are generated in `metadata/`.

The configured contrast is Cont_a/b/c (A) versus DKO_a/b/c (B). Positive delta means higher inclusion/site use/retention in DKO. The labels do not establish independent animals, pooling, pairing or batches.

If using an existing CRG R 4.5.2 instead of rig, activate that installation and set `RSCRIPT` to its executable. Run `"$RSCRIPT" scripts/setup.R` followed by `scripts/run_R.sh scripts/01_filter.R`. For a non-portable R build, `RENV_CONFIG_REPOS_OVERRIDE=https://cloud.r-project.org` selects source packages during restore; system build dependencies may be required.

## Submit simulations

After preparation succeeds, from the login/submission node:

```bash
scripts/submit_crg.sh results/spire_N10_nsim1000_run01
```

This invokes **sbatch** for the simulation array and a dependent collection job. It prints job IDs and a monitoring command. Once submitted, jobs continue after disconnection. Existing run directories are refused; use a new name for each analysis.

Defaults in `config/analysis.R`: 500 beta points, 1,000 null simulations, equal sample weighting, fixed betAS internal seed (`seed=TRUE`, internally 21122023), and 500 events per chunk. The pinned Fstat wrapper hardcodes 500 points, so npoints is validated as 500. All coverage-passing events are simulated. At 289,583 passing events, this produces 580 chunks, with at most 20 concurrent tasks per array.

Simulation jobs request one CPU, 8 GB and `vlong`/48 hours; collection requests one CPU, 16 GB and `short`/one hour. These QoS names follow the reference CRG scripts but have not been checked against live CRG policy. Adjust the SBATCH directives if your account differs. No account or partition is assumed. `MAX_CONCURRENT` and `ARRAY_WINDOW` configure submission; keep ARRAY_WINDOW within the site's MaxArraySize.

Jobs never install packages. Each writes an atomic RDS chunk. Collection requires every expected chunk and verifies event IDs/order, configuration and input-object/lockfile hashes before exporting. Run configuration, sample mapping, lockfile and session information accompany the outputs. A `COMPLETE` marker is written only after successful collection. Do not alter configuration or processed inputs during an active run.

For a failed chunk, rerun `scripts/run_R.sh scripts/02_simulate.R CHUNK_ID results/RUN` inside an allocation, using the 1-based global chunk ID, then rerun collection with `scripts/run_R.sh scripts/03_collect.R results/RUN`.

## Results and statistical interpretation

Under `results/RUN/csv/`:

- `spire_all_events.csv`: all passing events, individual PSI and new statistics.
- `spire_statistics_{exons,introns,alt,ann}.csv`: complete family-specific results.
- `differential_spire_*.csv`: upstream FDR <=0.05 and |delta PSI| >=0.10.
- `spire_*_combination.csv`: the same limits plus Pdiff >=0.95.

These are explicit new selection defaults, not a reconstruction of the old unknown thresholds. Full results retain all unselected events. No heatmap clusters or protein-impact predictions are generated.

Pdiff, Fstat and upstream FDR are joined by EVENT after checking that their delta PSI agrees. Delta PSI is on a 0–1 scale; 0.10 means 10 percentage points. The betAS variable named `FDR` is a per-event simulated null tail fraction, **not a BH-adjusted value**; its empirical resolution with 1,000 simulations is 0.001. See the [pinned source](https://github.com/DiseaseTranscriptomicsLab/betAS/blob/cd37737c75ab17a64652baf2fdb00d1b7c67c294/R/jointBetas_nofitting.R).

The previous local end-to-end check passed on 19 real events spanning all nine passing types, including deterministic reruns and rejection of missing chunks/corrupted event IDs. Those local data and test artifacts are not part of this cluster-sync commit. Live CRG execution remains untested.
