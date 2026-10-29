dataset <- Sys.getenv('BETAS_DATASET', unset = 'spire')
if (!dataset %in% c('spire', 'spire_fastq', 'fmn2'))
  stop('BETAS_DATASET must be spire, spire_fastq or fmn2')
config <- list(
  input = Sys.getenv('SPIRE_INCLUSION_TABLE', unset =
    '/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/notebooks/inclusion_tables/spiredko_INCLUSION_LEVELS_FULL-mm10.tab'),
  samples = 'metadata/samples.csv',
  # Expected columns and group assignment; biological replicate identity is unverified.
  sample_groups = setNames(rep(c('Cont', 'DKO'), each = 3L),
    paste0('Oocytes_FG_Spire12_', rep(c('Cont', 'DKO'), each = 3L), '_', rep(letters[1:3], 2))),
  filtered = 'data/processed/spire_N10.rds',
  min_reads = 10L,
  analysis_types = c('S','C1','C2','C3','MIC','Alt3','Alt5','IR','ANN'),
  group_a = 'Cont', group_b = 'DKO',
  npoints = 500L, nsim = 1000L,
  # betAS seed=TRUE fixes its internal seed to 21122023 (not set.seed(42)).
  seed = TRUE, coverage_weight = FALSE,
  chunk_size = 500L, # Legacy filter-manifest setting; unused by single-job submission.
  pdiff_min = 0.95, fdr_max = 0.05, abs_delta_min = 0.10,
  betas_sha = 'cd37737c75ab17a64652baf2fdb00d1b7c67c294'
)

if (dataset == 'spire_fastq') {
  samples <- read.delim('config/spire_fastq.tsv', check.names = FALSE)
  config$input <- Sys.getenv('SPIRE_FASTQ_INCLUSION_TABLE', unset =
    'results/spire_fastq/vast/inclusion_tables/spire_fastq_INCLUSION_LEVELS_FULL-mm10.tab')
  config$samples <- 'metadata/spire_fastq/samples.csv'
  config$sample_groups <- setNames(samples$condition, samples$sample)
  config$filtered <- 'data/processed/spire_fastq/spire_fastq_N10.rds'
}

if (dataset == 'fmn2') {
  config$input <- Sys.getenv('FMN2_INCLUSION_TABLE', unset =
    '/users/aaljord/aaljord/Al Jord Lab - CRG/Lab Members/Andres Gordo Ortiz/Data/Oocyte/FMN2KO and SPIRE/FMN2KO SPIRE sent to MHV 14 02 2025/oocyte_fmn2_spire_report_andresgordo/raw_data/fmndko_INCLUSION_LEVELS_FULL-mm10.tab')
  config$samples <- 'metadata/fmn2/samples.csv'
  # GSE103718: Fmn7/Fmn8 are +/- pools; Fmn9/Fmn10 are -/- pools.
  config$sample_groups <- setNames(c('Het','Het','KO','KO'), c(
    'SRR6026682_SRR6026683_SRR6026684_merged',
    'SRR6026685_SRR6026686_SRR6026687_merged',
    'SRR6026688_SRR6026689_SRR6026690_merged',
    'SRR6026691_SRR6026692_SRR6026693_merged'))
  config$filtered <- 'data/processed/fmn2/fmn2_N10.rds'
  config$group_a <- 'Het'
  config$group_b <- 'KO'
}
