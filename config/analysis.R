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
