source('config/analysis.R')
stopifnot(config$npoints == 500L, config$nsim >= 1L, config$nsim <= 10000L,
          config$nsim == as.integer(config$nsim), config$chunk_size >= 1L,
          config$chunk_size == as.integer(config$chunk_size), config$min_reads >= 0,
          config$fdr_max >= 0, config$fdr_max <= 1,
          config$pdiff_min >= 0, config$pdiff_min <= 1,
          config$abs_delta_min >= 0, config$abs_delta_min <= 1)
stopifnot(as.character(getRversion()) == trimws(readLines('.R-version')))
stopifnot(requireNamespace('betAS', quietly = TRUE))
desc <- packageDescription('betAS')
stopifnot(identical(desc$RemoteSha, config$betas_sha))
atomic_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, '.', Sys.getpid(), '.tmp')
  saveRDS(x, tmp)
  if (!file.rename(tmp, path)) stop('Cannot publish ', path)
}
write_csv <- function(x, path) write.csv(x, path, row.names = FALSE, na = 'NA')
family <- function(events, complex) {
  ans <- ifelse(grepl('^MmuEX', events), 'exons',
    ifelse(grepl('^MmuINT', events), 'introns',
      ifelse(grepl('^MmuALT[AD]', events), 'alt', NA_character_)))
  ans[complex == 'ANN'] <- 'ann'
  if (anyNA(ans)) stop('Unexpected event family')
  ans
}
