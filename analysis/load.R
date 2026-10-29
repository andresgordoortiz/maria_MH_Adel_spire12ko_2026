#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop('Usage: BETAS_DATASET=spire|spire_fastq|fmn2 scripts/run_R.sh analysis/run.R legacy|RUN/csv/DATASET_all_events.csv OUTPUT_DIR')
.libPaths(c(normalizePath('analysis/renv/library', mustWork = FALSE), .libPaths()))
source('analysis/common.R')
source('config/analysis.R')
stopifnot(as.character(getRversion()) == trimws(readLines('.R-version')))
dataset_title <- if (dataset == 'fmn2') 'FMN2 KO' else 'Spire1/2 DKO'
contrast_label <- paste(config$group_b, '−', config$group_a)
out <- args[2]
if (file.exists(out)) stop('Output already exists; choose a new directory: ', out)
input <- readRDS(config$filtered)
samples <- names(config$sample_groups)
group_a_columns <- which(config$sample_groups == config$group_a)
group_b_columns <- which(config$sample_groups == config$group_b)
sample_labels <- setNames(paste0(config$sample_groups, '_', ave(seq_along(samples), config$sample_groups, FUN = seq_along)), samples)
stopifnot(length(group_a_columns) >= 2, length(group_b_columns) >= 2,
  identical(names(input$PSI)[-(1:6)], samples), !anyDuplicated(input$PSI$EVENT),
  nrow(input$PSI) > 0, all(is.finite(as.matrix(input$PSI[, samples]))),
  all(as.matrix(input$PSI[, samples]) >= 0), all(as.matrix(input$PSI[, samples]) <= 100))
source_files <- config$filtered
raw_path <- if (dataset == 'spire_fastq') config$input else
  file.path('final_input_tables', if (dataset == 'spire') 'spiredko_INCLUSION_LEVELS_FULL-mm10.tab' else 'fmndko_INCLUSION_LEVELS_FULL-mm10.tab')
if (file.exists(raw_path)) {
  stopifnot(identical(input$provenance$input_md5, unname(tools::md5sum(raw_path))))
  source_files <- c(source_files, raw_path)
}
stats_cols <- c('FDR', 'invertedFDR', 'betasPsiA', 'betasPsiB', 'deltapsi')
if (args[1] == 'legacy') {
  stopifnot(dataset == 'spire')
  events <- input$PSI
  ref_stats <- NULL
  alignment <- list()
  for (kind in c('exons', 'introns', 'alt')) {
    path <- paste0('legacy/input_tables/spire_pdiff_', kind, '.csv')
    old <- as.data.frame(fread(path, check.names = FALSE))
    identified <- !is.na(old$EVENT) & nzchar(old$EVENT)
    stopifnot(nrow(old) == nrow(events), all(stats_cols %in% names(old)), any(identified))
    for (col in names(events)) stopifnot(isTRUE(all.equal(old[identified, col],
      events[identified, col], check.attributes = FALSE)))
    if (is.null(ref_stats)) ref_stats <- old[, stats_cols]
    else stopifnot(identical(ref_stats, old[, stats_cols]))
    alignment[[kind]] <- data.frame(family = kind, rows = nrow(old), identified = sum(identified), mismatches = 0)
    source_files <- c(source_files, path)
  }
  events <- cbind(events, ref_stats)
  provenance_label <- 'PROVISIONAL: legacy betAS; Spire1/2 DKO vs Cont'
} else {
  path <- args[1]
  run_dir <- dirname(dirname(path))
  if (!file.exists(file.path(run_dir, 'COMPLETE'))) stop('New simulation must have its COMPLETE marker')
  run <- readRDS(file.path(run_dir, 'run.rds'))
  stopifnot(identical(run$mode, 'full-table'), run$rows == nrow(input$PSI),
    identical(run$config$sample_groups, config$sample_groups), run$config$min_reads == 10,
    identical(run$config$group_a, config$group_a), identical(run$config$group_b, config$group_b))
  events <- as.data.frame(fread(path, check.names = FALSE))
  stopifnot(all(c(names(input$PSI), stats_cols) %in% names(events)),
    nrow(events) == nrow(input$PSI), identical(events$EVENT, input$PSI$EVENT))
  for (col in names(input$PSI)) stopifnot(isTRUE(all.equal(events[[col]], input$PSI[[col]], check.attributes = FALSE)))
  binary_match <- identical(run$filtered_md5, unname(tools::md5sum(config$filtered)))
  checkpoint <- file.path(run_dir, 'statistics.rds')
  if (!binary_match && !file.exists(checkpoint)) stop('Different filtered-object hash requires the completed statistics checkpoint')
  if (file.exists(checkpoint)) {
    saved <- readRDS(checkpoint)
    stopifnot(identical(saved$run, run), !anyDuplicated(saved$results$EVENT))
    for (col in names(events)) stopifnot(isTRUE(all.equal(events[[col]], saved$results[[col]],
      check.attributes = FALSE, tolerance = 1e-12)))
    stopifnot(identical(run$lock_md5, unname(tools::md5sum(file.path(run_dir, 'renv.lock')))))
    sheet <- read.csv(file.path(run_dir, 'samples.csv'))
    stopifnot(identical(sheet$sample, samples), identical(sheet$group, unname(config$sample_groups)))
    source_files <- c(source_files, checkpoint, file.path(run_dir, c('renv.lock', 'samples.csv')))
    rm(saved); invisible(gc())
  }
  source_files <- c(source_files, path, file.path(run_dir, c('COMPLETE', 'run.rds')))
  provenance_label <- paste(dataset_title, 'vs', config$group_a, '; final full-table betAS run')
}
events$direction <- factor(classify(events$FDR, events$deltapsi), levels = c('up', 'down', 'ndiff'))
events$event_type <- factor(ifelse(grepl('^MmuINT', events$EVENT), 'IR',
  ifelse(grepl('^MmuALTA', events$EVENT), 'Alt3', ifelse(grepl('^MmuALTD', events$EVENT), 'Alt5', 'EX'))),
  levels = c('EX', 'IR', 'Alt5', 'Alt3'))
stopifnot(!anyNA(events$event_type), all(abs(events$deltapsi - (events$betasPsiB - events$betasPsiA)) < 1e-8))
dir.create(file.path(out, 'figures'), recursive = TRUE)
dir.create(file.path(out, 'tables'))
dir.create(file.path(out, 'external_inputs'))
if (exists('binary_match')) write_table(data.frame(local_filtered_md5 = unname(tools::md5sum(config$filtered)),
  cluster_filtered_md5 = run$filtered_md5, binary_match = binary_match,
  validation = if (binary_match) 'identical filtered object + exact PSI/annotations' else
    'different serialized bytes; exact PSI/annotations + completed statistics/run/lock/sample checkpoints'), 'input_validation')
write_table(events, 'all_events')
write_table(events[events$direction != 'ndiff', ], 'differential_events')
if (exists('alignment')) write_table(bind_rows(alignment), 'legacy_alignment')
writeLines(dataset, file.path(out, 'DATASET'))
