message('Preparing Matt and rMAPS2 inputs')
parse_range <- function(x) {
  parts <- stringr::str_match(x, '^([^:]+):([0-9]+)-([0-9]+)$')
  data.frame(chr = parts[, 2], start = as.integer(parts[, 3]), end = as.integer(parts[, 4]))
}
external_counts <- list()
for (kind in levels(events$event_type)) {
  ix <- which(events$event_type == kind)
  ann <- annotation[ix, ]
  ann$DATASET <- events$direction[ix]
  ann$GENE <- events$GENE[ix]
  a <- parse_range(ann$CO_A); c1 <- parse_range(ann$CO_C1); c2 <- parse_range(ann$CO_C2)
  valid <- complete.cases(a, c1, c2) & ann$STRAND %in% c('+', '-') &
    a$start <= a$end & c1$start <= c1$end & c2$start <= c2$end &
    a$chr == c1$chr & a$chr == c2$chr
  valid[is.na(valid)] <- FALSE
  # Preserve the paper's custom rMAPS2 table format and VASTDB 1-based coordinates.
  rmaps <- data.frame(chr = a$chr, strand = ann$STRAND, exonStart = a$start, exonEnd = a$end,
    firstExonStart = pmin(c1$start, c2$start),
    firstExonEnd = ifelse(c1$start < c2$start, c1$end, c2$end),
    secondExonStart = pmax(c1$start, c2$start),
    secondExonEnd = ifelse(c1$start < c2$start, c2$end, c1$end))
  for (group in c('up', 'down', 'ndiff')) {
    rows <- which(valid & ann$DATASET == group)
    filename <- paste0(dataset, '_', kind, '_', group, '.txt')
    fwrite(rmaps[rows, ], file.path(out, 'external_inputs', filename), sep = '\t')
    write_table(ann[rows, c('EVENT', 'GENE', 'DATASET')], paste0('rmaps_', kind, '_', group, '_event_manifest'))
    external_counts[[length(external_counts) + 1L]] <- data.frame(tool = 'rMAPS2', type = kind, group = group,
      available = sum(ann$DATASET == group), exported = length(rows))
  }
  if (!kind %in% c('EX', 'IR')) next
  matt_ok <- valid & !is.na(ann$GENEID_ENSEMBL) & nzchar(ann$GENEID_ENSEMBL)
  rows <- which(matt_ok & ann$DATASET != 'ndiff')
  bg <- which(matt_ok & ann$DATASET == 'ndiff')
  set.seed(42)
  if (length(bg) > 10000) bg <- sample(bg, 10000)
  rows <- c(rows, bg)
  matt <- ann[rows, c('EVENT', 'START', 'END', 'SCAFFOLD', 'STRAND', 'UPSTRM_EX_BORDER',
    'DOSTRM_EX_BORDER', 'GENEID_ENSEMBL', 'DATASET')]
  stopifnot(!anyNA(matt[, c('EVENT', 'START', 'END', 'SCAFFOLD', 'STRAND', 'GENEID_ENSEMBL', 'DATASET')]),
    !anyDuplicated(matt$EVENT), all(matt$START <= matt$END))
  fwrite(matt, file.path(out, 'external_inputs', paste0(dataset, '_', kind, '_matt.tab')), sep = '\t', na = 'NA')
  for (group in c('up', 'down', 'ndiff')) external_counts[[length(external_counts) + 1L]] <-
    data.frame(tool = 'Matt', type = kind, group = group,
      available = sum(ann$DATASET == group), exported = sum(matt$DATASET == group))
}
write_table(bind_rows(external_counts), 'external_input_coverage')
input_names <- list.files(file.path(out, 'external_inputs'), pattern = '[.]txt$')
write.csv(data.frame(file = input_names, sha256 = vapply(input_names, function(n)
  digest::digest(file = file.path(out, 'external_inputs', n), algo = 'sha256'), character(1))),
  file.path(out, 'external_inputs', 'input_checksums.csv'), row.names = FALSE)
status <- rbind(status, data.frame(step = 'RBP_rMAPS2_inputs', status = 'ready; external server analysis needs user intervention'))

matt_dir <- Sys.getenv('MATT_OUTPUT_DIR', paste0('data/processed/', dataset, '_matt_final'))
if (!isTRUE(get0('prepare_only')) && file.exists(file.path(matt_dir, 'COMPLETE'))) {
  checksums <- read.csv(file.path(matt_dir, 'input_checksums.csv'))
  for (i in seq_len(nrow(checksums))) stopifnot(identical(checksums$sha256[i],
    digest::digest(file = file.path(out, 'external_inputs', checksums$file[i]), algo = 'sha256')))
  source('analysis/matt_figures.R')
  status <- rbind(status, data.frame(step = 'Matt_features_and_GC_gradients', status = 'complete; pinned paper container'))
} else status <- rbind(status, data.frame(step = 'Matt_features_and_GC_gradients',
  status = 'ready: bash analysis/run_matt.sh ANALYSIS_DIR NEW_MATT_DIR; then rerender with MATT_OUTPUT_DIR'))

rmaps_dir <- Sys.getenv('RMAPS_OUTPUT_DIR', paste0('data/processed/', dataset, '_rmaps'))
if (!isTRUE(get0('prepare_only')) && dir.exists(rmaps_dir)) {
  source('analysis/rbp_figures.R')
} else status <- rbind(status, data.frame(step = 'RBP_heatmaps_and_region_summaries',
  status = 'pending external rMAPS2 results; set RMAPS_OUTPUT_DIR and rerender'))
