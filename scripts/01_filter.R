source('scripts/common.R')
if (!file.exists(config$input)) stop('Inclusion table not found: ', config$input)
metadata_dir <- dirname(config$samples)
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
metadata <- data.frame(sample = names(config$sample_groups),
  group = unname(config$sample_groups),
  replicate_label = if (dataset == 'fmn2') c('1','2','1','2') else
    sub('.*_', '', names(config$sample_groups)),
  metadata_source = if (dataset == 'fmn2') 'GSE103718_and_reference_run_manifest' else
    'configured_from_inclusion_header', stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(metadata$sample),
          identical(metadata$sample, names(config$sample_groups)),
          all(metadata$group %in% c(config$group_a, config$group_b)),
          all(table(factor(metadata$group, levels=c(config$group_a,config$group_b))) >= 2L))
raw <- betAS::getDataset(config$input, tool = 'vast-tools')
expected <- c('GENE','EVENT','COORD','LENGTH','FullCO','COMPLEX',
              as.vector(rbind(metadata$sample, paste0(metadata$sample, '.Q'))))
stopifnot(identical(names(raw), expected), !anyNA(raw$EVENT), !anyDuplicated(raw$EVENT))
write_csv(metadata, config$samples)
events <- betAS::getEvents(raw, tool = 'vast-tools')
stopifnot(identical(events$Samples, metadata$sample), identical(events$PSI$EVENT, events$Qual$EVENT))
filtered <- betAS::filterEvents(events, N = config$min_reads)
stopifnot(!anyNA(filtered$PSI), identical(filtered$PSI$EVENT, filtered$Qual$EVENT),
          !anyDuplicated(filtered$PSI$EVENT), nrow(filtered$PSI) > 0)
# Independently check the coverage filter, including sentinel/negative counts.
counts <- lapply(events$Qual[-(1:6)], function(x) {
  pairs <- strsplit(sub('.*@', '', x), ',', fixed = TRUE)
  matrix(as.numeric(unlist(pairs)), ncol = 2, byrow = TRUE)
})
coverage <- do.call(cbind, lapply(counts, rowSums))
valid_counts <- Reduce(`&`, lapply(counts, function(x) rowSums(!is.finite(x) | x < 0) == 0))
complete <- complete.cases(events$PSI[-(1:6)])
keep <- complete & valid_counts & rowSums(coverage >= config$min_reads, na.rm = TRUE) == nrow(metadata)
stopifnot(identical(filtered$PSI$EVENT, events$PSI$EVENT[keep]))
n10_counts <- table(filtered$PSI$COMPLEX)
all_n10_ids <- filtered$PSI$EVENT
# Keep ANN annotation-module exons too; export separately during collection.
selected <- filtered$PSI$COMPLEX %in% config$analysis_types
filtered$PSI <- filtered$PSI[selected, , drop = FALSE]
filtered$Qual <- filtered$Qual[selected, , drop = FALSE]
filtered$EventsPerType <- table(filtered$PSI$COMPLEX)
filtered$metadata <- metadata
filtered$provenance <- list(input_md5 = unname(tools::md5sum(config$input)),
  config = config, R = R.version.string, betAS = as.character(packageVersion('betAS')),
  created_utc = format(Sys.time(), tz = 'UTC', usetz = TRUE))
atomic_rds(filtered, config$filtered)
raw_filtered <- raw[match(all_n10_ids, raw$EVENT), ]
# Restore original quality header spelling for a genuine FULL table.
names(raw_filtered) <- sub('[.]Q$', '-Q', names(raw_filtered))
write.table(raw_filtered, file.path(dirname(config$filtered), paste0(dataset, '_N10_INCLUSION_LEVELS_FULL-mm10.tab')),
            sep = '\t', quote = FALSE, row.names = FALSE, na = 'NA')
qc <- data.frame(sample = metadata$sample, group = metadata$group,
  total_events = nrow(raw), missing_PSI = colSums(is.na(events$PSI[-(1:6)])),
  coverage_ge10 = colSums(coverage >= config$min_reads, na.rm = TRUE))
write_csv(qc, file.path(metadata_dir, 'sample_qc.csv'))
write_csv(data.frame(COMPLEX = names(n10_counts), retained_N10 = as.integer(n10_counts),
                    included_in_simulation = names(n10_counts) %in% config$analysis_types), file.path(metadata_dir, 'filter_counts.csv'))
writeLines(capture.output(sessionInfo()), file.path(metadata_dir, 'filter_sessionInfo.txt'))
cat(nrow(raw), 'input events;', length(all_n10_ids), 'pass N10;', nrow(filtered$PSI), 'AS events\n')
