source('scripts/common.R')
a <- commandArgs(TRUE)
if (length(a) != 1L) stop('Usage: 03_collect.R RUN_DIRECTORY')
run_dir <- a[1]
if (file.exists(file.path(run_dir, 'COMPLETE'))) unlink(file.path(run_dir, 'COMPLETE'))
run <- readRDS(file.path(run_dir, 'run.rds'))
stopifnot(identical(run$config, config), identical(run$filtered_md5, unname(tools::md5sum(config$filtered))),
          identical(run$lock_md5, unname(tools::md5sum('renv.lock'))))
x <- readRDS(config$filtered)
stopifnot(identical(run$mode, 'full-table'))
p <- file.path(run_dir, 'statistics.rds')
if (!file.exists(p)) stop('Missing full-table statistics.rds; simulation must finish first')
z <- readRDS(p)
stopifnot(identical(z$run, run))
all <- z$results; rownames(all) <- NULL
stopifnot(nrow(all) == run$rows, !anyDuplicated(all$EVENT),
          identical(all$EVENT, x$PSI$EVENT))
out <- file.path(run_dir, 'csv'); dir.create(out, showWarnings = FALSE)
write_csv(all, file.path(out, paste0(dataset, '_all_events.csv')))
for (kind in c('exons','introns','alt','ann')) {
  z <- all[family(all$EVENT, all$COMPLEX) == kind, , drop = FALSE]
  write_csv(z, file.path(out, paste0(dataset, '_statistics_', kind, '.csv')))
  differential <- z[z$FDR <= config$fdr_max & abs(z$deltapsi) >= config$abs_delta_min, , drop = FALSE]
  write_csv(differential, file.path(out, paste0('differential_', dataset, '_', kind, '.csv')))
  selected <- z[z$Pdiff >= config$pdiff_min & z$FDR <= config$fdr_max &
                abs(z$deltapsi) >= config$abs_delta_min,
              c('GENE','EVENT','COORD','deltapsi','Pdiff','Fstat','FDR'), drop = FALSE]
  write_csv(selected, file.path(out, paste0(dataset, '_', kind, '_combination.csv')))
}
writeLines(capture.output(sessionInfo()), file.path(run_dir, 'collection_sessionInfo.txt'))
writeLines('Full-table results verified; CSV export complete.', file.path(run_dir, 'COMPLETE'))
cat('Exported', nrow(all), 'events to', out, '\n')
