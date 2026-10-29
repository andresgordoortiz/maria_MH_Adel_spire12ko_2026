# Exercise the new-result loader without biological data or a simulation.
.libPaths(c(normalizePath('analysis/renv/library'), .libPaths()))
source('config/analysis.R')
root <- normalizePath('.')
work <- tempfile('spire-ingest-'); dir.create(work)
dir.create(file.path(work, 'analysis')); dir.create(file.path(work, 'config'))
files <- list.files('analysis', full.names = TRUE)
files <- files[!file.info(files)$isdir]
invisible(file.copy(files, file.path(work, 'analysis')))
invisible(file.copy(c('.R-version', 'renv.lock'), work))
writeLines(c(paste0('dataset <- "', dataset, '"'), 'config <-', capture.output(dput(config))), file.path(work, 'config/analysis.R'))
psi <- data.frame(GENE = paste0('Gene', 1:4),
  EVENT = c('MmuEX0000001', 'MmuINT0000001', 'MmuALTD0000001-1/2', 'MmuALTA0000001-1/2'),
  COORD = 'chr1:100-199', LENGTH = 100, FullCO = 'chr1:1,100-199,300', COMPLEX = c('S', 'IR', 'Alt5', 'Alt3'))
values <- rbind(c(20, 21, 19, 40, 41, 39), c(60, 59, 61, 40, 39, 41),
  c(40, 39, 41, 50, 49, 51), c(50, 51, 49, 50, 51, 49))
columns <- if (dataset == 'fmn2') c(1, 2, 4, 5) else 1:6
for (i in seq_along(config$sample_groups)) psi[[names(config$sample_groups)[i]]] <- values[, columns[i]]
dir.create(file.path(work, dirname(config$filtered)), recursive = TRUE)
saveRDS(list(PSI = psi), file.path(work, config$filtered))
dir.create(file.path(work, 'results/new/csv'), recursive = TRUE)
run <- list(mode = 'full-table', rows = 4L, config = config,
  filtered_md5 = unname(tools::md5sum(file.path(work, config$filtered))))
saveRDS(run, file.path(work, 'results/new/run.rds'))
results <- psi
results$FDR <- c(.01, 0, .05, .4); results$invertedFDR <- 1 - results$FDR
results$betasPsiA <- c(.2, .6, .4, .5); results$betasPsiB <- c(.4, .4, .5, .5)
results$deltapsi <- c(.2, -.2, .1, 0)
csv <- file.path(work, 'results/new/csv/spire_all_events.csv')
write.csv(results, csv, row.names = FALSE)
setwd(work)
execute <- function(output, expected) {
  code <- system2(file.path(R.home('bin'), 'Rscript'), c('--vanilla', 'analysis/run.R',
    'results/new/csv/spire_all_events.csv', output),
    env = paste0('R_LIBS=', shQuote(paste(.libPaths(), collapse = .Platform$path.sep))),
    stdout = file.path(work, paste0(basename(output), '.log')),
    stderr = file.path(work, paste0(basename(output), '.log')))
  stopifnot(code == expected)
}
execute('results/no-marker', 1L)
writeLines('Synthetic complete fixture', 'results/new/COMPLETE')
execute('results/valid', 0L)
stopifnot(file.exists('results/valid/COMPLETE'), file.exists('results/valid/index.html'),
  nrow(read.csv('results/valid/tables/differential_events.csv')) == 3)
run$filtered_md5 <- strrep('0', 32)
run$lock_md5 <- unname(tools::md5sum('renv.lock'))
saveRDS(run, 'results/new/run.rds')
saveRDS(list(run = run, results = results), 'results/new/statistics.rds')
invisible(file.copy('renv.lock', 'results/new/renv.lock'))
write.csv(data.frame(sample = names(config$sample_groups), group = unname(config$sample_groups)), 'results/new/samples.csv', row.names = FALSE)
execute('results/metadata-serialization-differs', 0L)
corrupt <- results; corrupt$FDR[1] <- .9
saveRDS(list(run = run, results = corrupt), 'results/new/statistics.rds')
execute('results/corrupt-checkpoint', 1L)
saveRDS(list(run = run, results = results), 'results/new/statistics.rds')
write.csv(results[4:1, ], csv, row.names = FALSE)
execute('results/wrong-order', 1L)
results[[names(config$sample_groups)[1]]][1] <- 22
write.csv(results, csv, row.names = FALSE)
execute('results/wrong-PSI', 1L)
setwd(root)
unlink(work, recursive = TRUE)
cat(dataset, ': final-result ingest, sample layout, metadata-only hash differences and corruption rejection passed.\n')
