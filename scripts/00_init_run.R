source('scripts/common.R')
a <- commandArgs(TRUE)
if (length(a) != 1L) stop('Usage: 00_init_run.R NEW_RUN_DIRECTORY')
if (dir.exists(a[1])) stop('Run directory already exists; use a new name')
x <- readRDS(config$filtered)
stopifnot(identical(x$provenance$config, config),
  identical(x$provenance$input_md5, unname(tools::md5sum(config$input))),
  nrow(x$PSI) > 0L, identical(x$PSI$EVENT, x$Qual$EVENT),
  !anyNA(x$PSI$EVENT), !anyDuplicated(x$PSI$EVENT))
run <- list(config = config, mode = 'full-table', rows = nrow(x$PSI),
            filtered_md5 = unname(tools::md5sum(config$filtered)),
            lock_md5 = unname(tools::md5sum('renv.lock')),
            created_utc = format(Sys.time(), tz = 'UTC', usetz = TRUE))
atomic_rds(run, file.path(a[1], 'run.rds'))
file.copy('renv.lock', file.path(a[1], 'renv.lock'))
file.copy(config$samples, file.path(a[1], 'samples.csv'))
writeLines(capture.output(dput(config)), file.path(a[1], 'config.R'))
cat('Initialized full-table run with', nrow(x$PSI), 'events in', a[1], '\n')
