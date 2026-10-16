source('scripts/common.R')
a <- commandArgs(TRUE)
if (length(a) != 1L) stop('Usage: 00_init_run.R NEW_RUN_DIRECTORY')
if (dir.exists(a[1])) stop('Run directory already exists; use a new name')
x <- readRDS(config$filtered)
stopifnot(identical(x$provenance$config, config),
  identical(x$provenance$input_md5, unname(tools::md5sum(config$input))))
manifest <- read.csv('data/processed/chunks.csv')
stopifnot(identical(manifest$chunk, seq_len(nrow(manifest))), manifest$start[1] == 1L,
          tail(manifest$end, 1) == nrow(x$PSI), all(manifest$start[-1] == head(manifest$end, -1) + 1L))
run <- list(config = config, manifest = manifest,
            filtered_md5 = unname(tools::md5sum(config$filtered)),
            lock_md5 = unname(tools::md5sum('renv.lock')),
            created_utc = format(Sys.time(), tz = 'UTC', usetz = TRUE))
atomic_rds(run, file.path(a[1], 'run.rds'))
file.copy('renv.lock', file.path(a[1], 'renv.lock'))
file.copy(config$samples, file.path(a[1], 'samples.csv'))
writeLines(capture.output(dput(config)), file.path(a[1], 'config.R'))
cat('Initialized', nrow(manifest), 'chunks in', a[1], '\n')
