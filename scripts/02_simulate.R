source('scripts/common.R')
args <- commandArgs(TRUE)
if (length(args) != 2L) stop('Usage: 02_simulate.R CHUNK_ID RUN_DIRECTORY')
id <- as.integer(args[1]); run_dir <- args[2]
run <- readRDS(file.path(run_dir, 'run.rds'))
stopifnot(identical(run$config, config), identical(run$filtered_md5, unname(tools::md5sum(config$filtered))),
          identical(run$lock_md5, unname(tools::md5sum('renv.lock'))))
x <- readRDS(config$filtered)
manifest <- run$manifest
stopifnot(length(id) == 1L, !is.na(id), id %in% manifest$chunk)
row <- manifest[manifest$chunk == id, ]
ix <- seq.int(row$start, row$end)
psi <- x$PSI[ix, , drop = FALSE]; qual <- x$Qual[ix, , drop = FALSE]
stopifnot(identical(psi$EVENT, qual$EVENT), !anyNA(psi$EVENT))
A <- x$metadata$sample[x$metadata$group == config$group_a]
B <- x$metadata$sample[x$metadata$group == config$group_b]
reference <- new.env()
data('maxDevSimulationN100', package = 'betAS', envir = reference)
common <- list(psitable = psi, qualtable = qual, npoints = config$npoints,
  colsA = match(A, names(psi)), colsB = match(B, names(psi)),
  labA = config$group_a, labB = config$group_b,
  basalColor = '#777777', interestColor = '#CC3300',
  maxDevTable = reference$maxDevSimulationN100,
  seed = config$seed, CoverageWeight = config$coverage_weight)
fdr <- do.call(betAS::prepareTableVolcanoFDR, c(common, list(nsim = config$nsim)))
pdiff <- do.call(betAS::prepareTableVolcano, common)
fstat <- do.call(betAS::prepareTableVolcanoFstat, common)
stopifnot(identical(fdr$EVENT, psi$EVENT), identical(pdiff$EVENT, psi$EVENT),
          identical(fstat$EVENT, psi$EVENT),
          isTRUE(all.equal(fdr$deltapsi, pdiff$deltapsi, tolerance = 1e-12)),
          isTRUE(all.equal(fdr$deltapsi, fstat$deltapsi, tolerance = 1e-12)))
# Join explicitly by EVENT, never by original row names or global vector position.
fdr$Pdiff <- pdiff$Pdiff[match(fdr$EVENT, pdiff$EVENT)]
fdr$Fstat <- fstat$Fstat[match(fdr$EVENT, fstat$EVENT)]
stopifnot(all(is.finite(fdr$deltapsi)), all(is.finite(fdr$FDR)),
          all(fdr$FDR >= 0 & fdr$FDR <= 1), all(is.finite(fdr$Pdiff)),
          all(fdr$Pdiff >= 0 & fdr$Pdiff <= 1), !anyNA(fdr$Fstat))
# Infinite Fstat is possible if simulated within-group variation is zero.
atomic_rds(list(results = fdr, chunk = id, run = run,
               session = capture.output(sessionInfo())),
           file.path(run_dir, 'chunks', sprintf('%05d.rds', id)))
cat('Finished chunk', id, ':', nrow(fdr), 'events\n')
