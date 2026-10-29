source('scripts/common.R')
args <- commandArgs(TRUE)
if (length(args) != 1L) stop('Usage: 02_simulate.R RUN_DIRECTORY')
run_dir <- args[1]
run <- readRDS(file.path(run_dir, 'run.rds'))
stopifnot(identical(run$mode, 'full-table'), identical(run$config, config),
          identical(run$filtered_md5, unname(tools::md5sum(config$filtered))),
          identical(run$lock_md5, unname(tools::md5sum('renv.lock'))))
x <- readRDS(config$filtered)
stopifnot(identical(x$provenance$input_md5, unname(tools::md5sum(config$input))))
psi <- x$PSI; qual <- x$Qual
stopifnot(nrow(psi) == run$rows, identical(psi$EVENT, qual$EVENT),
          !anyNA(psi$EVENT), !anyDuplicated(psi$EVENT))
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
# Each upstream call sees the entire filtered table. Checkpoints are per
# statistic, not per subset of rows, and must belong to this exact run.
stage <- function(name, calculate) {
  path <- file.path(run_dir, paste0(name, '.rds'))
  if (file.exists(path)) {
    message('Loading completed ', name, ' stage')
    saved <- readRDS(path)
    stopifnot(identical(saved$run, run), identical(saved$results$EVENT, psi$EVENT))
    return(saved$results)
  }
  message(format(Sys.time(), tz = 'UTC', usetz = TRUE), ': starting ', name,
          ' on all ', nrow(psi), ' events')
  result <- calculate()
  stopifnot(identical(result$EVENT, psi$EVENT))
  atomic_rds(list(results = result, run = run,
    session = capture.output(sessionInfo())), path)
  message('Completed and saved ', name)
  invisible(gc())
  result
}
fdr <- stage('FDR', function() do.call(betAS::prepareTableVolcanoFDR,
                                     c(common, list(nsim = config$nsim))))
pdiff <- stage('Pdiff', function() do.call(betAS::prepareTableVolcano, common))
stopifnot(isTRUE(all.equal(fdr$deltapsi, pdiff$deltapsi, tolerance = 1e-12)))
fdr$Pdiff <- pdiff$Pdiff[match(fdr$EVENT, pdiff$EVENT)]
rm(pdiff); invisible(gc())
fstat <- stage('Fstat', function() do.call(betAS::prepareTableVolcanoFstat, common))
stopifnot(isTRUE(all.equal(fdr$deltapsi, fstat$deltapsi, tolerance = 1e-12)))
fdr$Fstat <- fstat$Fstat[match(fdr$EVENT, fstat$EVENT)]
rm(fstat); invisible(gc())
stopifnot(all(is.finite(fdr$deltapsi)), all(is.finite(fdr$FDR)),
          all(fdr$FDR >= 0 & fdr$FDR <= 1), all(is.finite(fdr$Pdiff)),
          all(fdr$Pdiff >= 0 & fdr$Pdiff <= 1), !anyNA(fdr$Fstat))
# Infinite Fstat is possible if simulated within-group variation is zero.
atomic_rds(list(results = fdr, run = run,
               session = capture.output(sessionInfo())),
           file.path(run_dir, 'statistics.rds'))
cat('Finished full-table simulation:', nrow(fdr), 'events\n')
