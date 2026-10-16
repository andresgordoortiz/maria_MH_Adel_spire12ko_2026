.libPaths(c(normalizePath('analysis/renv/library', mustWork = FALSE), .libPaths()))
source('analysis/common.R')
stopifnot(identical(classify(c(.05, .05, .051, 0, 0), c(.1, -.1, .9, 0, .099)),
  c('up', 'down', 'ndiff', 'ndiff', 'ndiff')))
for (bad in list(c(NA, .2), c(-.1, .2), c(1.1, .2)))
  stopifnot(inherits(try(classify(bad, c(0, 0)), silent = TRUE), 'try-error'))
mat <- rbind(c(1, 1, 1, 5, 6, 7), c(3, 3, 3, 3, 3, 3), c(8, 7, 6, 1, 1, 1))
scaled <- robust_scale(mat, 1)
reference <- t(apply(mat, 1, function(x) {
  spread <- mad(x); if (spread == 0) spread <- 1
  pmax(pmin((x - median(x)) / spread, 1), -1)
}))
stopifnot(identical(dim(scaled), dim(mat)), isTRUE(all.equal(scaled, reference)),
  all(scaled[2, ] == 0), length(cluster_rows(scaled[-2, ])$order) == 2)
stopifnot(inherits(try(cluster_rows(scaled), silent = TRUE), 'try-error'))

out <- tempfile('spire-analysis-check-')
dir.create(file.path(out, 'external_inputs'), recursive = TRUE)
dir.create(file.path(out, 'tables'))
dir.create(file.path(out, 'figures'))
main_plots <- list(); figure_names <- character(); source_files <- character()
provenance_label <- 'Synthetic check; not biological results'
events <- data.frame(EVENT = c('MmuEXtest1', 'MmuEXtest2', 'MmuEXmissing', 'MmuINTtest1'),
  GENE = c('A', 'B', 'C', 'D'), event_type = factor(c('EX', 'EX', 'EX', 'IR'), levels = c('EX', 'IR')),
  direction = factor(c('up', 'down', 'ndiff', 'ndiff'), levels = c('up', 'down', 'ndiff')))
annotation <- data.frame(EVENT = events$EVENT, GENE = events$GENE,
  CO_A = c('chr1:100-199', 'chr2:100-199', NA, 'chr3:20-29'),
  CO_C1 = c('chr1:10-19', 'chr2:300-399', NA, 'chr3:10-19'),
  CO_C2 = c('chr1:300-399', 'chr2:10-19', NA, 'chr3:30-39'),
  STRAND = c('+', '-', NA, '+'), SCAFFOLD = c('chr1', 'chr2', NA, 'chr3'),
  START = c(100, 100, NA, 20), END = c(199, 199, NA, 29),
  GENEID_ENSEMBL = c('ENSMUSG1', 'ENSMUSG2', NA, 'ENSMUSG3'),
  UPSTRM_EX_BORDER = c(19, 19, NA, 19), DOSTRM_EX_BORDER = c(300, 300, NA, 30))
old_env <- Sys.getenv(c('MATT_OUTPUT_DIR', 'RMAPS_OUTPUT_DIR'), unset = NA)
Sys.setenv(MATT_OUTPUT_DIR = file.path(out, 'no-matt'), RMAPS_OUTPUT_DIR = file.path(out, 'no-rmaps'))
status <- data.frame(step = character(), status = character())
source('analysis/external.R')
down <- fread(file.path(out, 'external_inputs', 'spire_EX_down.txt'))
stopifnot(nrow(down) == 1, down$strand == '-', down$firstExonStart == 10,
  down$firstExonEnd == 19, down$secondExonStart == 300, down$secondExonEnd == 399,
  down$exonStart == 100, down$exonEnd == 199)
coverage <- read.csv(file.path(out, 'tables', 'external_input_coverage.csv'))
stopifnot(coverage$exported[coverage$tool == 'rMAPS2' & coverage$type == 'EX' & coverage$group == 'ndiff'] == 0)
matt <- fread(file.path(out, 'external_inputs', 'spire_EX_matt.tab'))
stopifnot(identical(matt$EVENT, c('MmuEXtest1', 'MmuEXtest2')), all(matt$UPSTRM_EX_BORDER < matt$START),
  all(matt$DOSTRM_EX_BORDER > matt$END))
rmaps_dir <- file.path(out, 'rmaps')
dir.create(file.path(rmaps_dir, 'EX'), recursive = TRUE)
file.copy(file.path(out, 'external_inputs', 'input_checksums.csv'), rmaps_dir)
fixture <- data.frame(RBP = c('RBP1.motif1', 'RBP1.motif2', 'RBP2.motif1'))
for (i in 1:8) fixture[[paste0('R', i)]] <- c(.01, .02, .5)
fixture$R1[1] <- 0
for (file in c('pVal.up.vs.bg.RNAmap.txt', 'pVal.dn.vs.bg.RNAmap.txt'))
  fwrite(fixture, file.path(rmaps_dir, 'EX', file), sep = '\t')
source('analysis/rbp_figures.R')
rbp <- read.csv(file.path(out, 'tables', 'RBP_minimum_motif_pvalues.csv'))
stopifnot(nrow(rbp) == 32, rbp$pvalue[rbp$RBP == 'RBP1' & rbp$region == 'R1' & rbp$direction == 'up'] == 0,
  file.exists(file.path(out, 'figures', 'Fig2H_S6_RBP_EX.pdf')))
sequence <- paste0(strrep('A', 80), strrep('GC', 80), strrep('G', 80))
ex <- data.frame(EVENT = paste0('EX', 1:3), DATASET = c('up', 'down', 'ndiff'),
  SEQ_UPEXON = strrep('ATGC', 100), SEQ_UPINTRON = sequence, SEQ_EXON = strrep('ATGC', 30),
  SEQ_DOINTRON = sequence, SEQ_DOEXON = strrep('ATGC', 100))
intr <- data.frame(EVENT = paste0('IR', 1:3), DATASET = c('up', 'down', 'ndiff'),
  SEQ_LONGESTUPEXON = strrep('ATGC', 100), SEQ_INTRON = sequence, SEQ_LONGESTDOEXON = strrep('ATGC', 100))
matt_tables <- list(EX = ex, IR = intr)
source('analysis/gc_gradient.R')
gradient <- build_event_intron_gradient(strrep('A', 400), sequence, strrep('A', 400), make_intron_layout('fine'))
stopifnot(all(gradient[17:24] == 0), all(gradient[65:72] == 100))
for (name in names(old_env)) if (is.na(old_env[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, setNames(list(old_env[[name]]), name))
unlink(out, recursive = TRUE)
cat('Analysis thresholds, scaling, clustering, strand/export, RBP parser and GC-end checks passed.\n')
