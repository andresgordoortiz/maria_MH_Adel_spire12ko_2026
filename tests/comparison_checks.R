.libPaths(c(normalizePath('analysis/renv/library'), .libPaths()))
source('analysis/common.R')
source('analysis/comparison_helpers.R')
universe <- as.character(1:100)
a <- as.character(1:5); b <- as.character(3:7)
x <- overlap_test(a, b, universe)
stopifnot(x$shared == 3, x$expected == .25,
  isTRUE(all.equal(x$p.value, phyper(2, 5, 95, 5, lower.tail = FALSE))))
stopifnot(inherits(try(overlap_test(c('outside', a), b, universe), silent = TRUE), 'try-error'))
fixture <- data.frame(GENE = c('A', 'B', 'C', 'D'), EVENT = c('EX1', 'INT1', 'ALT1-1/2', 'ALT1-2/2'),
  COORD = 'chr1:1-10', COMPLEX = c('S', 'IR', 'Alt3', 'Alt3'), event_type = c('EX', 'IR', 'Alt3', 'Alt3'),
  FDR = c(.01, .02, .5, .01), betasPsiA = .4, betasPsiB = .6, deltapsi = c(.2, -.2, .2, -.2))
other <- fixture[c(4, 2, 1), ];other$deltapsi <- c(.2, -.3, .01)
j <- join_common_events(fixture, other)
stopifnot(identical(j$EVENT, c('EX1', 'INT1', 'ALT1-2/2')), nrow(j) == 3,
  identical(j$comparison, c('Spire only', 'Both: same direction', 'Both: opposite direction')),
  j$locus[3] == 'ALT1')
other$GENE[1] <- 'wrong'
stopifnot(inherits(try(join_common_events(fixture, other), silent = TRUE), 'try-error'))
cat('Comparison universe, overlap probability, ID joining, concordance and locus checks passed.\n')
