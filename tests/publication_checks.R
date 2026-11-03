source('analysis/publication_helpers.R')
a <- rbind(c(10, 20, 30), c(10, 20, 30), c(20, 20, 20), c(80, 90, 100))
b <- rbind(c(60, 70, 80), c(0, 0, 100), c(20, 20, 20), c(10, 20, 30))
x <- psi_sensitivity(a, b, c(.5, .1, 0, -.7))
stopifnot(identical(x$complete_separation, c(TRUE, FALSE, FALSE, TRUE)),
  identical(x$loo_same_direction, c(TRUE, FALSE, FALSE, TRUE)),
  all.equal(x$raw_delta_pp, c(50, 100/3 - 20, 0, -70)),
  identical(x$loo_min_signed_pp, c(45, -20, 0, 65)))
# Swapping groups and fitted signs preserves sensitivity, but reverses effects.
y <- psi_sensitivity(b, a, c(-.5, -.1, 0, .7))
stopifnot(all.equal(y$raw_delta_pp, -x$raw_delta_pp),
  identical(y$loo_same_direction, x$loo_same_direction),
  identical(y$loo_min_signed_pp, x$loo_min_signed_pp))
stopifnot(inherits(try(psi_sensitivity(a, b * NA, rep(.1, 4)), silent = TRUE), 'try-error'))
z <- degree_matched_overlap(c(TRUE, FALSE, TRUE), c(TRUE, TRUE, FALSE), c(1, 1, 4), 20L)
stopifnot(z$observed == 1L, all(z$null == 1), z$pvalue == 1)
z <- degree_matched_overlap(rep(FALSE, 3), rep(FALSE, 3), c(1, 1, 4), 20L)
stopifnot(all(z$null == 0), z$pvalue == 1)
cat('PASS: leave-one-sample sensitivity, sign reversal, ties and invalid PSI.\n')
m <- rbind(a = c(0, 1), b = c(0, 1.1), c = c(9, 10), d = c(9, 10.1))
tree <- row_clustering(m)
stopifnot(identical(sort(tree$order), seq_len(nrow(m))),
  abs(tree$y[1] - tree$y[2]) == 1, abs(tree$y[3] - tree$y[4]) == 1,
  nrow(tree$segments) == 3*(nrow(m)-1), all(tree$segments$xend >= tree$segments$x),
  all(tree$segments$y >= 1 & tree$segments$y <= nrow(m)))
stopifnot(inherits(try(row_clustering(m * NA), silent = TRUE), 'try-error'))
cat('PASS: row-cluster adjacency, dendrogram geometry and invalid matrices.\n')
stopifnot(identical(sequence_gc(c('ACGT', 'ggNNcc', 'AATT')), c(.5, 1, 0)),
  inherits(try(sequence_gc('NNN'), silent = TRUE), 'try-error'))
cat('PASS: sequence GC fractions, ambiguous bases and empty sequence rejection.\n')
