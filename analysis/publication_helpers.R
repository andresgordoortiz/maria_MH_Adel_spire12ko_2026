# Raw-PSI sensitivity is descriptive; it does not refit or recalibrate betAS.
psi_sensitivity <- function(a, b, fitted_delta) {
  stopifnot(is.matrix(a), is.matrix(b), nrow(a) == nrow(b),
    ncol(a) >= 2, ncol(b) >= 2, length(fitted_delta) == nrow(a),
    all(is.finite(c(a, b, fitted_delta))), all(c(a, b) >= 0 & c(a, b) <= 100))
  delta <- rowMeans(b) - rowMeans(a)
  leave_one <- cbind(vapply(seq_len(ncol(a)), function(i)
    rowMeans(b) - rowMeans(a[, -i, drop = FALSE]), numeric(nrow(a))),
    vapply(seq_len(ncol(b)), function(i)
      rowMeans(b[, -i, drop = FALSE]) - rowMeans(a), numeric(nrow(a))))
  signed <- leave_one * sign(fitted_delta)
  data.frame(raw_delta_pp = delta,
    raw_agrees_with_fit = delta * fitted_delta > 0,
    loo_same_direction = rowSums(signed > 0) == ncol(signed),
    loo_min_signed_pp = apply(signed, 1, min),
    loo_min_pp = apply(leave_one, 1, min), loo_max_pp = apply(leave_one, 1, max),
    complete_separation = ifelse(fitted_delta > 0,
      apply(b, 1, min) > apply(a, 1, max), apply(b, 1, max) < apply(a, 1, min)))
}

degree_matched_overlap <- function(spire_selected, fmn2_selected, tested_events, iterations = 10000L) {
  stopifnot(is.logical(spire_selected), is.logical(fmn2_selected),
    length(spire_selected) == length(fmn2_selected), length(tested_events) == length(spire_selected),
    !anyNA(c(spire_selected, fmn2_selected)), all(is.finite(tested_events) & tested_events >= 1),
    iterations >= 1)
  strata <- split(seq_along(tested_events), floor(log2(tested_events)))
  set.seed(42)
  null <- replicate(iterations, sum(vapply(strata, function(ix) {
    n <- sum(fmn2_selected[ix])
    sum(spire_selected[ix[sample.int(length(ix), n)]])
  }, integer(1))))
  observed <- sum(spire_selected & fmn2_selected)
  list(observed = observed, null = null, expected = mean(null),
    pvalue = (1 + sum(null >= observed)) / (iterations + 1))
}

sequence_gc <- function(x) {
  stopifnot(is.character(x), !anyNA(x))
  clean <- gsub('[^ACGT]', '', toupper(x))
  stopifnot(all(nchar(clean) > 0))
  (nchar(clean) - nchar(gsub('[GC]', '', clean))) / nchar(clean)
}

go_counts <- function(universe, query, sets) {
  stopifnot(!anyDuplicated(universe), !anyDuplicated(query), all(query %in% universe),
    length(universe) > length(query), length(query) > 0, !anyDuplicated(names(sets)))
  sets <- lapply(sets, intersect, universe)
  n <- lengths(sets); k <- lengths(lapply(sets, intersect, query))
  data.frame(ID = names(sets), universe_genes = length(universe), query_genes = length(query),
    background_term_genes = n, Count = k, expected = length(query)*n/length(universe),
    FoldEnrichment = (k/length(query))/(n/length(universe)),
    pvalue = phyper(k-1, n, length(universe)-n, length(query), lower.tail = FALSE))
}

matched_category_counts <- function(selected, membership, degree, iterations = 10000L) {
  stopifnot(is.logical(selected), is.matrix(membership), is.logical(membership),
    nrow(membership) == length(selected), length(degree) == length(selected),
    !anyNA(c(selected, membership, degree)), all(degree >= 1), iterations > 0)
  strata <- split(seq_along(degree), floor(log2(degree)))
  set.seed(42)
  null <- matrix(replicate(iterations, {
    sampled <- unlist(lapply(strata, function(ix) ix[sample.int(length(ix), sum(selected[ix]))]))
    colSums(membership[sampled, , drop = FALSE])
  }), nrow = iterations, ncol = ncol(membership), byrow = TRUE)
  observed <- colSums(membership[selected, , drop = FALSE])
  list(summary = data.frame(ID = colnames(membership), observed = observed, expected_matched = colMeans(null),
    pvalue = (1 + colSums(sweep(null, 2, observed, '>=')))/(iterations+1)), null = null)
}

# Numeric leaf positions let dendrograms and heatmaps share exactly the same y axis.
row_clustering <- function(m) {
  stopifnot(is.matrix(m), nrow(m) >= 2, all(is.finite(m)), !anyDuplicated(rownames(m)))
  hc <- hclust(dist(m), method = 'complete')
  n <- nrow(m)
  y <- numeric(n - 1L); segments <- list()
  leaf_y <- numeric(n); leaf_y[hc$order] <- rev(seq_len(n))
  for (i in seq_len(n - 1L)) {
    children <- hc$merge[i, ]
    cy <- vapply(children, function(j) if (j < 0) leaf_y[-j] else y[j], numeric(1))
    cx <- vapply(children, function(j) if (j < 0) 0 else hc$height[j], numeric(1))
    y[i] <- mean(cy)
    segments[[i]] <- data.frame(x = c(cx, hc$height[i]), xend = rep(hc$height[i], 3),
      y = c(cy, cy[1]), yend = c(cy, cy[2]))
  }
  list(order = hc$order, y = leaf_y, segments = do.call(rbind, segments))
}

publication_theme <- function() {
  ggplot2::theme_bw(base_size = 8, base_family = 'Liberation Sans') +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = '#EBEBEB', linewidth = .25),
      panel.border = ggplot2::element_rect(color = '#777777', linewidth = .35),
      plot.title = ggplot2::element_text(size = 9, face = 'bold', margin = ggplot2::margin(b = 7)),
      plot.subtitle = ggplot2::element_text(size = 7, color = '#555555', margin = ggplot2::margin(b = 5)),
      plot.tag = ggplot2::element_text(size = 12, face = 'bold'),
      axis.text = ggplot2::element_text(size = 7, color = '#222222'),
      axis.title = ggplot2::element_text(size = 8, face = 'bold'),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_line(linewidth = .3),
      strip.background = ggplot2::element_rect(fill = '#F5F5F5', color = NA),
      strip.text = ggplot2::element_text(size = 7, face = 'bold'),
      legend.position = 'bottom', legend.title = ggplot2::element_text(size = 7),
      legend.text = ggplot2::element_text(size = 7),
      legend.key.height = grid::unit(3, 'mm'), legend.key.width = grid::unit(4, 'mm'),
      plot.margin = ggplot2::margin(6, 7, 6, 7))
}
