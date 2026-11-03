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
