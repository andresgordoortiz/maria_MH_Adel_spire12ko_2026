# Paper Figure 2G, feature boxplots and the repository's GC gradient analysis.
matt_stats <- list()
matt_tables <- list()
for (kind in c('EX', 'IR')) {
  suffix <- if (kind == 'EX') 'efeatures' else 'ifeatures'
  path <- file.path(matt_dir, kind, paste0('spire_', kind, '_matt_with_', suffix, '.tab'))
  df <- as.data.frame(fread(path, na.strings = c('NA', '', 'nan', '-nan', 'NaN')))
  source_files <- c(source_files, path, file.path(matt_dir, 'COMPLETE'), file.path(matt_dir, 'input_checksums.csv'))
  stopifnot(nrow(df) > 0, all(c('EVENT', 'DATASET') %in% names(df)), !anyDuplicated(df$EVENT),
    all(df$EVENT %in% events$EVENT), all(as.character(df$DATASET) == as.character(events$direction[match(df$EVENT, events$EVENT)])))
  expected <- as.data.frame(fread(file.path(out, 'external_inputs', paste0('spire_', kind, '_matt.tab'))))
  stopifnot(all(df$EVENT %in% expected$EVENT))
  for (column in names(expected)) stopifnot(isTRUE(all.equal(df[[column]],
    expected[[column]][match(df$EVENT, expected$EVENT)], check.attributes = FALSE)))
  matt_tables[[kind]] <- df
  write_table(df %>% count(DATASET, name = 'retained_after_Matt_GTF_and_duplicate_filters'), paste0('Matt_', kind, '_retained_counts'))
  feature_cols <- setdiff(names(df)[vapply(df, is.numeric, logical(1))],
    c('START', 'END', 'UPSTRM_EX_BORDER', 'DOSTRM_EX_BORDER'))
  for (feature in feature_cols) for (pair in list(c('up', 'ndiff'), c('down', 'ndiff'), c('up', 'down'))) {
    a <- df[[feature]][df$DATASET == pair[1]]; a <- a[is.finite(a)]
    b <- df[[feature]][df$DATASET == pair[2]]; b <- b[is.finite(b)]
    p <- if (length(a) >= 3 && length(b) >= 3 && length(a) + length(b) >= 10)
      suppressWarnings(wilcox.test(a, b, exact = FALSE)$p.value) else NA_real_
    matt_stats[[length(matt_stats) + 1L]] <- data.frame(event_type = kind, feature = feature,
      group1 = pair[1], group2 = pair[2], n1 = length(a), n2 = length(b),
      median1 = if (length(a)) median(a) else NA_real_, median2 = if (length(b)) median(b) else NA_real_,
      p.value = p, difference = if (length(a) && length(b)) median(a) - median(b) else NA_real_)
  }
}
matt_stats <- bind_rows(matt_stats)
matt_stats$p.BH <- p.adjust(matt_stats$p.value, 'BH')
write_table(matt_stats, 'Matt_all_feature_Wilcoxon_statistics')
keep <- list(EX = c('UPINTRON_MEDIANLENGTH', 'DOINTRON_MEDIANLENGTH', 'MEDIAN_EXON_NUMBER'),
  IR = c('INTRON_LENGTH', 'INTRON_GCC', 'MAXENTSCR_HSAMODEL_5SS'))
stars <- function(p) ifelse(is.na(p), 'NA', ifelse(p <= .0001, '****', ifelse(p <= .001, '***',
  ifelse(p <= .01, '**', ifelse(p <= .05, '*', 'ns')))))
for (kind in c('EX', 'IR')) {
  summary <- matt_stats %>% filter(event_type == kind, feature %in% keep[[kind]], group2 == 'ndiff') %>%
    mutate(significance = stars(p.value), effect = ifelse(difference >= 0, 'Increase', 'Decrease'))
  save_plot(ggplot(summary, aes(group1, feature)) + geom_point(aes(color = effect, size = -log10(pmax(p.value, 1e-300))), alpha = .8) +
    geom_text(aes(label = significance), nudge_x = .2, size = 3) +
    scale_color_manual(values = c(Increase = '#C85450', Decrease = '#5B8FA3')) +
    labs(x = 'Versus ndiff', y = NULL, title = paste(kind, 'sequence-feature summary'),
      subtitle = 'Paper features; stars use raw Mann–Whitney p-values', size = '-log10 p'), paste0('Fig2G_Matt_', kind), 9, 4)
  sig_features <- unique(matt_stats$feature[matt_stats$event_type == kind & !is.na(matt_stats$p.value) & matt_stats$p.value <= .05])
  selected_features <- union(keep[[kind]], sig_features)
  df <- matt_tables[[kind]]
  for (feature in selected_features) {
    stopifnot(feature %in% names(df))
    plot_df <- data.frame(direction = factor(df$DATASET, levels = c('up', 'down', 'ndiff')), value = df[[feature]])
    plot_df <- plot_df[is.finite(plot_df$value), ]
    if (!nrow(plot_df)) next
    bounds <- quantile(plot_df$value, c(.25, .75))
    view <- c(bounds[1] - 3 * diff(bounds), bounds[2] + 3 * diff(bounds))
    tests <- matt_stats[matt_stats$event_type == kind & matt_stats$feature == feature, ]
    subtitle <- paste(paste(tests$group1, 'vs', tests$group2, ':', stars(tests$p.value)), collapse = '; ')
    p <- ggplot(plot_df, aes(direction, value, fill = direction)) + geom_boxplot(outlier.shape = NA, width = .6) +
      stat_summary(fun = mean, geom = 'point', shape = 18, size = 2) + scale_fill_manual(values = direction_colors) +
      labs(x = NULL, y = gsub('_', ' ', feature), title = paste(kind, gsub('_', ' ', feature)), subtitle = subtitle)
    if (diff(view) > 0) p <- p + coord_cartesian(ylim = view)
    save_plot(p, paste0('Supplement_Matt_', kind, '_', feature), 8, 4.5)
  }
}
source('analysis/gc_gradient.R')
