feature_paths <- vapply(studies, function(p) file.path(p, 'tables', 'Matt_all_feature_Wilcoxon_statistics.csv'), character(1))
stopifnot(all(file.exists(feature_paths)))
feature_results <- lapply(feature_paths, function(p) as.data.frame(fread(p)))
joined_features <- merge(feature_results$spire, feature_results$fmn2,
  by = c('event_type', 'feature', 'group1', 'group2'), suffixes = c('_spire', '_fmn2'))
joined_features$same_effect_direction <- joined_features$difference_spire * joined_features$difference_fmn2 > 0
joined_features$raw_significant_both <- joined_features$p.value_spire <= .05 & joined_features$p.value_fmn2 <= .05
joined_features$BH_significant_both <- joined_features$p.BH_spire <= .05 & joined_features$p.BH_fmn2 <= .05
write_table(joined_features, 'Matt_feature_comparison_all_tests')
write_table(joined_features[which(joined_features$same_effect_direction & joined_features$raw_significant_both), ],
  'Matt_shared_feature_changes_raw_p')
write_table(joined_features[which(joined_features$same_effect_direction & joined_features$BH_significant_both), ],
  'Matt_shared_feature_changes_BH')
main_features <- c('UPINTRON_MEDIANLENGTH', 'DOINTRON_MEDIANLENGTH', 'MEDIAN_EXON_NUMBER',
  'INTRON_LENGTH', 'INTRON_GCC', 'MAXENTSCR_HSAMODEL_5SS')
rows <- joined_features[joined_features$group2 == 'ndiff' & (joined_features$feature %in% main_features |
  joined_features$raw_significant_both %in% TRUE), ]
if (nrow(rows)) {
  x <- bind_rows(lapply(c('spire', 'fmn2'), function(ds) {
    data.frame(dataset = ds, event_type = rows$event_type, feature = rows$feature, group = rows$group1,
      effect = rows[[paste0('difference_', ds)]] / rows[[paste0('background_IQR_', ds)]],
      p.value = rows[[paste0('p.value_', ds)]], p.BH = rows[[paste0('p.BH_', ds)]])
  }))
  x$effect[!is.finite(x$effect)] <- NA_real_
  x$label <- paste(x$event_type, x$feature, x$group, sep = ': ')
  write_table(x, 'Matt_comparison_plotted_effects')
  save_plot(ggplot(x, aes(dataset, label)) + geom_point(aes(color = effect, size = -log10(pmax(p.value, 1e-300)), shape = p.BH <= .05)) +
    scale_color_gradient2(low = '#5B8FA3', mid = 'white', high = '#C85450', midpoint = 0) +
    labs(x = NULL, y = NULL, color = 'Median shift / ndiff IQR', size = '-log10 raw p', shape = 'BH p ≤0.05',
      title = 'Sequence-feature responses in both knockouts', subtitle = 'Each study uses its own tested-event background.\nShared feature signatures do not establish an RBP mechanism.') + theme(legend.box = 'vertical'),
    'Comparison_Matt_feature_signatures', 12, max(8, 2 + length(unique(x$label)) * .25))
}
if (all(vapply(studies, function(p) file.exists(file.path(p, 'tables', 'protein_impact_summary.csv')), logical(1)))) {
  x <- bind_rows(lapply(names(studies), function(ds) {
    z <- fread(file.path(studies[[ds]], 'tables', 'protein_impact_summary.csv'));z$dataset <- ds;z
  }))
  write_table(x, 'protein_impact_comparison')
  save_plot(ggplot(x, aes(dataset, 100 * prop, fill = Impact)) + geom_col(width = .65) + scale_fill_brewer(palette = 'Set2') +
    labs(x = NULL, y = '% of selected events', title = 'Predicted functional impact in each full dataset',
      subtitle = 'Descriptive full-dataset proportions; the tested event universes differ'), 'Comparison_protein_impact', 10, 5)
}
