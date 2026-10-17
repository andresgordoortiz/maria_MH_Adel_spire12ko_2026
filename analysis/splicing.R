# Paper Figure 2C/D/F and supplementary splicing volcanoes/heatmaps.
message('Core splicing plots')
counts <- events %>% count(event_type, direction, .drop = FALSE)
write_table(counts, 'event_counts')
write_table(events %>% count(COMPLEX, direction), 'event_counts_by_complexity')
selected <- events$direction != 'ndiff'
write_table(data.frame(tested_without_gene = sum(!events$GENE %in% valid_genes(events$GENE)),
  selected_without_gene = sum(selected & !events$GENE %in% valid_genes(events$GENE))), 'gene_annotation_coverage')
write_table(events %>% filter(selected, GENE %in% valid_genes(GENE)) %>% group_by(GENE) %>% summarise(
  n_events = n(), types = paste(sort(unique(event_type)), collapse = ';'),
  max_abs_deltapsi = max(abs(deltapsi)), .groups = 'drop'), 'differential_genes')

v <- events
positive <- v$FDR[v$FDR > 0]
floor_fpr <- if (length(positive)) min(positive) / 10 else 1e-6
v$neglog <- -log10(pmax(v$FDR, floor_fpr))
ceiling_fpr <- max(3.01, max(v$neglog))
v$display_y <- ifelse(v$neglog > 3, 3 + .5 * (v$neglog - 3) / (ceiling_fpr - 3), v$neglog)
set.seed(42)
v$display_y[v$FDR == 0] <- v$display_y[v$FDR == 0] + runif(sum(v$FDR == 0), -.04, .04)
v <- v[order(v$direction != 'ndiff'), ]
p <- ggplot(v, aes(deltapsi, display_y)) +
  geom_point(aes(color = direction), size = .4, alpha = .45) +
  scale_color_manual(values = direction_colors) +
  geom_vline(xintercept = c(-.1, .1), linetype = 'dashed', linewidth = .3) +
  geom_hline(yintercept = -log10(.05), linetype = 'dashed', linewidth = .3) +
  scale_y_continuous(breaks = c(0:3, 3.5), labels = c('0', '1', '2', '3', paste0('≥', round(ceiling_fpr, 1)))) +
  labs(x = paste0('ΔPSI (', contrast_label, ')'), y = '−log10 empirical FPR (compressed above 3)',
    title = paste(dataset_title, 'differential splicing'), color = 'Direction')
labels <- v %>% filter(direction != 'ndiff') %>% arrange(desc(abs(deltapsi))) %>% slice_head(n = 8)
p <- p + ggrepel::geom_text_repel(data = labels, aes(label = GENE), size = 3, seed = 42, max.overlaps = Inf)
save_plot(p, 'Fig2C_splicing_volcano')
save_plot(p + facet_wrap(~event_type, ncol = 2), 'Supplement_splicing_volcano_by_type', 9, 7)
write_table(data.frame(zero_FPR_floor_for_display = floor_fpr, upper_y_compression = 3,
  zero_FPR_events = sum(events$FDR == 0), zero_FPR_abs_delta_below_0_1 = sum(events$FDR == 0 & abs(events$deltapsi) < .1)), 'volcano_display_and_boundary_diagnostic')

save_plot(ggplot(filter(counts, direction != 'ndiff'), aes(event_type, n, fill = direction)) +
  geom_col(position = 'stack', width = .7) + geom_text(aes(label = n), position = position_stack(vjust = .5), size = 3) +
  scale_fill_manual(values = direction_colors) + labs(x = 'Event type', y = 'Selected events', fill = 'Direction',
    title = 'Differential events by type'), 'Fig2D_event_counts', 6, 4)
direction_counts <- counts %>% filter(direction != 'ndiff') %>% group_by(event_type) %>% mutate(percent = 100 * n / sum(n))
write_table(direction_counts, 'inclusion_skipping_percentages')
save_plot(ggplot(direction_counts, aes(event_type, percent, fill = direction)) +
  geom_col(width = .65) + geom_text(aes(label = sprintf('%.1f%%', percent)), position = position_stack(vjust = .5), size = 3) +
  scale_fill_manual(values = direction_colors) + labs(x = 'Event type', y = '% of selected events', fill = 'Direction',
    title = 'Inclusion/site use/retention balance'), 'Supplement_direction_balance', 6, 4)

for (kind in c('EX', 'IR')) {
  z <- events[events$event_type == kind & selected, ]
  if (!nrow(z)) {
    save_plot(empty_plot(paste('No selected', kind, 'events')), paste0('Fig2F_', kind, '_miniheatmap'))
    next
  }
  mat <- as.matrix(z[, samples]); rownames(mat) <- z$EVENT
  varying <- apply(mat, 1, sd) > 0
  write_table(z[!varying, c('EVENT', 'GENE', 'deltapsi', 'FDR')], paste0(kind, '_constant_PSI_excluded_from_clustering'))
  mat <- mat[varying, , drop = FALSE]
  if (!nrow(mat)) next
  mini <- robust_scale(mat, 1)
  cl <- cluster_rows(mini)
  collapsed <- data.frame(EVENT = rownames(mat), Cluster = cl$cluster,
    A = rowMeans(mini[, group_a_columns, drop = FALSE]), B = rowMeans(mini[, group_b_columns, drop = FALSE]))
  write_table(collapsed, paste0(kind, '_miniheatmap_event_clusters'))
  summary <- collapsed %>% group_by(Cluster) %>% summarise(n = n(), A = mean(A), B = mean(B), .groups = 'drop')
  write_table(summary, paste0(kind, '_miniheatmap_summary'))
  mini_plot <- pivot_longer(summary, c(A, B), names_to = 'Group', values_to = 'scaled_PSI')
  mini_plot$Group <- ifelse(mini_plot$Group == 'A', config$group_a, config$group_b)
  mini_plot$Pattern <- paste0('Pattern ', mini_plot$Cluster, ' (n=', mini_plot$n, ')')
  save_plot(ggplot(mini_plot, aes(factor(Group, levels = c(config$group_a, config$group_b)), Pattern, fill = scaled_PSI)) +
    geom_tile(width = .95, height = .95) + geom_text(aes(label = sprintf('%.2f', scaled_PSI)), size = 4) +
    scale_fill_gradientn(colors = signed_colors, limits = c(-1, 1)) + coord_equal() +
    labs(x = NULL, y = NULL, fill = 'MAD-scaled PSI', title = paste(kind, 'collapsed patterns')),
    paste0('Fig2F_', kind, '_miniheatmap'), 6, 3.5)
  full <- robust_scale(mat, 2)
  # The published large heatmap uses ComplexHeatmap's default Euclidean clustering.
  hc <- if (nrow(full) > 1) hclust(dist(full)) else NULL
  order_rows <- if (is.null(hc)) 1L else hc$order
  full_cluster <- if (is.null(hc)) 1L else cutree(hc, 2)
  # Cluster samples within each condition, as in the reference column_split.
  order_cols <- unlist(lapply(list(group_a_columns, group_b_columns), function(ix) ix[hclust(dist(t(full[, ix, drop = FALSE])))$order]))
  ordered <- full[order_rows, order_cols, drop = FALSE]
  hm <- as.data.frame(as.table(ordered)); names(hm) <- c('EVENT', 'sample', 'value')
  hm$row <- match(hm$EVENT, rownames(ordered))
  hm$sample <- factor(hm$sample, levels = colnames(ordered), labels = unname(sample_labels[colnames(ordered)]))
  write_table(data.frame(EVENT = rownames(full), Cluster = full_cluster, row_order = match(seq_len(nrow(full)), order_rows)),
    paste0(kind, '_full_heatmap_event_clusters'))
  write_table(data.frame(sample = samples, cor(full, method = 'spearman'), check.names = FALSE), paste0(kind, '_sample_correlations'))
  save_plot(ggplot(hm, aes(sample, row, fill = value)) + geom_raster() +
    scale_fill_gradientn(colors = signed_colors, limits = c(-2, 2)) +
    labs(x = NULL, y = paste0(nrow(full), ' events'), fill = 'MAD-scaled PSI', title = paste(kind, 'sample heatmap')) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)),
    paste0('Supplement_', kind, '_full_heatmap'), 6, 7)
}

sample_mat <- as.matrix(events[, samples])
cor_mat <- cor(sample_mat, method = 'spearman')
cor_df <- as.data.frame(as.table(cor_mat)); names(cor_df) <- c('SampleA', 'SampleB', 'rho')
cor_df$SampleA <- unname(sample_labels[as.character(cor_df$SampleA)])
cor_df$SampleB <- unname(sample_labels[as.character(cor_df$SampleB)])
write_table(data.frame(sample = samples, cor_mat, check.names = FALSE), 'QC_sample_Spearman')
save_plot(ggplot(cor_df, aes(SampleA, SampleB, fill = rho)) + geom_tile() +
  geom_text(aes(label = sprintf('%.3f', rho)), size = 3, color = 'white') +
  scale_fill_viridis_c(option = sequential_palette, direction = -1, limits = c(0, 1)) +
  coord_equal() + labs(x = NULL, y = NULL, title = 'N10 PSI sample correlations') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)), 'QC_sample_correlations', 6, 5)
