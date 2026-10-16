message('Sequence features, protein impact and GO')
annotation <- as.data.frame(fread(annotation_path, na.strings = c('', 'NA')))
stopifnot(all(c('EVENT', 'GENE', 'COORD_o') %in% names(annotation)), !anyDuplicated(annotation$EVENT))
matched <- events$EVENT %in% annotation$EVENT
annotation <- annotation[match(events$EVENT, annotation$EVENT), ]
stopifnot(identical(annotation$EVENT[matched], events$EVENT[matched]))
changed_coords <- matched & !is.na(annotation$COORD_o) & annotation$COORD_o != events$COORD
write_table(data.frame(EVENT = events$EVENT[changed_coords], input_coord = events$COORD[changed_coords],
  reference_coord = annotation$COORD_o[changed_coords]), 'VASTDB_changed_event_coordinates')
annotation[changed_coords, setdiff(names(annotation), c('EVENT', 'GENE'))] <- NA
renamed <- matched & !is.na(annotation$GENE) & annotation$GENE != events$GENE
write_table(data.frame(EVENT = events$EVENT[renamed], input_gene = events$GENE[renamed],
  reference_gene = annotation$GENE[renamed], direction = events$direction[renamed]), 'VASTDB_gene_name_changes')
write_table(events[!matched, c('EVENT', 'GENE', 'COMPLEX', 'direction')], 'missing_VASTDB_event_annotations')
annotation$EVENT <- events$EVENT
annotation$reference_GENE <- annotation$GENE
annotation$GENE <- events$GENE
source_files <- c(source_files, annotation_path, 'data/reference/PROT_IMPACT-mm10-v3.tab.gz')
features <- cbind(events[, c('EVENT', 'GENE', 'event_type', 'direction', 'deltapsi', 'FDR')],
  annotation[, c('TARGET_LENGTH', 'TARGET_GC', 'C1_GC', 'C2_GC')])
write_table(features, 'VASTDB_sequence_features')
feature_stats <- list()
for (kind in c('EX', 'IR')) for (feature in c('TARGET_LENGTH', 'TARGET_GC', 'C1_GC', 'C2_GC')) {
  df <- features[features$event_type == kind, ]
  set.seed(42)
  bg <- which(df$direction == 'ndiff')
  if (length(bg) > 10000) df <- df[c(which(df$direction != 'ndiff'), sample(bg, 10000)), ]
  for (pair in list(c('up', 'ndiff'), c('down', 'ndiff'), c('up', 'down'))) {
    a <- df[[feature]][df$direction == pair[1]]; a <- a[is.finite(a)]
    b <- df[[feature]][df$direction == pair[2]]; b <- b[is.finite(b)]
    p <- if (min(length(a), length(b)) >= 2 && length(a) + length(b) >= 10)
      suppressWarnings(wilcox.test(a, b, exact = FALSE)$p.value) else NA_real_
    feature_stats[[length(feature_stats) + 1L]] <- data.frame(event_type = kind, feature = feature,
      group1 = pair[1], group2 = pair[2], n1 = length(a), n2 = length(b),
      median1 = if (length(a)) median(a) else NA_real_, median2 = if (length(b)) median(b) else NA_real_, p.value = p)
  }
  plot <- ggplot(df, aes(direction, .data[[feature]], fill = direction)) + geom_boxplot(outlier.shape = NA, width = .6) +
    stat_summary(fun = mean, geom = 'point', shape = 18, size = 2) +
    scale_fill_manual(values = direction_colors) + labs(x = NULL, y = if (grepl('GC', feature)) 'GC (%)' else 'Length (nt)',
      title = paste(kind, gsub('_', ' ', feature)), subtitle = 'Wilcoxon tests use all values; ndiff background capped at 10,000 (seed 42)')
  if (feature == 'TARGET_LENGTH') plot <- plot + scale_y_log10()
  save_plot(plot, paste0('Supplement_VAST_', kind, '_', feature), 7, 4.5)
}
feature_stats <- bind_rows(feature_stats)
feature_stats$p.BH <- p.adjust(feature_stats$p.value, 'BH')
write_table(feature_stats, 'VASTDB_feature_Wilcoxon_statistics')
status <- rbind(status, data.frame(step = 'VASTDB_target_length_GC', status = 'complete; not a replacement for Matt features'))

impact <- as.data.frame(fread(cmd = 'gzip -dc data/reference/PROT_IMPACT-mm10-v3.tab.gz'))
stopifnot(all(c('EventID', 'ONTO') %in% names(impact)), !anyDuplicated(impact$EventID))
onto <- impact$ONTO[match(events$EVENT, impact$EventID)]
simplified <- rep('Unknown/noncoding', length(onto))
onto[is.na(onto)] <- ''
simplified[grepl('inclusion', onto, ignore.case = TRUE)] <- 'ORF disrupted upon inclusion'
simplified[grepl('exclusion', onto, ignore.case = TRUE)] <- 'ORF disrupted upon exclusion'
simplified[grepl('alternative|isoform', onto, ignore.case = TRUE)] <- 'Alternative isoform'
simplified[grepl('UTR', onto, ignore.case = TRUE)] <- 'Regulatory'
impact_levels <- c('ORF disrupted upon inclusion', 'ORF disrupted upon exclusion', 'Alternative isoform', 'Regulatory', 'Unknown/noncoding')
impact_events <- data.frame(EVENT = events$EVENT, GENE = events$GENE, direction = events$direction,
  deltapsi = events$deltapsi, FDR = events$FDR, ONTO = onto, Impact = factor(simplified, levels = impact_levels),
  predicted_disruption_in_DKO = (simplified == 'ORF disrupted upon inclusion' & events$direction == 'up') |
    (simplified == 'ORF disrupted upon exclusion' & events$direction == 'down'))
write_table(impact_events, 'protein_impact_all_events')
impact_selected <- impact_events[events$direction != 'ndiff', ]
write_table(impact_selected, 'protein_impact_differential_events')
impact_summary <- impact_selected %>% count(Impact) %>% mutate(prop = n / sum(n))
write_table(impact_summary, 'protein_impact_summary')
write_table(data.frame(tested = nrow(events), matched = sum(nzchar(onto)),
  selected = nrow(impact_selected), selected_matched = sum(nzchar(impact_selected$ONTO))), 'protein_impact_annotation_coverage')
if (nrow(impact_summary)) {
  cells <- floor(impact_summary$prop * 100)
  remainder <- 100 - sum(cells)
  if (remainder > 0) cells[order(impact_summary$prop * 100 - cells, decreasing = TRUE)[seq_len(remainder)]] <-
    cells[order(impact_summary$prop * 100 - cells, decreasing = TRUE)[seq_len(remainder)]] + 1L
  stopifnot(sum(cells) == 100)
  waffle <- data.frame(x = rep(1:10, 10), y = rep(10:1, each = 10), Impact = rep(impact_summary$Impact, cells))
  save_plot(ggplot(waffle, aes(x, y, fill = Impact)) + geom_tile(color = 'white', linewidth = .4) + coord_equal() +
    scale_fill_brewer(palette = 'Set2') + theme_void() + labs(title = 'Predicted functional impact of differential events',
      subtitle = 'VASTDB PROT_IMPACT mm10 v3; each square ≈1%'), 'Fig3A_protein_impact', 8, 5)
  save_plot(ggplot(impact_selected %>% count(Impact, direction), aes(Impact, n, fill = direction)) + geom_col() +
    coord_flip() + scale_fill_manual(values = direction_colors) + labs(x = NULL, y = 'Events', title = 'Functional annotations and splicing direction'),
    'Supplement_protein_impact_direction', 8, 5)
}
status <- rbind(status, data.frame(step = 'protein_impact', status = 'complete; available VASTDB v3; original v2.3 not archived in paper repo'))

if (requireNamespace('clusterProfiler', quietly = TRUE) && requireNamespace('org.Mm.eg.db', quietly = TRUE)) {
  suppressPackageStartupMessages(library(org.Mm.eg.db))
  gene_map <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db,
    keys = unique(events$GENE), keytype = 'SYMBOL', columns = 'ENTREZID'))
  write_table(gene_map, 'GO_gene_mapping')
  map_genes <- function(g) unique(na.omit(gene_map$ENTREZID[gene_map$SYMBOL %in% g]))
  sets <- list(all = unique(impact_selected$GENE),
    predicted_ORF_disruption = unique(impact_selected$GENE[impact_selected$predicted_disruption_in_DKO]),
    EX = unique(events$GENE[events$direction != 'ndiff' & events$event_type == 'EX']),
    IR = unique(events$GENE[events$direction != 'ndiff' & events$event_type == 'IR']))
  go_results <- list()
  for (ont in c('BP', 'CC', 'MF')) for (background in c('paper_OrgDb', 'tested_N10')) {
    for (set in names(sets)) {
      if (set != 'all' && ont != 'BP') next
      message('GO ', background, ' ', set, ' ', ont)
      ids <- map_genes(sets[[set]])
      universe <- if (background == 'tested_N10') map_genes(events$GENE) else NULL
      ego <- if (length(ids) >= 5) suppressMessages(clusterProfiler::enrichGO(gene = ids,
        OrgDb = org.Mm.eg.db, keyType = 'ENTREZID', ont = ont, universe = universe,
        pAdjustMethod = 'BH', pvalueCutoff = .05, qvalueCutoff = .2,
        minGSSize = 10, maxGSSize = 500, readable = TRUE)) else NULL
      df <- if (is.null(ego)) data.frame() else as.data.frame(ego)
      name <- paste('GO', background, set, ont, sep = '_')
      write_table(df, name)
      saveRDS(ego, file.path(out, 'tables', paste0(name, '.rds')))
      go_results[[name]] <- data.frame(background = background, gene_set = set, ontology = ont,
        input_genes = length(sets[[set]]), mapped_genes = length(ids), enriched_terms = nrow(df))
      if (set == 'all' || set == 'predicted_ORF_disruption') {
        if (nrow(df)) {
          df <- head(df[order(-df$Count, df$p.adjust), ], 20)
          df$ratio <- vapply(strsplit(df$GeneRatio, '/', fixed = TRUE), function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1))
          df$Description <- factor(df$Description, levels = rev(df$Description))
          p <- ggplot(df, aes(ratio, Description)) + geom_segment(aes(x = 0, xend = ratio, yend = Description), color = 'grey70') +
            geom_point(aes(size = Count, color = p.adjust)) + scale_color_gradient(low = '#d73027', high = '#4575b4', trans = 'log10') +
            labs(x = 'Gene ratio', y = NULL, title = paste('Spire GO', ont, '-', set), subtitle = paste('Universe:', background))
        } else p <- empty_plot(paste('No enriched GO terms:', set, ont, background))
        save_plot(p, name, 10, 7)
      }
    }
  }
  write_table(bind_rows(go_results), 'GO_analysis_summary')
  status <- rbind(status, data.frame(step = 'GO_BP_CC_MF', status = 'complete: paper OrgDb universe + tested-N10 sensitivity analysis'))
} else status <- rbind(status, data.frame(step = 'GO_BP_CC_MF', status = 'pending: scripts/run_R.sh analysis/setup.R'))
