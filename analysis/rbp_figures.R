# Paper Figure 2H / S6: minimum motif p-value per RBP and region.
checks <- read.csv(file.path(rmaps_dir, 'input_checksums.csv'))
expected <- read.csv(file.path(out, 'external_inputs', 'input_checksums.csv'))
stopifnot(identical(checks, expected))
rbp_all <- list()
for (kind in c('EX', 'IR', 'Alt5', 'Alt3')) {
  region_count <- if (kind == 'EX') 8L else if (kind == 'IR') 4L else 6L
  for (direction in c('up', 'down')) {
    name <- if (direction == 'up') 'pVal.up.vs.bg.RNAmap.txt' else 'pVal.dn.vs.bg.RNAmap.txt'
    path <- file.path(rmaps_dir, kind, name)
    if (!file.exists(path)) next
    df <- as.data.frame(fread(path))
    stopifnot(ncol(df) == region_count + 1L)
    names(df) <- c('RBP', paste0('R', seq_len(region_count)))
    for (column in names(df)[-1]) if (is.logical(df[[column]]) && all(is.na(df[[column]]))) df[[column]] <- as.numeric(df[[column]])
    stopifnot(!anyNA(df$RBP), all(vapply(df[-1], is.numeric, logical(1))))
    df <- pivot_longer(df, -RBP, names_to = 'region', values_to = 'pvalue')
    if (any(!is.na(df$pvalue) & (!is.finite(df$pvalue) | df$pvalue < 0 | df$pvalue > 1))) stop('Invalid rMAPS2 p-value')
    df$RBP <- sub('[.].*', '', df$RBP)
    df <- df %>% filter(!is.na(pvalue)) %>% group_by(RBP, region) %>% summarise(pvalue = min(pvalue), .groups = 'drop')
    df$type <- kind; df$direction <- direction
    rbp_all[[length(rbp_all) + 1L]] <- df
    source_files <- c(source_files, path)
  }
}
if (!length(rbp_all)) stop('RMAPS_OUTPUT_DIR exists but contains no supported output tables')
rbp_all <- bind_rows(rbp_all)
write_table(rbp_all, 'RBP_minimum_motif_pvalues')
for (kind in unique(rbp_all$type)) {
  df <- filter(rbp_all, type == kind)
  selected_rbps <- df %>% group_by(RBP) %>% summarise(minp = min(pvalue), .groups = 'drop') %>%
    filter(minp <= .05) %>% arrange(minp) %>% slice_head(n = 50) %>% pull(RBP)
  if (!length(selected_rbps)) {
    save_plot(empty_plot(paste('No enriched RBPs:', kind)), paste0('Fig2H_S6_RBP_', kind))
    next
  }
  df <- df[df$RBP %in% selected_rbps, ]
  positive <- df$pvalue[df$pvalue > 0]
  floorp <- if (length(positive)) min(positive) / 10 else 1e-300
  df$score <- -log10(pmax(df$pvalue, max(floorp, 1e-300)))
  wide <- df %>% mutate(column = paste(region, direction, sep = '_')) %>% select(RBP, column, score) %>%
    pivot_wider(names_from = column, values_from = score, values_fill = 0)
  stopifnot(all(is.finite(as.matrix(wide[-1]))))
  order_rbps <- if (nrow(wide) > 1) wide$RBP[hclust(dist(as.matrix(wide[-1])))$order] else wide$RBP
  df$RBP <- factor(df$RBP, levels = order_rbps)
  df$column <- factor(paste(df$region, df$direction, sep = '_'),
    levels = as.vector(t(outer(sort(unique(df$region)), c('up', 'down'), paste, sep = '_'))))
  save_plot(ggplot(df, aes(column, RBP, fill = score)) + geom_tile(color = 'white', linewidth = .1) +
    scale_fill_viridis_c() + labs(x = 'Genomic region × direction', y = NULL, fill = '-log10 p',
      title = paste(kind, 'RBP motif enrichment'), subtitle = 'Minimum motif p per RBP/region; up to 50 RBPs; raw p ≤0.05 in any region') +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)), paste0('Fig2H_S6_RBP_', kind), 10, max(5, length(selected_rbps) * .16))
  counts <- df %>% filter(pvalue <= .05) %>% count(region, direction, name = 'n_enriched_RBPs')
  write_table(counts, paste0('RBP_', kind, '_region_counts'))
  save_plot(ggplot(counts, aes(region, n_enriched_RBPs, fill = direction)) +
    geom_col(position = 'dodge') + scale_fill_manual(values = direction_colors) +
    labs(x = 'Region', y = 'Enriched RBPs', title = paste(kind, 'motif enrichment by region')),
    paste0('Supplement_RBP_', kind, '_region_counts'), 7, 4)
}
status <- rbind(status, data.frame(step = 'RBP_heatmaps_and_region_summaries', status = 'complete from supplied rMAPS2 outputs; raw minimum p-values'))
