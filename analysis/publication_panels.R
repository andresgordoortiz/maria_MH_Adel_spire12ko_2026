# Main panels answer the biological question; three shared supplements hold supporting evidence.
polish <- function(p, title) p + labs(title = title, subtitle = NULL, caption = NULL)
heat_theme <- theme(panel.grid = element_blank(), panel.grid.major = element_blank(),
  panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
  axis.ticks = element_blank(), panel.border = element_blank())
usage_labels <- c(down = 'Decreased', ndiff = 'Background', up = 'Increased')
study_labels <- c(spire = 'SPIRE', fmn2 = 'FMN2')
tree_plot <- function(tree, n) ggplot(tree$segments) +
  geom_segment(aes(x, y, xend = xend, yend = yend), linewidth = .18, color = '#555555') +
  scale_x_reverse(expand = expansion(mult = c(.02, 0))) +
  scale_y_continuous(limits = c(.5, n + .5), expand = c(0, 0)) + theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

psi_heatmap <- function(ds) {
  x <- sensitivity[[ds]]
  m <- as.matrix(x[, names(groups[[ds]])]); rownames(m) <- x$EVENT
  m <- m - apply(m, 1, median)
  tree <- row_clustering(m)
  h <- sample_psi(x, ds); h$row <- tree$y[match(h$EVENT, x$EVENT)]
  h$centered_PSI <- as.vector(t(m))
  h$sample <- factor(h$sample, names(groups[[ds]]))
  h$group <- factor(h$group, unique(h$group), if (ds == 'spire') c('Control', 'Spire1/2 DKO') else c('Fmn2 Het', 'Fmn2 KO'))
  stopifnot(nrow(h) == length(m), all(h$PSI - rep(apply(as.matrix(x[, names(groups[[ds]])]), 1, median), each = ncol(m)) == h$centered_PSI))
  export(h, paste0(toupper(ds), '_clustered_PSI'), 'source_data')
  export(data.frame(EVENT = x$EVENT[tree$order], GENE = x$GENE[tree$order], event_type = x$event_type[tree$order]),
    paste0(toupper(ds), '_PSI_row_order'), 'source_data')
  p <- ggplot(h, aes(sample, row, fill = centered_PSI)) + geom_raster() +
    facet_grid(~group, scales = 'free_x', space = 'free_x') +
    scale_y_continuous(limits = c(.5, nrow(x)+.5), expand = c(0, 0)) +
    scale_x_discrete(labels = setNames(seq_along(groups[[ds]]), names(groups[[ds]]))) +
    scale_fill_gradient2(low = signed_colors[1], mid = signed_colors[2], high = signed_colors[3],
      limits = c(-50, 50), oob = scales::squish, breaks = c(-50, 0, 50), name = 'PSI − row median (pp)') +
    labs(x = 'Biological replicate', y = NULL, title = paste(format(nrow(x), big.mark = ','), 'splicing-event profiles')) +
    heat_theme + theme(axis.text.y = element_blank(), legend.position = 'bottom')
  wrap_elements(full = tree_plot(tree, nrow(x)) + p + plot_layout(widths = c(.13, 1)))
}
event_counts <- function(ds) {
  x <- sensitivity[[ds]] %>% count(event_type, direction)
  x$event_type <- factor(x$event_type, c('EX', 'IR', 'Alt5', 'Alt3'))
  export(x, paste0(toupper(ds), '_selected_event_counts'), 'source_data')
  ggplot(x, aes(event_type, n, fill = direction)) + geom_col(position = position_dodge(.8), width = .65) +
    geom_text(aes(label = n), position = position_dodge(.8), vjust = -.4) +
    scale_fill_manual(values = direction_colors, labels = usage_labels, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, .15))) +
    labs(x = 'Event class', y = 'Selected events', title = 'Changes across splicing classes')
}
ma <- dge[is.finite(dge$baseMean) & dge$baseMean > 0 & is.finite(dge$log2FoldChange), ]
ma <- ma[order(ma$selected), ]
export(ma, 'SPIRE_expression_MA', 'source_data')
p_ma <- ggplot(ma, aes(baseMean, log2FoldChange)) +
  geom_hline(yintercept = 0, color = '#AAAAAA', linewidth = .3) +
  geom_point(aes(color = direction), size = .65, alpha = .65) +
  ggrepel::geom_text_repel(data = ma[ma$gene_name %in% c('Spire1', 'Spire2', 'Mt1'), ],
    aes(label = gene_name), seed = 42, size = 2.5, fontface = 'italic', family = 'Liberation Sans',
    min.segment.length = 0, max.overlaps = Inf) +
  scale_x_log10(breaks = c(1, 100, 10000), labels = scales::label_number()) + direction_scale() +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(x = 'Mean normalized counts', y = expression(log[2]~'fold change'), title = '111 differentially expressed genes')

# The email asks about division-related splicing, so show all annotated events, without a top-gene ranking.
division_plot_data <- inner_join(division_selected[!division_selected$GENE %in% c('Spire1', 'Spire2'), ],
  unique(division[, c('SYMBOL', 'GOALL', 'term')]), by = c('GENE' = 'SYMBOL'), relationship = 'many-to-many')
division_plot_data$term <- factor(division_plot_data$term, rev(unname(division_terms)))
export(division_plot_data, 'SPIRE_division_event_annotations', 'source_data')
stopifnot(setequal(division_plot_data$EVENT, division_selected$EVENT))
category_counts <- division_plot_data %>% group_by(term) %>% summarise(events = n(), genes = n_distinct(GENE), .groups = 'drop')
labels <- setNames(paste0(category_counts$term, '\n', category_counts$events, ' events / ', category_counts$genes, ' genes'), category_counts$term)
p_division <- ggplot(division_plot_data, aes(100*deltapsi, term, color = event_type)) +
  geom_vline(xintercept = 0, color = '#999999', linewidth = .3) +
  geom_point(position = position_jitter(width = 0, height = .16, seed = 42), size = 1.35, alpha = .65) +
  scale_y_discrete(labels = labels) +
  scale_color_manual(values = c(EX = '#6A3D9A', IR = '#1B9E77', Alt5 = '#D95F02', Alt3 = '#7570B3'), name = 'Event class') +
  labs(x = 'Fitted ΔPSI (DKO − control; pp)', y = NULL,
    title = 'Splicing events in genes annotated to division processes')


# Six features fixed by the reference analysis, including nonsignificant results.
feature_spec <- data.frame(event_type = rep(c('EX', 'IR'), each = 3),
  feature = c('UPINTRON_MEDIANLENGTH', 'DOINTRON_MEDIANLENGTH', 'MEDIAN_EXON_NUMBER',
    'INTRON_LENGTH', 'INTRON_GCC', 'MAXENTSCR_HSAMODEL_5SS'),
  label = c('Upstream intron length', 'Downstream intron length', 'Exon number',
    'Retained intron length', 'Intron GC content', '5′ splice-site strength'),
  axis = c('Length (nt)', 'Length (nt)', 'Median exon number', 'Length (nt)', 'GC (%)', 'MaxEnt score'))
feature_tests <- inner_join(matt, feature_spec, by = c('event_type', 'feature')) %>% filter(group2 == 'ndiff')
feature_tests$scaled_difference <- feature_tests$difference / feature_tests$background_IQR
feature_tests$label <- factor(paste(feature_tests$event_type, feature_tests$label, sep = ' · '),
  rev(paste(feature_spec$event_type, feature_spec$label, sep = ' · ')))
feature_tests$dataset <- factor(feature_tests$dataset, c('spire', 'fmn2'))
feature_tests$group1 <- factor(feature_tests$group1, c('down', 'up'))
stopifnot(nrow(feature_tests) == 24, all(is.finite(feature_tests$scaled_difference)))
export(feature_tests, 'Reference_feature_effects', 'source_data')
feature_values <- list()
for (ds in names(studies)) for (kind in c('EX', 'IR')) {
  prefix <- if (ds == 'spire') 'spire_fastq' else 'fmn2'
  root <- if (ds == 'spire') 'data/processed/spire_fastq/matt' else 'data/processed/fmn2_matt_final'
  path <- file.path(root, kind, paste0(prefix, '_', kind, '_matt_with_',
    if (kind == 'EX') 'efeatures' else 'ifeatures', '.tab'))
  input_files <- union(input_files, path)
  x <- as.data.frame(fread(path, select = c('EVENT', 'DATASET', feature_spec$feature[feature_spec$event_type == kind]),
    na.strings = c('NA', '', 'nan', '-nan', 'NaN')))
  stopifnot(!anyDuplicated(x$EVENT), all(x$EVENT %in% events[[ds]]$EVENT),
    all(x$DATASET == events[[ds]]$direction[match(x$EVENT, events[[ds]]$EVENT)]))
  x$dataset <- ds; x$event_type <- kind
  feature_values[[paste(ds, kind)]] <- pivot_longer(x, all_of(feature_spec$feature[feature_spec$event_type == kind]),
    names_to = 'feature', values_to = 'value')
}
feature_values <- bind_rows(feature_values)
feature_values <- feature_values[is.finite(feature_values$value), ]
feature_values$DATASET <- factor(feature_values$DATASET, c('ndiff', 'down', 'up'))
feature_values$dataset <- factor(feature_values$dataset, c('spire', 'fmn2'))
# Verify plotting distributions against the exact populations used by the saved tests.
for (i in seq_len(nrow(feature_tests))) {
  z <- feature_tests[i, ]; x <- subset(feature_values, dataset == z$dataset & feature == z$feature & as.character(DATASET) == as.character(z$group1))
  stopifnot(nrow(x) == z$n1, isTRUE(all.equal(median(x$value), z$median1)))
  x <- subset(feature_values, dataset == z$dataset & feature == z$feature & DATASET == 'ndiff')
  stopifnot(nrow(x) == z$n2, isTRUE(all.equal(median(x$value), z$median2)))
}
export(feature_values, 'S7_six_feature_measurements')
# Reconcile whole-intron Matt GC with the event subset used for spatial profiles.
gc_audit <- list()
for (ds in names(studies)) {
  path <- if (ds == 'spire') 'data/processed/spire_fastq/matt/IR/spire_fastq_IR_matt_with_ifeatures.tab' else
    'data/processed/fmn2_matt_final/IR/fmn2_IR_matt_with_ifeatures.tab'
  x <- as.data.frame(fread(path, select = c('EVENT', 'DATASET', 'INTRON_GCC', 'SEQ_INTRON')))
  x$sequence_GC <- sequence_gc(x$SEQ_INTRON); x$SEQ_INTRON <- NULL
  stopifnot(all(abs(x$INTRON_GCC - x$sequence_GC) < 1e-12))
  manifest <- tab(studies[[ds]], 'GC_gradient_IR_event_manifest')
  x$in_spatial_profile <- x$EVENT %in% manifest$EVENT; x$dataset <- ds
  gc_audit[[ds]] <- x
}
gc_audit <- bind_rows(gc_audit)
export(gc_audit, 'S7_intron_GC_sequence_audit')
export(gc_audit %>% group_by(dataset, DATASET) %>% summarise(n = n(),
  Matt_mean_percent = 100*mean(INTRON_GCC), Matt_median_percent = 100*median(INTRON_GCC),
  profile_subset_n = sum(in_spatial_profile),
  profile_subset_mean_percent = 100*mean(sequence_GC[in_spatial_profile]),
  max_absolute_disagreement = max(abs(INTRON_GCC - sequence_GC)), .groups = 'drop'), 'S7_intron_GC_audit_summary')

feature_distribution <- function(ds, feature, title) {
  x <- feature_values[feature_values$dataset == ds & feature_values$feature == feature, ]
  stats <- feature_tests[feature_tests$dataset == ds & feature_tests$feature == feature, ]
  sizes <- table(x$DATASET)
  group_labels <- setNames(paste0(usage_labels[names(sizes)], '\nn=', format(sizes, big.mark = ',', trim = TRUE)), names(sizes))
  q <- paste0(ifelse(stats$group1 == 'down', 'Down', 'Up'), ': q=', formatC(stats$p.BH, format = 'g', digits = 2))
  summary <- x %>% group_by(DATASET) %>% summarise(n = n(), median = median(value),
    q25 = quantile(value, .25), q75 = quantile(value, .75), mean = mean(value),
    minimum = min(value), maximum = max(value), .groups = 'drop')
  stopifnot(sum(summary$n) == nrow(x), all(summary$q25 <= summary$median & summary$median <= summary$q75))
  summary$dataset <- ds; summary$feature <- feature
  export(summary, paste0(toupper(ds), '_', feature, '_display_summary'), 'source_data')
  p <- ggplot(x, aes(DATASET, value, fill = DATASET, color = DATASET)) +
    geom_hline(yintercept = summary$median[summary$DATASET == 'ndiff'], linewidth = .3,
      linetype = 'dotted', color = '#999999') +
    geom_point(data = x[x$DATASET != 'ndiff', ],
      position = position_jitter(width = .14, height = 0, seed = 42), size = .6, alpha = .25) +
    geom_boxplot(width = .48, alpha = .22, outlier.shape = NA, linewidth = .4) +
    geom_point(data = summary, aes(y = median), shape = 21, fill = 'white', size = 1.8, stroke = .5) +
    scale_x_discrete(limits = names(sizes), labels = group_labels) +
    scale_color_manual(values = direction_colors, guide = 'none') +
    scale_fill_manual(values = direction_colors, guide = 'none') +
    labs(x = NULL, y = feature_spec$axis[match(feature, feature_spec$feature)], title = title,
      subtitle = paste(q, collapse = '   '))
  if (grepl('LENGTH|NUMBER', feature)) p <- p + scale_y_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale()))
  if (feature == 'INTRON_GCC') p <- p + scale_y_continuous(labels = scales::label_number(scale = 100))
  if (feature == 'MAXENTSCR_HSAMODEL_5SS') {
    tail <- x %>% group_by(DATASET) %>% summarise(n = n(), below = sum(value < -2), above = sum(value > 12), .groups = 'drop')
    export(tail, paste0(toupper(ds), '_splice_site_display_tails'), 'source_data')
    p <- p + coord_cartesian(ylim = c(-2, 12)) + labs(caption = paste0('Display −2 to 12; ', sum(tail$below + tail$above),
      ' values outside view retained in summaries and tests.')) + theme(plot.caption = element_text(size = 6))
  }
  p
}

rbp_heatmap <- function(ds, kind) {
  all_motifs <- tab('rbp', 'RBP_all_motif_pvalues')
  df <- rbp_panel_data(all_motifs[all_motifs$dataset == ds, ], kind)
  collapsed <- df %>% group_by(RBP, direction, region) %>%
    summarise(pvalue = min(pvalue), score = max(score), .groups = 'drop')
  top <- collapsed %>% group_by(RBP) %>% summarise(score = max(score), .groups = 'drop') %>%
    arrange(desc(score), RBP) %>% slice_head(n = if (kind == 'EX') 30 else 15) %>% pull(RBP)
  regions <- rbp_regions(kind)
  h <- expand_grid(RBP = top, direction = c('up', 'down'), region = regions$region) %>%
    left_join(collapsed, by = c('RBP', 'direction', 'region')) %>% mutate(score = replace_na(score, 0), pvalue = replace_na(pvalue, 1),
      column = paste(region, direction), position = match(direction, c('up', 'down')))
  columns <- as.vector(t(outer(regions$region, c('up', 'down'), paste)))
  wide <- h %>% dplyr::select(RBP, column, score) %>% pivot_wider(names_from = column, values_from = score)
  m <- as.matrix(as.data.frame(wide)[, columns]); rownames(m) <- wide$RBP
  deviation <- apply(m, 1, sd); deviation[deviation == 0] <- 1
  tree <- row_clustering((m - rowMeans(m))/deviation)
  h$row <- tree$y[match(h$RBP, rownames(m))]; h$region <- factor(h$region, regions$region)
  stopifnot(nrow(h) == length(top)*nrow(regions)*2, !anyNA(h$row),
    all(h$pvalue >= 0 & h$pvalue <= 1), all(abs(h$score + log10(pmax(h$pvalue, 1e-300))) < 1e-12))
  export(h, paste0(toupper(ds), '_RBP_', kind, '_heatmap'), 'source_data')
  export(data.frame(RBP = rownames(m)[tree$order]), paste0(toupper(ds), '_RBP_', kind, '_row_order'), 'source_data')
  ann <- distinct(h, position, region, direction)
  labs_region <- if (kind == 'EX') c('Up exon\n3′', 'Up intron\n5′', 'Up intron\n3′',
    'Exon\n5′', 'Exon\n3′', 'Down intron\n5′', 'Down intron\n3′', 'Down exon\n5′') else
    c('Up exon\n5′', 'Up exon\n3′', 'Intron\nbody', 'Down exon\n5′', 'Down exon\n3′')
  annotation <- ggplot(ann, aes(position, 1, fill = direction)) + geom_tile(color = 'white', linewidth = .15) +
    facet_grid(~region, labeller = as_labeller(setNames(labs_region, regions$region))) +
    scale_x_continuous(limits = c(.5, 2.5), expand = c(0, 0)) +
    scale_fill_manual(values = direction_colors, labels = usage_labels, name = NULL) +
    labs(title = paste(toupper(ds), kind, 'rMAPS3 motif P values'), subtitle = 'Minimum raw P; rows clustered') +
    theme_void(base_family = 'Liberation Sans', base_size = 7) +
    theme(strip.text = element_text(size = 6.5), panel.spacing.x = grid::unit(.7, 'mm'),
      legend.position = 'top', legend.key.size = grid::unit(2.5, 'mm'),
      plot.title = element_text(size = 9, face = 'bold'), plot.subtitle = element_text(size = 7))
  heat <- ggplot(h, aes(position, row, fill = score)) + geom_tile(color = 'white', linewidth = .12) +
    facet_grid(~region) + scale_x_continuous(limits = c(.5, 2.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(.5, length(top)+.5), expand = c(0, 0),
      breaks = tree$y, labels = rownames(m)) +
    scale_fill_viridis_c(option = 'cividis', direction = -1, limits = c(0, 8), oob = scales::squish,
      breaks = c(0, 4, 8), labels = c('0', '4', '≥8'), name = '−log10 raw P') + labs(x = NULL, y = NULL) + heat_theme +
    theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 6.5),
      strip.text = element_blank(), strip.background = element_blank(), legend.position = 'right',
      panel.spacing.x = grid::unit(.7, 'mm'), plot.margin = margin(0, 7, 0, 6))
  wrap_elements(full = annotation / heat +
    plot_layout(heights = grid::unit(c(3, 1), c('mm', 'null')))) +
    theme(plot.margin = margin(10, 0, 0, 0))
}

# Coarse IR bins and separate regions keep small event sets and mixed sequence scales readable.
gc_plot <- function(ds, kind, difference = FALSE) {
  level <- if (kind == 'IR') 'coarse' else 'fine'
  x <- tab(studies[[ds]], paste0('GC_gradient_', kind, '_', level, '_directions'))
  stopifnot(!anyDuplicated(x[, c('category', 'bin')]), all(is.finite(x$mean_GC)),
    all(is.finite(x$se) & x$se >= 0), all(x$n >= 2))
  regions <- x %>% group_by(region) %>% summarise(start = min(bin)-.5, end = max(bin)+.5, .groups = 'drop') %>% arrange(start)
  regions$label <- if (kind == 'EX') c('Upstream exon\nlast ≤200 nt', 'Upstream intron\nscaled', 'Target exon\nscaled',
    'Downstream intron\nscaled', 'Downstream exon\nfirst ≤200 nt') else
    c('Upstream exon\nlast ≤150 nt', 'Intron 5′ end\nfirst ≤80 nt', 'Intron body\nscaled',
      'Intron 3′ end\nlast ≤80 nt', 'Downstream exon\nfirst ≤150 nt')
  i <- match(x$region, regions$region)
  x$position <- (x$bin-regions$start[i])/(regions$end[i]-regions$start[i])
  x$layout <- level
  x$region <- factor(x$region, regions$region)
  stopifnot(all(x$position > 0 & x$position < 1), nrow(regions) == 5)
  export(x, paste0(toupper(ds), '_GC_', kind, '_profiles'), 'source_data')
  if (kind == 'IR') export(tab(studies[[ds]], 'GC_gradient_IR_fine_directions'),
    paste0(toupper(ds), '_GC_IR_fine_reference'), 'source_data')
  if (difference) {
    bg <- x[x$category == 'ndiff', c('bin', 'mean_GC', 'se')]; names(bg)[-1] <- c('background', 'background_se')
    x <- merge(x[x$category != 'ndiff', ], bg, by = 'bin')
    x$mean_GC <- x$mean_GC - x$background; x$se <- sqrt(x$se^2 + x$background_se^2)
    export(x, paste0(toupper(ds), '_GC_', kind, '_difference'), 'source_data')
  }
  sizes <- x %>% group_by(category) %>% summarise(low = min(n), high = max(n), .groups = 'drop')
  breaks <- intersect(c('ndiff', 'down', 'up'), sizes$category)
  labels <- setNames(paste0(usage_labels[sizes$category], ' (n=',
    format(sizes$low, big.mark = ',', trim = TRUE), '–', format(sizes$high, big.mark = ',', trim = TRUE), ')'), sizes$category)
  colors <- direction_colors; colors['ndiff'] <- '#555555'
  p <- ggplot(x, aes(position, mean_GC, color = category, fill = category)) +
    facet_grid(~region, labeller = as_labeller(setNames(regions$label, regions$region))) +
    geom_hline(yintercept = if (difference) 0 else numeric(), linewidth = .3, color = '#777777') +
    geom_ribbon(aes(ymin = mean_GC - 1.96*se, ymax = mean_GC + 1.96*se), alpha = .12, color = NA) +
    geom_line(aes(linetype = category), linewidth = .6) +
    scale_color_manual(values = colors, breaks = breaks, labels = labels[breaks], name = NULL) +
    scale_fill_manual(values = colors, breaks = breaks, labels = labels[breaks], name = NULL) +
    scale_linetype_manual(values = c(ndiff = 'dashed', down = 'solid', up = 'solid'),
      breaks = breaks, labels = labels[breaks], name = NULL) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, .5, 1), labels = c('0', '50', '100'), expand = c(0, 0)) +
    labs(x = 'Position within region (%)', y = if (difference) 'ΔGC versus background (pp)' else 'GC (%)',
      title = paste(toupper(ds), '·', kind, 'GC profile'),
      subtitle = 'Mean ±1.96 SE across events; n gives the range across bins') +
    theme(axis.text.x = element_text(size = 6), strip.text = element_text(size = 6.5),
      panel.spacing.x = grid::unit(4, 'mm'), legend.key.width = grid::unit(5, 'mm'))
  if (kind == 'IR' && !difference) {
    stopifnot(all(x$mean_GC - 1.96*x$se >= 25), all(x$mean_GC + 1.96*x$se <= 65))
    p <- p + scale_y_continuous(limits = c(25, 65), breaks = c(30, 40, 50, 60))
  }
  p
}
go_dotplot <- function() {
  x <- go_revised %>% filter(dataset == 'spire', background == 'tested_N10', selection == 'selected') %>%
    group_by(ontology) %>% arrange(q_BH, pvalue, ID, .by_group = TRUE) %>% slice_head(n = 3) %>% ungroup()
  x <- x[order(match(x$ontology, c('BP', 'CC', 'MF')), x$q_BH, x$pvalue, x$ID), ]
  x$GeneRatio <- x$Count/x$query_genes
  x$term <- factor(x$ID, rev(x$ID))
  x$ontology_label <- factor(x$ontology, c('BP', 'CC', 'MF'), c('Biological\nprocess', 'Cellular\ncomponent', 'Molecular\nfunction'))
  stopifnot(nrow(x) == 9, !anyDuplicated(x$ID), all(x$q_BH > .05))
  export(x, 'SPIRE_main_GO', 'source_data')
  ggplot(x, aes(GeneRatio, term, size = Count, color = q_BH)) + geom_point(alpha = .95) +
    facet_grid(ontology_label~., scales = 'free_y', space = 'free_y') +
    scale_y_discrete(labels = setNames(vapply(x$Description, function(label)
      paste(strwrap(label, 48), collapse = '\n'), character(1)), x$ID)) +
    scale_x_continuous(labels = scales::label_percent(accuracy = 1), expand = expansion(mult = c(.12, .15))) +
    scale_size_area(max_size = 4, breaks = c(10, 30, 50), name = 'Genes') +
    scale_color_gradient(low = '#6A3D9A', high = '#CFD3DB', trans = 'log10', limits = c(.05, 1),
      breaks = c(.05, .1, .5, 1), labels = c('0.05', '0.1', '0.5', '1'), name = 'BH q') +
    guides(color = guide_colorbar(order = 1, barheight = grid::unit(20, 'mm')),
      size = guide_legend(order = 2)) +
    labs(x = 'Gene ratio', y = NULL, title = 'Top-ranked GO terms in selected splicing genes',
      subtitle = 'Tested-gene background; top 3 terms per ontology; none pass BH q ≤0.05') +
    theme(legend.position = 'right', strip.text.y = element_text(angle = 0, size = 6.5),
      axis.text.y = element_text(size = 6.5), panel.spacing.y = grid::unit(1.5, 'mm'))
}

stability_plot <- function(ds) {
  x <- sensitivity[[ds]] %>% group_by(event_type) %>% summarise(
    `Same direction` = mean(loo_same_direction), `At least 10 pp` = mean(loo_min_signed_pp >= 10),
    .groups = 'drop') %>% pivot_longer(-event_type)
  x$event_type <- factor(x$event_type, c('EX', 'IR', 'Alt5', 'Alt3'))
  export(x, paste0(toupper(ds), '_display_stability'), 'source_data')
  ggplot(x, aes(event_type, value, fill = name)) + geom_col(position = 'dodge', width = .65) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .25), labels = scales::label_percent()) +
    scale_fill_manual(values = c('#999999', '#222222'), name = NULL) +
    labs(x = 'Event class', y = 'Selected events', title = paste(toupper(ds), 'leave-one-sample stability'))
}

p_spire_psi <- psi_heatmap('spire'); p_fmn2_psi <- psi_heatmap('fmn2')
p_spire_rbp <- rbp_heatmap('spire', 'EX'); p_fmn2_rbp <- rbp_heatmap('fmn2', 'EX')
save_figure(p_ma + event_counts('spire') + p_spire_psi + free(p_division, side = 'l') +
  plot_layout(design = 'AB\nCC\nDD', heights = c(1, 1.35, .85)), 'SPIRE_Figure_1',
  'SPIRE: expression and splicing changes, including events in division-annotated genes', 225)
save_figure(free(go_dotplot(), side = 'l') +
  feature_distribution('spire', 'UPINTRON_MEDIANLENGTH', 'EX: longer upstream introns') +
  feature_distribution('spire', 'INTRON_GCC', 'IR: higher whole-intron GC') + p_spire_rbp +
  plot_layout(design = 'AA\nBC\nDD', heights = c(1.1, 1, 1.65)), 'SPIRE_Figure_2',
  'SPIRE: functional context, sequence architecture and regional RBP motifs', 270)
save_figure(event_counts('fmn2') + p_fmn2_psi +
  feature_distribution('fmn2', 'UPINTRON_MEDIANLENGTH', 'EX: upstream intron length') +
  feature_distribution('fmn2', 'DOINTRON_MEDIANLENGTH', 'EX: downstream intron length') + p_fmn2_rbp +
  plot_layout(design = 'AB\nCD\nEE', heights = c(1, .9, 1.6)), 'FMN2_Figure_1',
  'FMN2: splicing profiles, exon-flanking architecture and regional RBP motifs', 245)
save_figure(polish(p2a, 'Overlap on the common tested universe') + polish(p2b, 'Global effect comparison') +
  free(polish(p2c, 'Effects at the 11 shared events'), side = 'l') +
  plot_layout(design = 'AB\nCC', heights = c(1, 1.05)), 'Convergence_Figure_1',
  'Limited shared targeting without a global concordant splicing response', 170)

# Three shared supplements: essential QC, complementary features, and intron sequence context.
expr_m <- as.matrix(expr[, names(groups$spire)]); rownames(expr_m) <- expr$gene_id
expr_tree <- row_clustering(expr_m - rowMeans(expr_m))
expr_long$gene_id <- factor(expr_long$gene_id, rev(rownames(expr_m)[expr_tree$order]))
export(expr_long, 'SPIRE_clustered_DEG_VST', 'source_data')
p_expr <- ggplot(expr_long, aes(sample, gene_id, fill = centered)) + geom_raster() +
  scale_fill_gradient2(low = signed_colors[1], mid = signed_colors[2], high = signed_colors[3],
    limits = c(-3, 3), oob = scales::squish, name = 'Centered VST') +
  scale_x_discrete(labels = function(s) sub('SPIRE_', '', s)) +
  labs(x = NULL, y = NULL, title = 'Expression of 111 selected genes') + heat_theme +
  theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
save_figure(polish(pca_plot('Expression_PCA', ''), 'SPIRE expression PCA') +
  polish(pca_plot('Splicing_PCA', ''), 'SPIRE splicing PCA') +
  (polish(ps1d, 'VAST support at target exons') + heat_theme +
    scale_x_discrete(labels = function(s) sub('SPIRE_(Cont|DKO)_', '', s)) +
    labs(x = 'Controls 1/2/7; DKO 4/5/6') + theme(axis.text.x = element_text(angle = 0))) +
  p_expr + stability_plot('spire') + stability_plot('fmn2') + plot_layout(ncol = 2),
  'Figure_S1', 'Sample quality, target-exon support, expression and splicing stability', 240)
save_figure(feature_distribution('spire', 'MEDIAN_EXON_NUMBER', 'SPIRE · EX exon number') +
  feature_distribution('fmn2', 'MEDIAN_EXON_NUMBER', 'FMN2 · EX exon number') +
  feature_distribution('spire', 'MAXENTSCR_HSAMODEL_5SS', 'SPIRE · IR 5′ splice-site strength') +
  feature_distribution('fmn2', 'MAXENTSCR_HSAMODEL_5SS', 'FMN2 · IR 5′ splice-site strength') +
  plot_layout(ncol = 2), 'Figure_S2', 'Additional exon architecture and intron splice-site context', 180)
save_figure(rbp_heatmap('spire', 'IR') + rbp_heatmap('fmn2', 'IR') +
  gc_plot('spire', 'IR') + gc_plot('fmn2', 'IR') +
  plot_layout(ncol = 1, heights = c(1, 1, .8, .8)),
  'Figure_S3', 'Intron rMAPS3 motif P values and GC profiles in SPIRE and FMN2', 270)
# Reference profiles remain in source tables without additional figures.
for (ds in names(studies)) for (kind in c('EX', 'IR')) invisible(gc_plot(ds, kind, TRUE))
