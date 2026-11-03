# Compositions reuse validated analysis outputs; no statistical selection is changed here.
polish <- function(p, title) p + labs(title = title, subtitle = NULL, caption = NULL)
heat_theme <- theme(panel.grid = element_blank(), panel.grid.major.y = element_blank(),
  panel.grid.major.x = element_blank(), axis.ticks = element_blank(),
  axis.line = element_blank(), panel.border = element_blank())
usage_labels <- c(down = 'Decreased', ndiff = 'Background', up = 'Increased')
study_labels <- c(spire = 'SPIRE', fmn2 = 'FMN2')

# Abundance exposes expression effect sizes without the single Mt1 P-value dominating.
ma <- dge[is.finite(dge$baseMean) & dge$baseMean > 0 & is.finite(dge$log2FoldChange), ]
ma <- ma[order(ma$selected), ]
export(ma, 'Figure1B_expression_MA', 'source_data')
p_ma <- ggplot(ma, aes(baseMean, log2FoldChange)) +
  geom_hline(yintercept = 0, color = '#AAAAAA', linewidth = .3) +
  geom_point(aes(color = direction), size = .65, alpha = .65) +
  ggrepel::geom_text_repel(data = ma[ma$gene_name %in% c('Spire1', 'Spire2', 'Mt1'), ],
    aes(label = gene_name), seed = 42, size = 2.5, fontface = 'italic', family = 'Liberation Sans',
    min.segment.length = 0, max.overlaps = Inf) +
  scale_x_log10(breaks = c(1, 100, 10000), labels = scales::label_number()) +
  direction_scale() + guides(color = guide_legend(override.aes = list(size = 2, alpha = 1), nrow = 1)) +
  labs(x = 'Mean normalized counts', y = expression(log[2]~'fold change'),
    title = '111 differentially expressed genes')

division_top <- head(eligible_examples[!duplicated(eligible_examples$GENE), ], 15)
division_heat <- sample_psi(division_top, 'spire')
division_heat$GENE <- factor(division_heat$GENE, rev(division_top$GENE))
division_heat$sample <- factor(division_heat$sample, names(groups$spire))
division_heat$group <- factor(division_heat$group, c('Cont', 'DKO'), c('Control', 'Spire1/2 DKO'))
export(division_top, 'Figure1D_candidate_selection', 'source_data')
export(division_heat, 'Figure1D_division_PSI', 'source_data')
stopifnot(nrow(division_heat) == 90, !anyDuplicated(division_top$GENE),
  all(division_top$loo_min_signed_pp > 0))
p_division <- ggplot(division_heat, aes(sample, GENE, fill = PSI)) +
  geom_tile(color = 'white', linewidth = .25) +
  facet_grid(~group, scales = 'free_x', space = 'free_x') +
  scale_fill_gradientn(colors = c('#F7FBFF', '#C6DBEF', '#6BAED6', '#2171B5', '#08306B'),
    limits = c(0, 100), breaks = c(0, 50, 100), name = 'PSI (%)') +
  scale_x_discrete(labels = function(s) sub('SPIRE_(Cont|DKO)_', '', s)) +
  labs(x = 'Animal', y = NULL, title = 'Division-associated candidates') + heat_theme +
  theme(axis.text.y = element_text(face = 'italic', size = 7),
    legend.position = 'right', panel.spacing.x = grid::unit(2, 'mm'))

example_psi$label <- factor(example_psi$GENE, examples$GENE,
  paste0(examples$GENE, '  ·  ', examples$event_type, '  ·  ΔPSI ', sprintf('%+.1f pp', 100*examples$deltapsi)))
p_examples <- psi_points(example_psi, 'Individual-animal profiles of prioritized candidates', NULL) +
  scale_x_discrete(labels = c(Cont = 'Control', DKO = 'Spire1/2 DKO')) +
  theme(strip.text = element_text(face = 'italic', size = 8))
p_target <- polish(p1a, 'RNA signal at the targeted loci') +
  scale_x_discrete(labels = c(Cont = 'Control', DKO = 'DKO')) +
  theme(strip.text = element_text(face = 'italic', size = 8))
p_counts <- polish(p1c, '1,547 selected splicing events') +
  scale_x_discrete(labels = c(EX = 'Exon', IR = 'Intron', Alt5 = "Alt 5′", Alt3 = "Alt 3′"))
save_figure(p_target + p_ma + p_counts + p_division + p_examples +
  plot_layout(design = 'AB\nCD\nEE', heights = c(1, 1.35, .95)),
  'Figure_1', 'SPIRE loss alters RNA expression and splicing, including division-associated genes', 210)

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
export(feature_tests, 'Figure2A_feature_effects', 'source_data')
p_features <- ggplot(feature_tests, aes(group1, label, fill = scaled_difference)) +
  geom_tile(color = 'white', linewidth = .6) +
  geom_text(aes(label = ifelse(p.BH <= .05, '•', '')), size = 4, color = '#111111') +
  facet_grid(~dataset, labeller = as_labeller(study_labels)) +
  scale_x_discrete(labels = usage_labels) +
  scale_fill_gradient2(low = signed_colors[1], mid = signed_colors[2], high = signed_colors[3],
    limits = c(-1, 1), oob = scales::squish, breaks = c(-1, 0, 1), name = 'Median shift\n/ background IQR') +
  labs(x = NULL, y = NULL, title = 'Sequence architecture of affected events') + heat_theme +
  theme(legend.position = 'right', strip.text = element_text(size = 8))

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
feature_violin <- function(feature, title) {
  x <- feature_values[feature_values$feature == feature, ]
  stats <- feature_tests[feature_tests$feature == feature, ]
  stats$DATASET <- stats$group1
  stats$annotation <- paste0('q=\n', formatC(stats$p.BH, format = 'g', digits = 2))
  n <- x %>% count(dataset, DATASET)
  is_log <- grepl('LENGTH|NUMBER', feature)
  n$y <- if (is_log) min(x$value)/2 else min(x$value) - .09*diff(range(x$value))
  p <- ggplot(x, aes(DATASET, value, fill = DATASET)) +
    geom_violin(scale = 'width', alpha = .65, linewidth = .3, color = '#777777', trim = TRUE) +
    geom_boxplot(width = .16, fill = 'white', outlier.shape = NA, linewidth = .3) +
    geom_text(data = n, aes(y = y, label = paste0('n=\n', n)), vjust = .5, size = 2.1) +
    geom_text(data = stats, aes(y = Inf, label = annotation), vjust = 1.3, size = 2.15) +
    facet_wrap(~dataset, nrow = 1, labeller = as_labeller(study_labels)) +
    scale_fill_manual(values = direction_colors, guide = 'none') +
    scale_x_discrete(labels = c(ndiff = 'Bg', down = 'Down', up = 'Up')) +
    labs(x = NULL, y = feature_spec$axis[match(feature, feature_spec$feature)], title = title)
  if (is_log) {
    stopifnot(all(x$value > 0))
    p <- p + scale_y_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale()), expand = expansion(mult = c(.1, .2)))
  } else p <- p + scale_y_continuous(expand = expansion(mult = c(.1, .2)),
    labels = scales::label_number(scale = if (feature == 'INTRON_GCC') 100 else 1))
  p
}

rbp_heatmap <- function(kind) {
  h <- tab('rbp', paste0('RBP_', kind, '_paper_heatmap_scores'))
  regions <- rbp_regions(kind)
  columns <- unlist(lapply(regions$region, function(r)
    paste(rep(c('spire', 'fmn2'), each = 2), rep(c('up', 'down'), 2), r, sep = ':')))
  wide <- h %>% dplyr::select(RBP, column, score) %>% pivot_wider(names_from = column, values_from = score)
  m <- as.matrix(as.data.frame(wide)[, columns]); rownames(m) <- wide$RBP
  stopifnot(!anyNA(m), all(is.finite(m)), nrow(m) == if (kind == 'EX') 30 else 15,
    all(abs(h$score + log10(pmax(h$pvalue, 1e-300))) < 1e-10))
  rbp_order <- rownames(m)[hclust(dist(m), method = 'complete')$order]
  h$RBP <- factor(h$RBP, rev(rbp_order)); h$column <- factor(h$column, columns)
  h$region <- factor(h$region, regions$region)
  h$position <- match(paste(h$dataset, h$direction), c('spire up', 'spire down', 'fmn2 up', 'fmn2 down'))
  stopifnot(!anyNA(h$position))
  export(h, paste0('RBP_', kind, '_heatmap'), 'source_data')
  ann <- h %>% distinct(column, position, region, dataset, direction)
  ann <- bind_rows(transmute(ann, column, position, region, track = 'Study', value = dataset),
    transmute(ann, column, position, region, track = 'Usage', value = direction))
  ann$track <- factor(ann$track, c('Usage', 'Study'))
  labs_region <- if (kind == 'EX') c('Up exon\n3′', 'Up intron\n5′', 'Up intron\n3′',
    'Exon\n5′', 'Exon\n3′', 'Down intron\n5′', 'Down intron\n3′', 'Down exon\n5′') else
    c('Up exon\n5′', 'Up exon\n3′', 'Intron\nbody', 'Down exon\n5′', 'Down exon\n3′')
  annotation <- ggplot(ann, aes(position, track, fill = value)) + geom_tile(color = 'white', linewidth = .15) +
    facet_grid(~region, labeller = as_labeller(setNames(labs_region, regions$region))) +
    scale_x_continuous(limits = c(.5, 4.5), expand = c(0, 0)) +
    scale_fill_manual(values = c(dataset_colors, direction_colors[c('up', 'down')]),
      breaks = c('spire', 'fmn2', 'up', 'down'), labels = c('SPIRE', 'FMN2', 'Increased', 'Decreased'), name = NULL) +
    labs(title = paste(if (kind == 'EX') 'Exon' else 'Intron', 'RBP motif landscape'), x = NULL, y = NULL) +
    theme_void(base_family = 'Liberation Sans', base_size = 7) +
    theme(axis.text.y = element_text(size = 6), strip.text = element_text(size = 6.5),
      panel.spacing.x = grid::unit(.7, 'mm'), legend.position = 'top',
      legend.key.size = grid::unit(2.5, 'mm'), plot.title = element_text(size = 9, face = 'bold'))
  heat <- ggplot(h, aes(position, RBP, fill = score)) + geom_tile(color = 'white', linewidth = .12) +
    facet_grid(~region) + scale_x_continuous(limits = c(.5, 4.5), expand = c(0, 0)) +
    scale_fill_viridis_c(option = 'cividis', direction = -1, limits = c(0, 8), oob = scales::squish,
      breaks = c(0, 2, 4, 6, 8), name = '−log10 raw P') +
    labs(x = NULL, y = NULL) + heat_theme +
    theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 6.5),
      strip.text = element_blank(), strip.background = element_blank(), legend.position = 'right',
      panel.spacing.x = grid::unit(.7, 'mm'), plot.margin = margin(0, 7, 5, 7))
  # wrap_elements keeps the annotation and matrix as one lettered manuscript panel.
  wrap_elements(full = wrap_plots(annotation, heat, ncol = 1, heights = c(.42, if (kind == 'EX') 6 else 3.8)))
}
p_ex_rbp <- rbp_heatmap('EX')
p_ir_rbp <- rbp_heatmap('IR')
save_figure(free(p_features, side = 'l') + feature_violin('UPINTRON_MEDIANLENGTH', 'Longer flanking introns') +
  feature_violin('MAXENTSCR_HSAMODEL_5SS', 'Weaker 5′ splice sites') + p_ex_rbp +
  plot_layout(design = 'AA\nBC\nDD', heights = c(.85, 1, 1.85)),
  'Figure_2', 'Affected events have distinctive sequence features and regional RBP motif patterns', 250)

save_figure(polish(pca_plot('Expression_PCA', ''), 'Expression PCA') +
  polish(pca_plot('Splicing_PCA', ''), 'Splicing PCA') + polish(ps1c, 'Read alignment') +
  (polish(ps1d, 'Support at targeted exons') + heat_theme +
    scale_x_discrete(labels = function(s) sub('SPIRE_(Cont|DKO)_', '', s)) +
    labs(x = 'Animal: controls 1/2/7; DKO 4/5/6') + theme(axis.text.x = element_text(angle = 0, hjust = .5))) + plot_layout(ncol = 2, heights = c(1, 1.25)),
  'Figure_S1', 'Sample quality and targeted-exon support', 185)
save_figure(polish(ps2a, 'Fitted and observed changes') + polish(ps2b, 'Leave-one-animal sensitivity') +
  polish(ps2c, 'Selection thresholds') + (polish(p1d, 'All selected splicing events') + heat_theme) +
  plot_layout(ncol = 2), 'Figure_S2', 'Splicing profiles and sensitivity to animals and thresholds', 185)
save_figure(polish(p2a, 'Shared and study-specific targets') + polish(p2b, 'Global raw-PSI comparison') +
  free(polish(p2c, 'Effects in the 11 shared events'), side = 'l') +
  free(polish(p2d, 'Individual samples at shared events') + heat_theme, side = 'l') +
  plot_layout(design = 'AB\nCC\nDD', heights = c(1, 1.1, 1.1)),
  'Figure_S3', 'SPIRE and FMN2 share targets with heterogeneous effects', 235)
save_figure(wrap_plots(lapply(seq_len(nrow(feature_spec)), function(i)
  feature_violin(feature_spec$feature[i], paste(feature_spec$event_type[i], feature_spec$label[i], sep = ' · '))),
  ncol = 2), 'Figure_S4', 'Full distributions for the six reference sequence features', 240)

gc_plot <- function(kind) {
  x <- bind_rows(lapply(names(studies), function(ds) {
    z <- tab(studies[[ds]], paste0('GC_gradient_', kind, '_fine_directions')); z$dataset <- ds; z
  }))
  x$dataset <- factor(x$dataset, c('spire', 'fmn2'))
  regions <- x %>% group_by(region) %>% summarise(start = min(bin) - .5, end = max(bin) + .5, .groups = 'drop') %>% arrange(start)
  regions$label <- if (kind == 'EX') c('Up exon', 'Up intron', 'Exon', 'Down intron', 'Down exon') else
    c('Up exon', '5′ flank', 'Intron body', '3′ flank', 'Down exon')
  regions$exon <- grepl('exon', regions$label, ignore.case = TRUE)
  stopifnot(all(x$mean_GC >= 0 & x$mean_GC <= 100), all(x$se >= 0),
    nrow(regions) == 5, max(x$bin) == if (kind == 'EX') 68 else 88)
  export(x, paste0('GC_', kind, '_profiles'), 'source_data')
  ggplot(x, aes(bin, mean_GC, color = category, fill = category)) +
    geom_rect(data = regions[regions$exon, ], aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE, fill = '#EEEEEE') +
    geom_vline(xintercept = head(regions$end, -1), color = '#AAAAAA', linetype = 'dotted', linewidth = .3) +
    geom_ribbon(aes(ymin = mean_GC - 1.96*se, ymax = mean_GC + 1.96*se), alpha = .14, color = NA) +
    geom_line(linewidth = .55) + facet_wrap(~dataset, ncol = 1, labeller = as_labeller(study_labels)) +
    scale_color_manual(values = direction_colors, labels = usage_labels, name = NULL) +
    scale_fill_manual(values = direction_colors, labels = usage_labels, name = NULL) +
    scale_x_continuous(breaks = (regions$start + regions$end)/2, labels = regions$label, expand = c(0, 0)) +
    labs(x = 'Scaled sequence regions', y = 'GC content (%)',
      title = paste(if (kind == 'EX') 'Exon' else 'Intron', 'sequence composition'))
}
save_figure(gc_plot('EX') / gc_plot('IR'), 'Figure_S5', 'GC composition across exon and intron neighborhoods', 235)

motif_distribution <- function(kind) {
  x <- tab('rbp', paste0('RBP_', kind, '_paper_motif_scores'))
  x$dataset <- factor(x$dataset, c('spire', 'fmn2'))
  x$direction <- factor(x$direction, c('down', 'up'))
  x$region <- factor(x$region, rbp_regions(kind)$region)
  export(x, paste0('RBP_', kind, '_regional_scores'), 'source_data')
  ggplot(x, aes(dataset, score, fill = direction)) +
    geom_point(aes(color = direction), position = position_jitterdodge(jitter.width = .18,
      dodge.width = .7, seed = 42), size = .4, alpha = .15) +
    geom_boxplot(width = .6, outlier.shape = NA, alpha = .7, linewidth = .25, position = position_dodge(.7)) +
    geom_hline(yintercept = -log10(.05), linetype = 'dotted', color = '#555555', linewidth = .3) +
    facet_wrap(~region, nrow = 1) +
    scale_fill_manual(values = direction_colors, labels = usage_labels, name = NULL) +
    scale_color_manual(values = direction_colors, guide = 'none') +
    scale_x_discrete(labels = study_labels) +
    labs(x = 'Regions correspond to the annotated heatmaps', y = '−log10 raw P',
      title = paste(kind, 'motif-score distributions')) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = 6))
}
save_figure(p_ir_rbp / motif_distribution('EX') / motif_distribution('IR') +
  plot_layout(heights = c(1.5, 1, 1)), 'Figure_S6', 'Intron RBP landscape and regional motif-score distributions', 250)
save_figure(polish(p1b, 'Expression significance and effect') +
  (polish(ps3a, 'Individual-animal expression') + heat_theme) +
  polish(ps3b, 'Expression and splicing targets') + free(polish(ps4b, 'Predicted transcript impact'), side = 'l') +
  polish(ps3c, 'GO background sensitivity') + free(polish(ps3d, 'Division-category tests'), side = 'l') +
  plot_layout(ncol = 2), 'Figure_S7', 'Expression profiles, transcript impact and functional interpretation', 240)
