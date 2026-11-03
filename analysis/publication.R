#!/usr/bin/env Rscript
# Assemble the manuscript from validated current results; never overwrite an analysis.
.libPaths(c(normalizePath('analysis/renv/library'), .libPaths()))
options(warn = 1)
source('analysis/common.R')
source('analysis/comparison_helpers.R')
source('analysis/publication_helpers.R')
source('analysis/rbp_helpers.R')
# Publication-only palette, adapted to the requested Lena/RBP visual references.
direction_colors <- c(down = '#2166AC', ndiff = '#BDBDBD', up = '#B2182B')
dataset_colors <- c(spire = '#6A3D9A', fmn2 = '#33A02C')
condition_colors <- c(Cont = '#777777', DKO = '#6A3D9A')
signed_colors <- c('#2166AC', '#F7F7F7', '#B2182B')
suppressPackageStartupMessages(library(patchwork))
args <- commandArgs(TRUE)
stopifnot(length(args) <= 1)
out <- if (length(args)) args[1] else 'results/publication'
if (file.exists(out)) stop('Choose a new output directory: ', out)
for (d in c('figures', 'tables', 'source_data', 'text')) dir.create(file.path(out, d), recursive = TRUE)
input_files <- character()
read_input <- function(path) {
  stopifnot(file.exists(path))
  input_files <<- union(input_files, path)
  as.data.frame(fread(path))
}
tab <- function(study, name) read_input(file.path('reports', study, 'tables', paste0(name, '.csv')))
export <- function(x, name, directory = 'tables') fwrite(x, file.path(out, directory, paste0(name, '.csv')), na = 'NA')
studies <- c(spire = 'spire_final', fmn2 = 'fmn2_final')
events <- groups <- list()
for (ds in names(studies)) {
  root <- file.path('reports', studies[[ds]])
  stopifnot(file.exists(file.path(root, 'COMPLETE')))
  path <- file.path(root, 'tables/analysis_parameters.json')
  input_files <- c(input_files, path, file.path(root, 'COMPLETE'))
  parameters <- jsonlite::fromJSON(path)
  stopifnot(parameters$dataset == if (ds == 'spire') 'spire_fastq' else 'fmn2',
    !parameters$legacy, parameters$fpr_max == .05, parameters$abs_deltapsi_min == .1,
    parameters$min_reads_all_samples == 10)
  groups[[ds]] <- unlist(parameters$samples)
  events[[ds]] <- tab(studies[[ds]], 'all_events')
  x <- events[[ds]]
  stopifnot(!anyDuplicated(x$EVENT), all(classify(x$FDR, x$deltapsi) == x$direction),
    all(abs(x$betasPsiB - x$betasPsiA - x$deltapsi) < 1e-8))
}
sheet <- read_input('config/spire_fastq.tsv')
stopifnot(identical(sheet$sample, names(groups$spire)), identical(sheet$condition, unname(groups$spire)))
sheet$experimental_unit <- 'independent biological replicate'
sheet$metadata_source <- 'Investigator confirmation, 2026-11-02'
export(sheet, 'S1_SPIRE_samples')
spire <- events$spire; fmn2 <- events$fmn2
original <- read_input('data/inputs/spire_fastq/betas_N10_nsim1000/csv/spire_fastq_all_events.csv')
stopifnot(identical(original$EVENT, spire$EVENT))
for (k in names(original)) stopifnot(isTRUE(all.equal(original[[k]], spire[[k]], tolerance = 1e-12)))
rm(original)
common <- join_common_events(spire, fmn2)
previous <- tab('comparison', 'common_tested_events')
stopifnot(setequal(common$EVENT, previous$EVENT))
for (k in c('spire_deltapsi', 'fmn2_deltapsi', 'spire_selected', 'fmn2_selected'))
  stopifnot(isTRUE(all.equal(common[[k]], previous[[k]][match(common$EVENT, previous$EVENT)])))
common <- previous
gene_degrees <- common[common$GENE %in% valid_genes(common$GENE), ] %>% group_by(GENE) %>%
  summarise(tested_events = n(), spire_selected = any(spire_selected), fmn2_selected = any(fmn2_selected), .groups = 'drop')
degree_overlap <- degree_matched_overlap(gene_degrees$spire_selected, gene_degrees$fmn2_selected, gene_degrees$tested_events)
export(gene_degrees, 'S5_gene_event_counts')
export(data.frame(observed = degree_overlap$observed, mean_null = degree_overlap$expected,
  pvalue = degree_overlap$pvalue, permutations = length(degree_overlap$null), seed = 42,
  strata = 'floor(log2(number of common tested events per gene))'), 'S5_event_count_matched_gene_overlap')
export(data.frame(shared_genes = degree_overlap$null), 'S5_event_count_matched_null')
dge <- tab('spire_review', 'DGE_all_genes')
original_dge <- read_input('data/inputs/spire_fastq/dge/tables/differential/rnaseq/DKO_vs_Cont_SPIRE_DKO.deseq2.results.tsv')
stopifnot(identical(original_dge$gene_id, dge$gene_id))
for (k in c('log2FoldChange', 'pvalue', 'padj')) stopifnot(isTRUE(all.equal(original_dge[[k]], dge[[k]])))
dge$selected <- !is.na(dge$padj) & dge$padj <= .05
dge$direction <- ifelse(dge$selected, ifelse(dge$log2FoldChange > 0, 'up', 'down'), 'ndiff')
dds_path <- 'data/inputs/spire_fastq/dge/other/deseq2/rnaseq/DKO_vs_Cont_SPIRE_DKO.dds.rld.rds'
dds <- readRDS(dds_path); input_files <- c(input_files, dds_path)
design <- paste(deparse(attr(dds, 'design')), collapse = '')
stopifnot(design == '~0 + condition', attr(dds, 'test') == 'Wald')
export(data.frame(design = design, contrast = 'DKO - Cont', test = attr(dds, 'test'),
  fit = attr(attr(dds, 'dispersionFunction'), 'fitType'), source = dds_path), 'S1_DESeq2_design')
rm(dds)

sensitivity <- lapply(names(events), function(ds) {
  x <- events[[ds]]; x <- x[x$direction != 'ndiff', ]
  g <- groups[[ds]]; control <- if (ds == 'spire') 'Cont' else 'Het'
  cbind(x, psi_sensitivity(as.matrix(x[, names(g)[g == control]]),
    as.matrix(x[, names(g)[g != control]]), x$deltapsi))
})
names(sensitivity) <- names(events)
robust <- sensitivity$spire
export(robust, 'S3_SPIRE_selected_with_sensitivity')
export(sensitivity$fmn2, 'S5_FMN2_selected_with_sensitivity')
export(spire, 'S3_SPIRE_all_tested_events')
export(fmn2, 'S5_FMN2_all_tested_events')
export(dge, 'S2_SPIRE_all_gene_expression')
export(dge[dge$selected, ], 'S2_SPIRE_selected_gene_expression')
export(common, 'S5_common_tested_events')

# Define division relevance from the question in the email; retain every matching event.
suppressPackageStartupMessages(library(org.Mm.eg.db))
division_terms <- c('GO:0051301' = 'cell division', 'GO:0051321' = 'meiotic cell cycle',
  'GO:0007051' = 'spindle organization', 'GO:0007059' = 'chromosome segregation')
division <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db, keys = names(division_terms),
  keytype = 'GOALL', columns = c('SYMBOL', 'EVIDENCEALL')))
division <- unique(division[!is.na(division$SYMBOL), ])
division$term <- unname(division_terms[division$GOALL])
export(division, 'S4_division_GO_membership')
division_genes <- unique(division$SYMBOL)
division_selected <- robust[robust$GENE %in% division_genes, ]
division_selected <- division_selected[order(-division_selected$loo_min_signed_pp,
  -abs(division_selected$raw_delta_pp), division_selected$EVENT), ]
export(division_selected, 'S4_division_selected_splicing')
export(dge[dge$gene_name %in% division_genes, ], 'S4_division_gene_expression')
counts <- spire %>% count(event_type, direction) %>% mutate(event_type = factor(event_type, c('EX', 'IR', 'Alt5', 'Alt3')))
totals <- spire %>% count(event_type, name = 'tested')
export(counts, 'SPIRE_event_counts', 'source_data')
robust_summary <- bind_rows(lapply(names(sensitivity), function(ds) {
  x <- sensitivity[[ds]]
  data.frame(dataset = ds, selected = nrow(x), raw_direction_agrees = sum(x$raw_agrees_with_fit),
    leave_one_direction_agrees = sum(x$loo_same_direction),
    leave_one_min_10pp = sum(x$loo_min_signed_pp >= 10), complete_separation = sum(x$complete_separation))
}))
export(robust_summary, 'S3_sensitivity_summary')
print(robust_summary)

theme_set(publication_theme())
update_geom_defaults('text', list(family = 'Liberation Sans', size = 2.5))
figures <- data.frame()
save_figure <- function(p, name, title, height) {
  p <- p + plot_annotation(tag_levels = 'A', theme = publication_theme())
  for (ext in c('pdf', 'png')) ggsave(file.path(out, 'figures', paste0(name, '.', ext)),
    p, width = 180, height = height, units = 'mm', dpi = 450, bg = 'white',
    device = if (ext == 'pdf') grDevices::cairo_pdf else 'png')
  figures <<- rbind(figures, data.frame(figure = name, title = title, width_mm = 180, height_mm = height))
}
direction_scale <- function() scale_color_manual(values = direction_colors,
  breaks = c('down', 'ndiff', 'up'), labels = c('Decreased', 'Unselected', 'Increased'), name = NULL)
pca_plot <- function(name, title) {
  if (name == 'Splicing_PCA') {
    m <- as.matrix(spire[, names(groups$spire)])
    variance <- apply(m, 1, var)
    ix <- which(variance >= 25 & spire$GENE %in% valid_genes(spire$GENE) & !spire$GENE %in% c('Spire1', 'Spire2'))
    ix <- ix[order(-variance[ix], spire$EVENT[ix])]
    ix <- ix[!duplicated(spire$GENE[ix])]
    pc <- prcomp(t(m[ix, ]), center = TRUE, scale. = FALSE)
    x <- data.frame(sample = names(groups$spire), group = unname(groups$spire), pc$x)
    v <- data.frame(component = seq_along(pc$sdev), fraction = pc$sdev^2/sum(pc$sdev^2))
    export(spire[ix, c('GENE', 'EVENT')], 'Splicing_PCA_features', 'source_data')
  } else {
    x <- tab('spire_review', paste0(name, '_coordinates'))
    v <- tab('spire_review', paste0(name, '_variance'))
  }
  export(x, paste0(name, '_coordinates'), 'source_data')
  export(v, paste0(name, '_variance'), 'source_data')
  ggplot(x, aes(PC1, PC2, shape = group, color = group)) + geom_point(size = 2.2) +
    ggrepel::geom_text_repel(aes(label = sub('SPIRE_', '', sample)), size = 2.3,
      seed = 42, max.overlaps = Inf, family = 'Liberation Sans', min.segment.length = 0, show.legend = FALSE) +
    scale_color_manual(values = condition_colors) + scale_shape_manual(values = c(Cont = 1, DKO = 16)) +
    labs(title = title, subtitle = '3 independent biological replicates / condition; Spire1/2 excluded',
      x = sprintf('PC1 (%.1f%%)', 100*v$fraction[1]), y = sprintf('PC2 (%.1f%%)', 100*v$fraction[2]),
      shape = NULL, color = NULL)
}
sample_psi <- function(rows, ds) {
  g <- groups[[ds]]
  x <- rows[, c('GENE', 'EVENT', names(g))] %>% pivot_longer(all_of(names(g)), names_to = 'sample', values_to = 'PSI')
  x$group <- factor(unname(g[x$sample]), unique(unname(g)))
  x$label <- paste(x$GENE, x$EVENT, sep = '\n')
  x
}
psi_points <- function(x, title, subtitle) {
  ggplot(x, aes(group, PSI, color = group)) +
    geom_point(aes(group = sample), position = position_dodge(width = .6), size = 1.8) +
    stat_summary(fun = mean, geom = 'point', shape = 95, size = 6, color = '#111111') +
    facet_wrap(~label, nrow = 1) + scale_color_manual(values = condition_colors, guide = 'none') +
    scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
    labs(x = NULL, y = 'PSI (%)', title = title, subtitle = subtitle)
}

target <- tab('spire_review', 'Target_gene_expression')
target <- target[target$gene_name %in% c('Spire1', 'Spire2'), ]
target$group <- factor(target$group, c('Cont', 'DKO'))
export(target, 'SPIRE_target_expression', 'source_data')
p1a <- ggplot(target, aes(group, normalized_count + 1, color = group)) +
  geom_point(aes(group = sample), position = position_dodge(width = .6), size = 2) +
  facet_wrap(~gene_name, nrow = 1) + scale_y_log10() +
  scale_color_manual(values = condition_colors, guide = 'none') +
  labs(x = NULL, y = 'Normalized counts + 1', title = 'Target-gene expression',
    subtitle = 'Each point is one sample; logarithmic scale')
volcano <- dge[is.finite(dge$padj) & is.finite(dge$log2FoldChange), ]
volcano <- volcano[order(volcano$selected), ]
label_genes <- unique(c('Spire1', 'Spire2', head(dge$gene_name[order(dge$padj)], 3)))
p1b <- ggplot(volcano, aes(log2FoldChange, -log10(pmax(padj, 1e-300)))) +
  geom_point(aes(color = direction), size = .6, alpha = .6) +
  geom_hline(yintercept = -log10(.05), linetype = 'dashed', linewidth = .3, color = '#777777') +
  ggrepel::geom_text_repel(data = volcano[volcano$gene_name %in% label_genes, ],
    aes(label = gene_name), size = 2.3, seed = 42, max.overlaps = Inf, family = 'Liberation Sans') +
  direction_scale() + labs(x = 'Shrunken log2 fold change (DKO / Cont)', y = expression(-log[10]('BH-adjusted P')),
    title = 'Selective expression changes', subtitle = sprintf('%s increased; %s decreased (BH P ≤ 0.05)',
      sum(dge$direction == 'up'), sum(dge$direction == 'down')))
selected_counts <- counts[counts$direction != 'ndiff', ]
p1c <- ggplot(selected_counts, aes(event_type, n, fill = direction)) +
  geom_col(position = position_dodge(.8), width = .65) +
  geom_text(aes(label = n), position = position_dodge(.8), vjust = -.4, size = 2.5) +
  scale_fill_manual(values = direction_colors, breaks = c('down', 'up'), labels = c('Decreased', 'Increased'), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, .15))) +
  labs(x = NULL, y = 'Selected events', title = 'Bidirectional splicing changes',
    subtitle = sprintf('%s events in %s genes; %s events tested', nrow(robust), length(valid_genes(robust$GENE)), format(nrow(spire), big.mark = ',')))
heat <- sample_psi(robust[order(robust$event_type, robust$deltapsi, robust$EVENT), ], 'spire')
heat$centered_PSI <- heat$PSI - ave(heat$PSI, heat$EVENT, FUN = median)
heat$sample <- factor(heat$sample, names(groups$spire))
heat$EVENT <- factor(heat$EVENT, rev(unique(heat$EVENT)))
export(heat, 'SPIRE_selected_PSI_original_order', 'source_data')
p1d <- ggplot(heat, aes(sample, EVENT, fill = centered_PSI)) + geom_raster() +
  scale_fill_gradient2(low = signed_colors[1], mid = signed_colors[2], high = signed_colors[3],
    limits = c(-50, 50), oob = scales::squish, name = 'PSI − median (pp)', breaks = c(-50, 0, 50)) +
  scale_x_discrete(labels = function(s) sub('SPIRE_', '', s)) +
  labs(x = NULL, y = NULL, title = 'Individual-sample splicing profiles', subtitle = 'All selected events; row median centered') +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1))
overlap <- tab('comparison', 'shared_and_dataset_specific_counts')
overlap_long <- pivot_longer(overlap, c(events, genes), names_to = 'unit', values_to = 'n')
overlap_long$category <- factor(overlap_long$category, c('Spire only', 'Shared', 'FMN2 only'))
p2a <- ggplot(overlap_long, aes(category, n, fill = category)) + geom_col(width = .65) +
  geom_text(aes(label = n), vjust = -.4) + facet_wrap(~unit, scales = 'free_y') +
  scale_fill_manual(values = c('Spire only' = dataset_colors[['spire']], Shared = '#222222', 'FMN2 only' = dataset_colors[['fmn2']]), guide = 'none') +
  scale_x_discrete(labels = c('SPIRE\nonly', 'Shared', 'FMN2\nonly')) +
  scale_y_continuous(expand = expansion(mult = c(0, .2))) +
  labs(x = NULL, y = 'Selected', title = 'A restricted set of shared targets',
    subtitle = '48,323 common events in 4,229 genes')
shared <- common[common$spire_selected & common$fmn2_selected, ]
export(shared, 'Convergence_shared_events', 'source_data')
p2b <- ggplot(common, aes(100*spire_raw_mean_deltapsi, 100*fmn2_raw_mean_deltapsi)) +
  geom_bin_2d(bins = 45) + scale_fill_gradient(low = '#E8E8E8', high = '#333333', trans = 'log10', name = 'Events / bin') +
  geom_hline(yintercept = 0, linewidth = .25, color = '#999999') +
  geom_vline(xintercept = 0, linewidth = .25, color = '#999999') +
  geom_point(data = shared, aes(shape = comparison), size = 1.6, color = '#6A3D9A') +
  scale_shape_manual(values = c('Both: same direction' = 1, 'Both: opposite direction' = 4),
    breaks = c('Both: same direction', 'Both: opposite direction'),
    labels = c('Concordant fit', 'Opposing fit'), name = NULL) +
  labs(x = 'SPIRE raw mean ΔPSI (pp)', y = 'FMN2 raw mean ΔPSI (pp)', title = 'No global directional similarity',
    subtitle = sprintf('All common events; Spearman ρ = %.4f', cor(common$spire_raw_mean_deltapsi, common$fmn2_raw_mean_deltapsi, method = 'spearman'))) +
  theme(legend.box = 'vertical')
shared$label <- paste(shared$GENE, shared$EVENT)
shared <- shared[order(shared$comparison, shared$GENE), ]
effect <- bind_rows(lapply(c('spire', 'fmn2'), function(ds)
  data.frame(EVENT = shared$EVENT, label = shared$label, dataset = ds,
    fitted_pp = shared[[paste0(ds, '_deltapsi')]]*100, comparison = shared$comparison)))
effect$label <- factor(effect$label, rev(shared$label))
export(effect, 'Convergence_shared_fitted_effects', 'source_data')
p2c <- ggplot(effect, aes(fitted_pp, label, color = dataset)) +
  geom_vline(xintercept = 0, linewidth = .3, color = '#999999') +
  geom_line(aes(group = EVENT), color = '#BBBBBB', linewidth = .4, orientation = 'y') +
  geom_point(aes(shape = dataset), size = 2) +
  scale_color_manual(values = dataset_colors, labels = c(spire = 'SPIRE', fmn2 = 'FMN2'), name = NULL) +
  scale_shape_manual(values = c(spire = 16, fmn2 = 17), labels = c(spire = 'SPIRE', fmn2 = 'FMN2'), name = NULL) +
  labs(x = 'Fitted ΔPSI (knockout − respective control; pp)', y = NULL,
    title = 'Shared targets do not imply matching direction', subtitle = '7 concordant and 4 opposing events; exploratory sign test P = 0.549')
shared_psi <- bind_rows(lapply(names(events), function(ds) {
  x <- sample_psi(events[[ds]][match(shared$EVENT, events[[ds]]$EVENT), ], ds)
  x$dataset <- if (ds == 'spire') 'SPIRE · 3 samples / group' else 'FMN2 · 2 pools / group'
  x$sample <- if (ds == 'spire') sub('SPIRE_', '', x$sample) else
    unname(setNames(c('Het_1', 'Het_2', 'KO_1', 'KO_2'), names(groups$fmn2))[x$sample])
  x
}))
shared_psi$sample <- factor(shared_psi$sample, c(sub('SPIRE_', '', names(groups$spire)), 'Het_1', 'Het_2', 'KO_1', 'KO_2'))
shared_psi$dataset <- factor(shared_psi$dataset, c('SPIRE · 3 samples / group', 'FMN2 · 2 pools / group'))
shared_psi$label <- factor(gsub('\n', ' ', shared_psi$label), rev(shared$label))
export(shared_psi, 'Convergence_shared_sample_PSI', 'source_data')
p2d <- ggplot(shared_psi, aes(sample, label, fill = PSI)) + geom_tile(color = 'white', linewidth = .2) +
  facet_grid(~dataset, scales = 'free_x', space = 'free_x') +
  scale_fill_viridis_c(option = 'cividis', direction = -1, limits = c(0, 100), name = 'PSI (%)') +
  labs(x = NULL, y = NULL, title = 'Raw measurements reveal heterogeneous responses',
    subtitle = 'All 11 shared events; none of the 7 concordant events separates both studies completely') +
  theme(axis.line = element_blank(), axis.ticks = element_blank())

qc <- tab('spire_review', 'RNAseq_alignment_QC')
qc$group <- unname(groups$spire[qc$Sample])
export(qc, 'S1_alignment_QC')
ps1c <- ggplot(qc, aes(factor(Sample, names(groups$spire)), uniquely_mapped_percent, fill = group)) +
  geom_col(width = .65) + geom_text(aes(label = sprintf('%.1f', uniquely_mapped_percent)), vjust = -.4) +
  scale_fill_manual(values = condition_colors, guide = 'none') + scale_y_continuous(limits = c(0, 100)) +
  scale_x_discrete(labels = function(s) sub('SPIRE_', '', s)) +
  labs(x = NULL, y = 'Uniquely mapped reads (%)', title = 'Alignment quality is similar across samples') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
exons <- tab('spire_review', 'Target_exon_support_all_samples')
exons <- exons[exons$GENE %in% c('Spire1', 'Spire2'), ]
exons$sample <- factor(exons$sample, names(groups$spire))
export(exons, 'S1_target_exon_support')
ps1d <- ggplot(exons, aes(sample, EVENT, fill = log10(corrected_total + 1))) +
  geom_tile() + facet_wrap(~GENE, scales = 'free_y') +
  scale_fill_viridis_c(option = 'cividis', direction = -1, name = 'log10(support + 1)') +
  scale_x_discrete(labels = function(s) sub('SPIRE_', '', s)) +
  labs(x = NULL, y = NULL, title = 'Exon-level loss at the targeted loci', subtitle = 'VAST corrected junction support before N10 filtering') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size = 6), axis.line = element_blank(), axis.ticks = element_blank())

export(robust[, c('GENE', 'EVENT', 'event_type', 'deltapsi', 'direction', 'raw_delta_pp', 'loo_min_pp', 'loo_max_pp', 'loo_same_direction')],
  'SPIRE_sensitivity', 'source_data')
ps2a <- ggplot(robust, aes(100*deltapsi, raw_delta_pp, color = direction)) +
  geom_abline(slope = 1, intercept = 0, linewidth = .3, color = '#999999') +
  geom_hline(yintercept = 0, linewidth = .25, color = '#BBBBBB') +
  geom_point(size = .7, alpha = .6) + direction_scale() +
  labs(x = 'Fitted ΔPSI (pp)', y = 'Raw mean ΔPSI (pp)', title = 'Fitted and observed effects',
    subtitle = 'All 1,547 selected SPIRE events; diagonal is identity')
stability <- robust %>% group_by(event_type) %>% summarise(selected = n(),
  `Same sign after each omission` = mean(loo_same_direction)*100,
  `At least 10 pp after each omission` = mean(loo_min_signed_pp >= 10)*100, .groups = 'drop') %>%
  pivot_longer(c(`Same sign after each omission`, `At least 10 pp after each omission`), names_to = 'criterion', values_to = 'percent')
export(stability, 'SPIRE_stability', 'source_data')
stability$event_type <- factor(stability$event_type, c('EX', 'IR', 'Alt5', 'Alt3'))
ps2b <- ggplot(stability, aes(event_type, percent, fill = criterion)) +
  geom_col(position = position_dodge(.75), width = .65) +
  scale_fill_manual(values = c('#999999', '#222222'), name = NULL) +
  scale_y_continuous(limits = c(0, 100)) + guides(fill = guide_legend(nrow = 2)) +
  labs(x = NULL, y = '% of selected events', title = 'Leave-one-sample sensitivity',
    subtitle = 'Raw PSI means; descriptive, without refitting betAS')
thresholds <- expand.grid(FPR = c(.01, .025, .05), delta_pp = c(10, 15, 20))
thresholds$events <- vapply(seq_len(nrow(thresholds)), function(i)
  sum(spire$FDR <= thresholds$FPR[i] & abs(spire$deltapsi)*100 >= thresholds$delta_pp[i]), integer(1))
export(thresholds, 'S3_threshold_sensitivity')
ps2c <- ggplot(thresholds, aes(factor(FPR), factor(delta_pp), fill = events)) + geom_tile(color = 'white') +
  geom_text(aes(label = events), color = '#111111') + scale_fill_gradient(low = '#FAFAFA', high = '#6A3D9A', guide = 'none') +
  labs(x = 'Empirical FPR cutoff', y = '|Fitted ΔPSI| cutoff (pp)', title = 'Dependence on selection thresholds',
    subtitle = 'FPR is per event, not a multiple-testing FDR') + theme(axis.line = element_blank(), axis.ticks = element_blank())
coverage <- bind_rows(lapply(names(events), function(ds) {
  x <- events[[ds]]
  data.frame(dataset = ds, scope = c('All tested', 'Common tested'),
    tested = c(nrow(x), nrow(common)), selected = c(sum(x$direction != 'ndiff'), sum(common[[paste0(ds, '_selected')]])))
}))
export(coverage, 'S5_tested_universes')
vst <- read_input('data/inputs/spire_fastq/dge/tables/processed_abundance/rnaseq/all.vst.tsv')
expr <- vst[match(dge$gene_id[dge$selected], vst$gene_id), ]
stopifnot(!anyNA(expr$gene_id))
expr$effect <- dge$log2FoldChange[match(expr$gene_id, dge$gene_id)]
expr <- expr[order(expr$effect), ]
expr_long <- pivot_longer(expr[, c('gene_id', names(groups$spire))], all_of(names(groups$spire)), names_to = 'sample', values_to = 'VST')
expr_long$centered <- expr_long$VST - ave(expr_long$VST, expr_long$gene_id, FUN = mean)
expr_long$gene_id <- factor(expr_long$gene_id, rev(expr$gene_id))
expr_long$sample <- factor(expr_long$sample, names(groups$spire))
export(expr_long, 'SPIRE_DEG_VST', 'source_data')
ps3a <- ggplot(expr_long, aes(sample, gene_id, fill = centered)) + geom_raster() +
  scale_fill_gradient2(low = signed_colors[1], mid = signed_colors[2], high = signed_colors[3],
    limits = c(-3, 3), oob = scales::squish, name = 'VST − gene mean') +
  scale_x_discrete(labels = function(s) sub('SPIRE_', '', s)) +
  labs(x = NULL, y = NULL, title = 'Expression across all selected genes', subtitle = '111 genes; ordered by shrunken log2 fold change') +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
joint_genes <- intersect(valid_genes(spire$GENE), valid_genes(dge$gene_name[is.finite(dge$padj)]))
selected_splicing <- intersect(joint_genes, valid_genes(robust$GENE))
selected_expression <- intersect(joint_genes, valid_genes(dge$gene_name[dge$selected]))
joint <- data.frame(GENE = joint_genes, expression_selected = joint_genes %in% selected_expression,
  splicing_selected = joint_genes %in% selected_splicing)
export(joint, 'S4_expression_splicing_gene_membership')
joint_counts <- data.frame(category = c('Expression\nonly', 'Both', 'Splicing\nonly'),
  genes = c(length(setdiff(selected_expression, selected_splicing)), length(intersect(selected_expression, selected_splicing)),
    length(setdiff(selected_splicing, selected_expression))))
joint_counts$category <- factor(joint_counts$category, joint_counts$category)
ps3b <- ggplot(joint_counts, aes(category, genes)) + geom_col(width = .6, fill = '#777777') +
  geom_text(aes(label = genes), vjust = -.4) + scale_y_continuous(expand = expansion(mult = c(0, .15))) +
  labs(x = NULL, y = 'Genes', title = 'Expression and splicing target sets', subtitle = paste(format(length(joint_genes), big.mark = ','), 'genes testable in both assays'))
go_summary <- tab('spire_final', 'GO_analysis_summary')
go <- go_summary[go_summary$gene_set == 'all', ]
go$background <- factor(go$background, c('paper_OrgDb', 'tested_N10'), c('All OrgDb', 'Tested genes'))
ps3c <- ggplot(go, aes(ontology, enriched_terms, fill = background)) +
  geom_col(position = position_dodge(.7), width = .6) +
  geom_text(aes(label = enriched_terms), position = position_dodge(.7), vjust = -.4) +
  scale_fill_manual(values = c('#999999', '#222222'), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, .18))) +
  labs(x = NULL, y = 'Enriched GO terms', title = 'Functional enrichment depends on background',
    subtitle = 'All selected SPIRE genes; no terms with tested-gene background')
division_tests <- bind_rows(lapply(names(division_terms), function(term) {
  universe <- valid_genes(spire$GENE); selected <- valid_genes(robust$GENE)
  annotated <- intersect(universe, division$SYMBOL[division$GOALL == term])
  test <- overlap_test(selected, annotated, universe)
  data.frame(GO = term, term = division_terms[[term]], universe = test$universe,
    selected_genes = test$spire, annotated_tested_genes = test$fmn2,
    selected_annotated_genes = test$shared, expected = test$expected,
    odds_ratio = test$odds_ratio, pvalue = test$p.value)
}))
division_tests$padj <- p.adjust(division_tests$pvalue, 'BH')
export(division_tests, 'S4_legacy_four_term_enrichment')
ps3d <- ggplot(division_tests, aes(-log10(padj), reorder(term, -log10(padj)))) +
  geom_segment(aes(x = 0, xend = -log10(padj), yend = reorder(term, -log10(padj))), color = '#BBBBBB') +
  geom_point(size = 2, color = '#222222') + geom_vline(xintercept = -log10(.05), linetype = 'dashed', linewidth = .3) +
  labs(x = expression(-log[10]('BH-adjusted P')), y = NULL, title = 'Division-category enrichment',
    subtitle = 'Four categories; tested-gene background') + theme(plot.title.position = 'plot')

matt <- bind_rows(lapply(names(studies), function(ds) {
  x <- tab(studies[[ds]], 'Matt_all_feature_Wilcoxon_statistics'); x$dataset <- ds; x
}))
impact <- tab('spire_final', 'protein_impact_differential_events')
impact_counts <- impact %>% count(Impact, direction)
ps4b <- ggplot(impact_counts, aes(n, Impact, fill = direction)) + geom_col() +
  scale_fill_manual(values = direction_colors, labels = c(down = 'Decreased', up = 'Increased'), name = NULL) +
  scale_y_discrete(labels = c('ORF disrupted upon inclusion' = 'ORF disrupted\non inclusion',
    'ORF disrupted upon exclusion' = 'ORF disrupted\non exclusion', 'Unknown/noncoding' = 'Unknown / noncoding')) +
  labs(x = 'Selected SPIRE events', y = NULL, title = 'Predicted transcript consequences', subtitle = 'VASTDB v3 annotation; not measured protein loss')
source('analysis/publication_inference.R')
source('analysis/publication_panels.R')
export(figures, 'figure_manifest')
for (name in c('Target_gene_DGE', 'Target_gene_expression', 'Target_exon_support_summary',
  'sample_correspondence', 'Expression_PCA_features')) export(tab('spire_review', name), paste0('S1_', name))
for (name in c('shared_selected_events', 'shared_selected_genes', 'overlap_enrichment',
  'direction_agreement', 'dataset_comparison_summary', 'common_universe_event_rates'))
  export(tab('comparison', name), paste0('S5_', name))
norm <- read_input('data/inputs/spire_fastq/dge/tables/processed_abundance/rnaseq/all.normalised_counts.tsv')
export(norm, 'S2_normalized_counts')
export(vst, 'S2_VST_expression')
export(matt, 'S7_Matt_all_feature_tests')
export(impact, 'S7_SPIRE_predicted_transcript_impact')
export(tab('fmn2_final', 'protein_impact_differential_events'), 'S7_FMN2_predicted_transcript_impact')
for (name in c('RBP_all_motif_pvalues', 'rMAPS3_shared_RBP_regions', 'rMAPS3_shared_RBP_summary'))
  export(tab('rbp', name), paste0('S8_', name))
for (ds in names(studies)) {
  export(tab(studies[[ds]], 'GO_analysis_summary'), paste0('S6_', ds, '_GO_summary'))
  paths <- list.files(file.path('reports', studies[[ds]], 'tables'), pattern = '^GO_.*[.]rds$', full.names = TRUE)
  for (path in paths) {
    input_files <- union(input_files, path)
    object <- readRDS(path)
    if (!is.null(object)) export(object@result, paste0('S6_', ds, '_', sub('[.]rds$', '', basename(path)), '_all_terms'))
  }
}
for (path in list.files('reports/comparison/tables', pattern = '^GO_common.*[.]csv$', full.names = TRUE))
  export(read_input(path), paste0('S6_comparison_', sub('[.]csv$', '', basename(path))))

summary_dge <- tab('spire_review', 'DGE_summary')
overlap_tests <- tab('comparison', 'overlap_enrichment')
gene_test <- overlap_tests[overlap_tests$unit == 'genes', ]
target_dge <- tab('spire_review', 'Target_gene_DGE')
division_non_target <- division_selected[!division_selected$GENE %in% c('Spire1', 'Spire2'), ]
for (name in c('RESULTS', 'STORY', 'INFERENCE_REVIEW')) {
  path <- paste0('docs/PUBLICATION_', name, '.md')
  input_files <- union(input_files, path)
  file.copy(path, file.path(out, 'text', paste0(name, '.md')))
}
for (name in c('METHODS', 'LEGENDS')) {
  path <- paste0('docs/PUBLICATION_', name, '.md')
  input_files <- c(input_files, path)
  lines <- readLines(path)
  lines <- gsub('[CLUSTER.md](CLUSTER.md)', 'the repository CLUSTER.md', lines, fixed = TRUE)
  lines <- gsub('[SOURCES.md](../analysis/SOURCES.md)', 'the repository SOURCES.md', lines, fixed = TRUE)
  lines <- gsub('[RMAPS3.md](../analysis/RMAPS3.md)', 'the repository RMAPS3.md', lines, fixed = TRUE)
  writeLines(lines, file.path(out, 'text', paste0(name, '.md')))
}
writeLines(c('# Supplementary tables and source data', '',
  'CSV files use comma separators, UTF-8, a header row and NA for missing values. Open them as data, preserving event/gene IDs as text. Tables are grouped by S1–S8 prefixes; split files are parts of the same logical table. Full tested universes and nonsignificant results are retained.', '',
  '| Table | Contents |', '|---|---|',
  '| S1 | Sample groups, saved DESeq2 design, inference-family ledger, alignment QC, target loci, historical-label correspondence and PCA feature IDs |',
  '| S2 | All and selected gene-expression statistics; normalized counts and VST values |',
  '| S3 | Every tested SPIRE event, selected events with raw-PSI sensitivity, and threshold sensitivity |',
  '| S4 | Division GO membership/evidence with unique-alias rescue, all mapped division events, legacy four-term reference and expression/splicing intersection |',
  '| S5 | FMN2 universe, common universe, shared targets, overlap tests and event-count-stratified permutation sensitivity |',
  '| S6 | Revised primary and reference GO, expression GO, all selection sensitivities, mapping/alias audits, event-opportunity matched tests/nulls, plus original GO objects exported as reference |',
  '| S7 | All Matt feature tests, individual feature measurements, a sequence-level GC audit, and predicted transcript consequences |',
  '| S8 | Primary exact regional RBP-presence counts/tests with full-family BH, plus original motif-window minima and shared-label summaries as reference only |', '',
  '`source_data/` contains individual plotted measurements, division-event annotations, feature effects, GC profiles and RBP matrices. The expression MA measurements come from S2; counts and effects can be reconstructed from S3/S5; functional and regional sequence panels use S6–S8. Current inference is identified by `q_BH`; original ontology-specific `p.adjust` and minimum-window raw P columns are reference outputs. All source input files and generated table checksums are listed separately.', '',
  '## Field definitions', '',
  '- `FDR`, `spire_FDR`, `fmn2_FDR`: upstream betAS empirical per-event FPR, **not** BH-adjusted FDR. `invertedFDR` is 1−FPR. Zero has finite simulation resolution.',
  '- `deltapsi`, `spire_deltapsi`, `fmn2_deltapsi`, `betasPsiA/B`: fractions on 0–1 (signed for delta). Delta is knockout minus respective control. Multiply by 100 for percentage points.',
  '- `PSI`, sample columns beginning SPIRE_/SRR, and fields ending `_pp`: percentages or percentage-point differences on a 0–100 scale, signed for differences.',
  '- `Pdiff`, `Fstat`: imported ancillary betAS statistics; neither is an additional publication selection filter.',
  '- `direction`: up/down if the stated FPR and effect criteria are met; ndiff means unselected, not demonstrated absence of change. `event_type`: EX (including ANN), IR, Alt5 or Alt3. `COMPLEX` retains the VAST annotation class.',
  '- `EVENT`: VASTDB event ID, including alternative-site suffixes. `GENE`/`gene_name`: mouse symbol; `gene_id`: Ensembl gene ID. Missing names do not remove events from event counts.',
  '- `COORD`, `FullCO`: original VAST GRCm38/mm10 coordinate annotations. Do not mix with GRCm39 annotations or assume they are BED zero-based intervals.',
  '- `baseMean`: mean normalized expression. `log2FoldChange` and `lfcSE`: imported shrunken effect and accompanying uncertainty; `pvalue`/`padj` in S2 are original DESeq2 Wald P/BH-adjusted P. Missing padj is not set to one.',
  '- `raw_delta_pp`: knockout mean PSI minus control mean PSI. `raw_agrees_with_fit`: strict sign agreement. `loo_same_direction`: every one-sample omission agrees with the fitted sign. `loo_min_signed_pp`: smallest omission effect after orienting to the fitted sign; may be negative. `loo_min_pp`/`loo_max_pp`: signed range of omission effects, not a confidence interval.',
  '- `complete_separation`: strict separation of every knockout/control PSI in the fitted direction. `selected`: uses the documented expression or splicing criterion for that table.',
  '- GO `p.adjust`: BH within the corresponding ontology/background/gene-set analysis; `qvalue`: exported q-value; `geneID`: contributing genes. Legacy focused `padj` is the historical four-term reference only. Current GO `q_BH` pools BP+CC+MF in the stated study/background/selection; `family_tests` gives its denominator. GOALL/EVIDENCEALL preserve inherited membership and evidence codes.',
  '- Overlap legacy `p.BH` covers the original seven tests. Primary S5_overlap_inference_BH `q_BH` covers all nine Fisher/direction/matched-gene tests. Permutation `pvalue` is (1+exceedances)/(1+iterations), not zero.',
  '- Matt `group1/group2`: comparison categories; `median1/median2`, `n1/n2`: actual group medians and event counts; `p.BH`: existing feature-family BH correction. Lengths are nucleotides.',
  '- Primary S8 regional-presence `pvalue` is an exact Fisher/hypergeometric tail; `q_BH` covers all 5278 tests per study. `positive/n` and `background_positive/background_n` are motif-positive/total event counts. `log2_OR` adds 0.5 only for display; tests use integer counts. Original minimum-window `pvalue` fields are reference outputs, not the primary test.',
  '- GO `odds_ratio`, `CI_low/CI_high` are conditional odds ratios and pointwise two-sided 95% exact intervals; they are not simultaneous FDR intervals. The main panel displays observed and expected counts instead.',
  '- `centered_PSI`: within-event median-centered PSI; `centered`: within-gene mean-centered VST. Plot saturation does not alter exported values.', ''), file.path(out, 'text/TABLES.md'))
writeLines(capture.output(sessionInfo()), file.path(out, 'sessionInfo.txt'))
export(AnnotationDbi::metadata(org.Mm.eg.db), 'annotation_metadata')
parameter_paths <- c('data/inputs/spire_fastq/dge/pipeline_info/params_2026-10-29_21-17-23.json',
  'data/inputs/spire_fastq/rnaseq/pipeline_info/params_2026-10-29_19-27-02.json',
  'data/inputs/spire_fastq/dge/pipeline_info/collated_versions.yml',
  'data/inputs/spire_fastq/dge/other/deseq2/rnaseq/DKO_vs_Cont_SPIRE_DKO.R_sessionInfo.log')
input_files <- unique(c(input_files, parameter_paths, 'analysis/publication.R', 'analysis/publication_helpers.R',
  'analysis/publication_panels.R', 'analysis/publication_inference.R', 'analysis/rbp_helpers.R', 'analysis/common.R', 'analysis/palettes.R', 'analysis/comparison_helpers.R', 'analysis/renv.lock',
  'analysis/package_publication.py', 'scripts/run_local_R.sh', 'renv.lock'))
export(data.frame(file = input_files, sha256 = vapply(input_files, function(path)
  digest::digest(file = path, algo = 'sha256'), character(1))), 'input_checksums')
stopifnot(nrow(robust) == sum(selected_counts$n), sum(counts$n) == nrow(spire),
  sum(dge$selected) == summary_dge$padj005, nrow(shared) == 11,
  all(division_tests$padj > .05), all(robust_summary$leave_one_min_10pp <= robust_summary$leave_one_direction_agrees),
  nrow(figures) == 16, all(file.exists(file.path(out, 'figures', paste0(figures$figure, '.pdf')))))
writeLines('Validated publication assembly; figures and manuscript tables generated.', file.path(out, 'COMPLETE'))
cat('Publication assembly complete:', out, '\n')
