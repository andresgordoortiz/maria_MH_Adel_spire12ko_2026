#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
if (length(args) != 3) stop('Usage: scripts/run_R.sh analysis/compare.R SPIRE_ANALYSIS FMN2_ANALYSIS NEW_OUTPUT_DIR')
.libPaths(c(normalizePath('analysis/renv/library'), .libPaths()))
source('analysis/common.R')
source('analysis/comparison_helpers.R')
out <- args[3]
if (file.exists(out)) stop('Comparison output already exists: ', out)
studies <- setNames(args[1:2], c('spire', 'fmn2'))
tables <- params <- list()
for (ds in names(studies)) {
  p <- studies[[ds]]
  stopifnot(file.exists(file.path(p, 'COMPLETE')))
  params[[ds]] <- jsonlite::fromJSON(file.path(p, 'tables', 'analysis_parameters.json'))
  stopifnot(params[[ds]]$dataset == ds, !params[[ds]]$legacy,
    params[[ds]]$min_reads_all_samples == 10, params[[ds]]$fpr_max == .05, params[[ds]]$abs_deltapsi_min == .1)
  tables[[ds]] <- as.data.frame(fread(file.path(p, 'tables', 'all_events.csv')))
}
dir.create(file.path(out, 'tables'), recursive = TRUE)
dir.create(file.path(out, 'figures'))
figure_names <- character(); main_plots <- list()
provenance_label <- 'Final betAS runs; FPR ≤0.05 and |ΔPSI| ≥0.10; common tested events'
colors <- c('Both: same direction' = '#268A6C', 'Both: opposite direction' = '#985277',
  'Spire only' = '#706993', 'FMN2 only' = '#70A0AF', Neither = '#BDBDBD')
common <- join_common_events(tables$spire, tables$fmn2)
common$comparison <- factor(common$comparison, levels = names(colors))
for (ds in names(studies)) {
  x <- tables[[ds]][match(common$EVENT, tables[[ds]]$EVENT), ]
  groups <- unlist(params[[ds]]$samples)
  a <- as.matrix(x[, names(groups)[groups == if (ds == 'spire') 'Cont' else 'Het'], drop = FALSE])
  b <- as.matrix(x[, names(groups)[groups == if (ds == 'spire') 'DKO' else 'KO'], drop = FALSE])
  common[[paste0(ds, '_raw_mean_deltapsi')]] <- (rowMeans(b) - rowMeans(a)) / 100
  common[[paste0(ds, '_complete_group_separation')]] <- ifelse(common[[paste0(ds, '_deltapsi')]] > 0,
    apply(b, 1, min) > apply(a, 1, max), apply(b, 1, max) < apply(a, 1, min))
}
write_table(common, 'common_tested_events')
both <- common$spire_selected & common$fmn2_selected
same <- both & common$spire_deltapsi * common$fmn2_deltapsi > 0
opposite <- both & !same
gene_universe <- valid_genes(common$GENE)
gene_sets <- list(spire = valid_genes(common$GENE[common$spire_selected]),
  fmn2 = valid_genes(common$GENE[common$fmn2_selected]))
shared_genes <- intersect(gene_sets$spire, gene_sets$fmn2)
candidates <- common[both, ]
candidates$min_abs_deltapsi <- pmin(abs(candidates$spire_deltapsi), abs(candidates$fmn2_deltapsi))
candidates <- candidates[order(candidates$comparison, -candidates$min_abs_deltapsi, candidates$EVENT), ]
write_table(candidates, 'shared_selected_events')
write_table(candidates[candidates$comparison == 'Both: same direction', ], 'concordant_event_candidates')
write_table(candidates[candidates$comparison == 'Both: opposite direction', ], 'opposing_event_candidates')
gene_summary <- bind_rows(lapply(shared_genes, function(g) {
  x <- common[common$GENE == g, ]
  data.frame(GENE = g, spire_events = sum(x$spire_selected), fmn2_events = sum(x$fmn2_selected),
    exact_shared_events = sum(x$spire_selected & x$fmn2_selected),
    concordant_events = sum(x$comparison == 'Both: same direction'),
    opposing_events = sum(x$comparison == 'Both: opposite direction'),
    spire_event_ids = paste(x$EVENT[x$spire_selected], collapse = ';'),
    fmn2_event_ids = paste(x$EVENT[x$fmn2_selected], collapse = ';'))
}))
write_table(gene_summary, 'shared_selected_genes')
write_table(data.frame(GENE = gene_universe), 'common_tested_gene_universe')
overlaps <- list(events = overlap_test(common$EVENT[common$spire_selected], common$EVENT[common$fmn2_selected], common$EVENT),
  event_loci = overlap_test(common$locus[common$spire_selected], common$locus[common$fmn2_selected], common$locus),
  genes = overlap_test(gene_sets$spire, gene_sets$fmn2, gene_universe))
for (kind in levels(factor(common$event_type))) {
  x <- common[common$event_type == kind, ]
  overlaps[[kind]] <- overlap_test(x$EVENT[x$spire_selected], x$EVENT[x$fmn2_selected], x$EVENT)
}
overlap_summary <- bind_rows(overlaps, .id = 'unit')
overlap_summary$p.BH <- p.adjust(overlap_summary$p.value, 'BH')
write_table(overlap_summary, 'overlap_enrichment')
direction_test <- if (sum(both)) binom.test(sum(same), sum(both), p = .5) else NULL
direction_summary <- data.frame(shared_events = sum(both), same_direction = sum(same), opposite_direction = sum(opposite),
  same_direction_fraction = if (sum(both)) mean(same[both]) else NA_real_,
  sign_test_p_vs_0_5 = if (is.null(direction_test)) NA_real_ else direction_test$p.value,
  Spearman_all_common = cor(common$spire_deltapsi, common$fmn2_deltapsi, method = 'spearman'),
  Spearman_selected_either = cor(common$spire_deltapsi[common$spire_selected | common$fmn2_selected],
    common$fmn2_deltapsi[common$spire_selected | common$fmn2_selected], method = 'spearman'),
  Spearman_raw_mean_delta = cor(common$spire_raw_mean_deltapsi, common$fmn2_raw_mean_deltapsi, method = 'spearman'))
write_table(direction_summary, 'direction_agreement')
summary <- bind_rows(lapply(names(tables), function(ds) {
  x <- tables[[ds]]; selected <- classify(x$FDR, x$deltapsi) != 'ndiff'
  data.frame(dataset = ds, tested_events = nrow(x), selected_events = sum(selected),
    selected_genes = length(valid_genes(x$GENE[selected])), common_tested_events = nrow(common),
    selected_in_common = sum(common[[paste0(ds, '_selected')]]), common_selected_genes = length(gene_sets[[ds]]),
    control = if (ds == 'spire') 'Cont' else 'Het', profiles_per_group = if (ds == 'spire') 3L else 2L)
}))
write_table(summary, 'dataset_comparison_summary')
counts <- bind_rows(lapply(names(studies), function(ds) {
  common %>% group_by(event_type) %>% summarise(tested = n(), selected = sum(.data[[paste0(ds, '_selected')]]), .groups = 'drop') %>%
    mutate(dataset = ds, selected_per_1000 = 1000 * selected / tested)
}))
write_table(counts, 'common_universe_event_rates')
save_plot(ggplot(counts, aes(event_type, selected_per_1000, fill = dataset)) + geom_col(position = 'dodge') +
  scale_fill_manual(values = c(spire = '#706993', fmn2 = '#70A0AF')) + labs(x = 'Event type', y = 'Selected per 1,000 commonly tested events',
    title = 'Differential splicing with matched detection opportunities'), 'Comparison_event_rates', 8, 5)
scatter <- common[order(common$comparison != 'Neither'), ]
label_rows <- candidates %>% arrange(desc(min_abs_deltapsi)) %>% distinct(GENE, .keep_all = TRUE) %>% slice_head(n = 8)
save_plot(ggplot(scatter, aes(spire_deltapsi, fmn2_deltapsi)) + geom_point(aes(color = comparison), size = .6, alpha = .5) +
  geom_vline(xintercept = c(-.1, .1), linetype = 'dashed', linewidth = .3) +
  geom_hline(yintercept = c(-.1, .1), linetype = 'dashed', linewidth = .3) + scale_color_manual(values = colors) +
  ggrepel::geom_text_repel(data = label_rows, aes(label = GENE), size = 3, seed = 42, max.overlaps = Inf) +
  labs(x = 'Spire ΔPSI (DKO − Cont)', y = 'FMN2 ΔPSI (KO − Het)', color = NULL,
    title = 'Shared events can change in the same or opposite direction'), 'Comparison_delta_PSI', 9, 7)
bars <- data.frame(category = c('Spire only', 'Shared', 'FMN2 only'),
  events = c(sum(common$spire_selected & !common$fmn2_selected), sum(both), sum(common$fmn2_selected & !common$spire_selected)),
  genes = c(length(setdiff(gene_sets$spire, gene_sets$fmn2)), length(shared_genes), length(setdiff(gene_sets$fmn2, gene_sets$spire))))
write_table(bars, 'shared_and_dataset_specific_counts')
save_plot(ggplot(pivot_longer(bars, c(events, genes), names_to = 'unit', values_to = 'n'), aes(category, n, fill = category)) +
  geom_col(width = .65) + geom_text(aes(label = n), vjust = -.4) + facet_wrap(~unit, scales = 'free_y') +
  scale_fill_manual(values = c('Spire only' = '#706993', Shared = '#268A6C', 'FMN2 only' = '#70A0AF')) +
  labs(x = NULL, y = 'Count', title = 'Shared affected events and genes within the common universe'), 'Comparison_overlap_counts', 9, 5)
if (nrow(candidates)) {
  psi <- bind_rows(lapply(names(tables), function(ds) {
    x <- tables[[ds]][match(candidates$EVENT, tables[[ds]]$EVENT), ]
    groups <- unlist(params[[ds]]$samples)
    x <- x[, c('GENE', 'EVENT', names(groups))]
    x <- pivot_longer(x, all_of(names(groups)), names_to = 'sample', values_to = 'PSI')
    labels <- setNames(paste0(groups, '_', ave(seq_along(groups), groups, FUN = seq_along)), names(groups))
    x$sample <- factor(paste(ds, labels[x$sample]), levels = paste(ds, labels))
    x$dataset <- ds; x
  }))
  psi$event_label <- factor(paste(psi$GENE, psi$EVENT), levels = rev(unique(paste(candidates$GENE, candidates$EVENT))))
  save_plot(ggplot(psi, aes(sample, event_label, fill = PSI)) + geom_tile() +
    facet_grid(~dataset, scales = 'free_x', space = 'free_x') + scale_fill_viridis_c(limits = c(0, 100)) +
    labs(x = NULL, y = NULL, title = 'Raw sample PSI for events selected in both studies') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)), 'Comparison_shared_event_sample_PSI', 13, max(6, nrow(candidates) * .3))
  write_table(psi, 'shared_event_sample_PSI')
}

source('analysis/compare_go.R')
source('analysis/compare_features.R')
writeLines(capture.output(sessionInfo()), file.path(out, 'sessionInfo.txt'))
sources <- c(list.files('analysis', pattern = '[.](R|py|sh|md|tsv)$', full.names = TRUE), 'analysis/renv.lock',
  unlist(lapply(studies, function(p) file.path(p, c('tables/all_events.csv', 'tables/analysis_parameters.json', 'COMPLETE')))))
write_table(data.frame(file = sources, sha256 = vapply(sources, function(p) digest::digest(file = p, algo = 'sha256'), character(1))), 'source_checksums')
report <- c('# Spire versus FMN2: final splicing comparison', '',
  sprintf('The common tested universe contains %s events in %s genes. Spire selects %s of these events; FMN2 selects %s.',
    nrow(common), length(gene_universe), sum(common$spire_selected), sum(common$fmn2_selected)), '',
  sprintf('Shared: %s event IDs at %s event loci; %s affected genes. Of the shared event IDs, %s change in the same direction and %s in opposite directions.',
    sum(both), length(unique(common$locus[both])), length(shared_genes), sum(same), sum(opposite)), '',
  sprintf('Gene overlap: %.2f expected by chance under a uniform common-gene universe; one-sided Fisher p = %.3g (BH across reported overlap tests = %.3g).',
    overlaps$genes$expected, overlaps$genes$p.value, overlap_summary$p.BH[overlap_summary$unit == 'genes']), '',
  sprintf('Direction agreement = %.1f%%; two-sided sign test against 50%%: p = %.3g. Shared targeting and concordant direction are separate findings.',
    100 * direction_summary$same_direction_fraction, direction_summary$sign_test_p_vs_0_5), '',
  paste('Concordant event genes:', paste(valid_genes(common$GENE[same]), collapse = ', ')), '',
  paste('Shared genes without the same selected event:', paste(gene_summary$GENE[gene_summary$exact_shared_events == 0], collapse = ', ')), '',
  'Highest shared concordant changes (percentage points):', '',
  vapply(seq_len(min(8, sum(candidates$comparison == 'Both: same direction'))), function(i) {
    x <- candidates[candidates$comparison == 'Both: same direction', ][i, ]
    sprintf('- %s (%s): Spire %+.1f pp; FMN2 %+.1f pp; complete raw-PSI group separation in both: %s.',
      x$GENE, x$EVENT, 100 * x$spire_deltapsi, 100 * x$fmn2_deltapsi,
      x$spire_complete_group_separation && x$fmn2_complete_group_separation)
  }, character(1)), '',
  'See tables/overlap_enrichment.csv, concordant_event_candidates.csv, opposing_event_candidates.csv and shared_selected_genes.csv for every candidate and test.', '',
  sprintf('GO terms significant in both separate study gene sets on the common universe: %s. Shared-gene-set GO: %s BP, %s CC, %s MF terms.',
    nrow(shared_terms), nrow(go_selected$shared_genes_BP), nrow(go_selected$shared_genes_CC), nrow(go_selected$shared_genes_MF)), '',
  sprintf('Matt feature comparisons versus ndiff with the same median-shift direction and significance in both: %s using raw p, %s using BH-adjusted p. Background feature distributions remain dataset-specific.',
    sum(joined_features$group2 == 'ndiff' & joined_features$same_effect_direction & joined_features$raw_significant_both, na.rm = TRUE),
    sum(joined_features$group2 == 'ndiff' & joined_features$same_effect_direction & joined_features$BH_significant_both, na.rm = TRUE)), '',
  'GO uses the same common-tested Entrez universe for both studies and for their shared-gene sets. A shared GO label alone does not show that the same genes or events drive enrichment. Complete results and overlapping term/driver tables are exported.', '',
  'Cautions: FMN2 has two biological pools per group with heterozygous controls; Spire has three profiles per group with animal/pooling identity unverified. These are separate experiments, not a factorial design. A missing selection is not proof of no effect. betAS FDR is an empirical per-event FPR, not a multiplicity-adjusted FDR. Complementary alternative-site rows and multiple events within genes are dependent; event/locus Fisher and sign tests are exploratory. Gene overlap uses a uniform gene null and does not fully correct gene length or event-count bias.', '',
  'RBP motifs still require the external rMAPS2 run. No motifs or gene-expression results are inferred here.')
writeLines(report, file.path(out, 'REPORT.md'))
html <- c('<!doctype html><html lang="en"><meta charset="utf-8"><title>Spire and FMN2 comparison</title>',
  '<style>body{font:16px system-ui;max-width:1100px;margin:30px auto;padding:0 20px}img{max-width:100%}article{border-top:1px solid #ccc;margin:30px 0}</style>',
  '<h1>Spire and FMN2: final splicing comparison</h1>',
  sprintf('<p><strong>%s shared selected event IDs; %s same direction, %s opposite; %s shared affected genes.</strong></p>', sum(both), sum(same), sum(opposite), length(shared_genes)),
  '<p>All overlap calculations use only events tested in both studies. Shared targets do not imply shared direction, mechanism or phenotype.</p>',
  '<p><a href="REPORT.md">Report and limitations</a> · <a href="tables/concordant_event_candidates.csv">Concordant candidates</a> · <a href="tables/shared_selected_genes.csv">Shared genes</a> · <a href="tables/overlap_enrichment.csv">Overlap tests</a></p>')
for (name in unique(figure_names)) html <- c(html, sprintf('<article><h2>%s</h2><a href="figures/%s.pdf">PDF</a><img loading="lazy" src="figures/%s.png" alt="%s"></article>', name, name, name, name))
writeLines(c(html, '</html>'), file.path(out, 'index.html'))
if (length(unique(dirname(c(unname(studies), out)))) == 1L) {
  links <- c(Spire = basename(studies[['spire']]), FMN2 = basename(studies[['fmn2']]), Comparison = basename(out))
  writeLines(c('<!doctype html><html lang="en"><meta charset="utf-8"><title>Final splicing analyses</title>',
    '<style>body{font:18px system-ui;max-width:800px;margin:40px auto;padding:0 20px}</style><h1>Final splicing analyses</h1>',
    sprintf('<p>%s shared selected events: %s concordant, %s opposing; %s shared affected genes.</p>', sum(both), sum(same), sum(opposite), length(shared_genes)),
    sprintf('<p><a href="%s/index.html">%s figures and report</a></p>', links, names(links)),
    '<p>Comparison uses the common tested universe. FMN2 controls are heterozygous; biological replication and coverage differ. RBP motifs require the external rMAPS2 analysis.</p></html>'),
    file.path(dirname(out), 'index.html'))
}
stopifnot(nrow(common) <= min(vapply(tables, nrow, integer(1))), sum(same) + sum(opposite) == sum(both))
writeLines('Final comparison complete; external RBP analysis remains separate.', file.path(out, 'COMPLETE'))
cat(length(unique(figure_names)), 'comparison figures generated in', out, '\n')
