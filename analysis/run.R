#!/usr/bin/env Rscript
source('analysis/load.R')
figure_names <- character()
main_plots <- list()
source('analysis/splicing.R')
status <- data.frame(step = 'core_splicing', status = 'complete')
annotation_path <- file.path('data/reference', 'event_annotations.tsv')
if (file.exists(annotation_path)) {
  source('analysis/annotations.R')
  source('analysis/external.R')
} else status <- rbind(status, data.frame(step = 'annotations_features_GO_RBP_inputs',
  status = 'pending: run analysis/prepare_references.py then rerun into a new directory'))
if (requireNamespace('patchwork', quietly = TRUE)) {
  names2 <- c('Fig2C_splicing_volcano', 'Fig2D_event_counts', 'Fig2F_EX_miniheatmap',
    'Fig2F_IR_miniheatmap', 'Fig2G_Matt_EX', 'Fig2G_Matt_IR')
  names2 <- intersect(names2, names(main_plots))
  combined <- patchwork::wrap_plots(main_plots[names2], ncol = 2) +
    patchwork::plot_annotation(title = paste(dataset_title, ': splicing landscape'), tag_levels = 'A')
  save_plot(combined, paste0('Figure2_', dataset, '_splicing_panels'), 14, if (length(names2) > 4) 15 else 10)
  names3 <- intersect(c('Fig3A_protein_impact', 'GO_paper_OrgDb_all_BP'), names(main_plots))
  if (length(names3) == 2) save_plot(patchwork::wrap_plots(main_plots[names3], ncol = 1) +
    patchwork::plot_annotation(title = paste(dataset_title, ': predicted functional impact'), tag_levels = 'A'),
    paste0('Figure3_', dataset, '_functional_panels'), 12, 14)
}
write_table(status, 'analysis_status')
source_files <- unique(c(source_files, list.files('analysis', pattern = '[.](R|py|sh|md|tsv)$', full.names = TRUE),
  'config/analysis.R', '.R-version', 'renv.lock', 'analysis/renv.lock',
  list.files('data/reference', pattern = 'sources.json$', full.names = TRUE)))
reference_sources <- list.files('data/reference', pattern = 'sources.json$', full.names = TRUE)
for (path in reference_sources) file.copy(path, file.path(out, 'tables', basename(path)))
file.copy('analysis/renv.lock', file.path(out, 'analysis_renv.lock'))
jsonlite::write_json(list(dataset = dataset, title = dataset_title, source = args[1], min_reads_all_samples = 10L,
  samples = as.list(config$sample_groups), fpr_max = .05, abs_deltapsi_min = .1,
  paper_commit = 'e231983fd82736971aef7fc938cf34bed5ce3121',
  legacy = args[1] == 'legacy', input_provenance = input$provenance),
  file.path(out, 'tables', 'analysis_parameters.json'), auto_unbox = TRUE, pretty = TRUE)
write_table(data.frame(file = source_files, sha256 = vapply(source_files,
  function(p) digest::digest(file = p, algo = 'sha256'), character(1))), 'source_checksums')
writeLines(capture.output(sessionInfo()), file.path(out, 'sessionInfo.txt'))
figs <- unique(figure_names)
overview <- grep('^Figure[23]_', figs, value = TRUE)
figs <- c(overview, setdiff(figs, overview))
html <- c(paste0('<!doctype html><html lang="en"><meta charset="utf-8"><title>', dataset_title, ' splicing analysis</title>'),
  '<style>body{font:16px system-ui;max-width:1100px;margin:30px auto;padding:0 20px}img{width:100%}article{border-top:1px solid #ccc;margin:30px 0}a{color:#514a77}</style>',
  paste0('<h1>', dataset_title, ' splicing</h1>'), paste0('<p><strong>', provenance_label, '</strong></p>'),
  sprintf('<p>N10: %s events; selected: %s events in %s genes. FPR ≤0.05, |ΔPSI| ≥0.10. Positive ΔPSI = greater inclusion/site use/retention in the KO group.</p>',
    nrow(events), sum(events$direction != 'ndiff'), length(valid_genes(events$GENE[events$direction != 'ndiff']))),
  '<p>betAS FDR is an empirical per-event FPR, not a BH-adjusted p-value. See REPORT.md, analysis/README.md and tables/analysis_status.csv for methods and pending panels.</p>',
  '<p><a href="tables/differential_events.csv">Differential events</a> · <a href="tables/analysis_status.csv">Panel status</a> · <a href="tables/source_checksums.csv">Source checksums</a></p>')
if (exists('go_results')) {
  gs <- bind_rows(go_results)
  bp <- gs[gs$gene_set == 'all' & gs$ontology == 'BP', ]
  html <- c(html, sprintf('<p><strong>GO background sensitivity:</strong> OrgDb background: %s enriched BP terms; tested-N10 background: %s. Interpret biological themes with this dependence in mind.</p>',
    bp$enriched_terms[bp$background == 'paper_OrgDb'], bp$enriched_terms[bp$background == 'tested_N10']))
}
if (exists('matched')) html <- c(html, sprintf('<p>VASTDB sequence annotations missing: %s tested events, including %s selected events. They remain in splicing results. Protein-impact annotation version: v3.</p>',
  sum(!matched), sum(!matched & events$direction != 'ndiff')))
if (exists('changed_coords')) html <- c(html, sprintf('<p>Reference coordinates changed for %s additional IDs; these are excluded from sequence and external-tool inputs and recorded separately.</p>', sum(changed_coords)))
html <- c(html, '<details><summary>Figure list</summary><ul>',
  paste0('<li><a href="#', figs, '">', figs, '</a></li>'), '</ul></details>')
for (name in figs) html <- c(html, sprintf('<article id="%s"><h2>%s</h2><a href="figures/%s.pdf">PDF</a><img loading="lazy" src="figures/%s.png" alt="%s"></article>', name, name, name, name, name))
writeLines(c(html, '</html>'), file.path(out, 'index.html'))
report <- c(paste('#', dataset_title, 'splicing'), '', provenance_label, '',
  sprintf('%s profiles (%s %s, %s %s). N10 retains %s events. The paper criteria select %s events in %s genes.',
    length(samples), length(group_a_columns), config$group_a, length(group_b_columns), config$group_b,
    nrow(events), sum(events$direction != 'ndiff'), length(valid_genes(events$GENE[events$direction != 'ndiff']))), '',
  sprintf('%s selected events lack gene labels; they remain in event-level analyses.',
    sum(events$direction != 'ndiff' & !events$GENE %in% valid_genes(events$GENE))), '',
  if (dataset != 'fmn2') 'Independent animals, pooling and batches remain unverified.' else
    'FMN2 controls are heterozygous pools, not wild type; two biological pools per group, with technical runs merged.', '',
  'Open [the figure gallery](index.html). PDFs and PNGs are in figures/; plotted data and full statistical tables are in tables/.', '',
  'betAS FDR is an empirical per-event FPR, not a BH-adjusted FDR. Small fitted differences at PSI boundaries can receive zero FPR; the paper effect-size filter is applied, but does not calibrate that statistic.', '',
  'The source code, sample mapping, parameters, source hashes and R package versions are recorded. Final exports are checked against local inclusion-table content and the completed run checkpoints. Different serialized RDS hashes are recorded explicitly.', '')
if (exists('go_results')) report <- c(report,
  sprintf('GO background sensitivity: %s BP terms with the paper OrgDb universe versus %s with the tested-N10 universe. The two analyses are exported separately.',
    bp$enriched_terms[bp$background == 'paper_OrgDb'], bp$enriched_terms[bp$background == 'tested_N10']), '')
if (exists('matched')) report <- c(report,
  sprintf('VASTDB EVENT_INFO lacks %s tested IDs (%s selected). An additional %s IDs have changed coordinates and are excluded from sequence/Matt/RBP inputs. They remain in splicing/GO results. Protein-impact annotations use available mm10 v3 rather than the unavailable v2.3.',
    sum(!matched), sum(!matched & events$direction != 'ndiff'), sum(changed_coords)), '')
report <- c(report, 'Panel status:', '', paste0('- ', status$step, ': ', status$status), '',
  'The companion comparison report evaluates shared events, directions, genes and GO using a common tested universe. Gene expression and experimental imaging/meiosis panels are outside the splicing-only analysis. RBP motif analysis can run locally with rMAPS3; see analysis/RMAPS3.md and the panel status above. No upload was performed.')
writeLines(report, file.path(out, 'REPORT.md'))
stopifnot(all(file.exists(file.path(out, 'figures', paste0(figs, '.pdf')))),
  sum(table(events$direction)) == nrow(events))
writeLines(provenance_label, file.path(out, 'COMPLETE'))
cat(length(figs), 'figures generated in', out, '\n')
