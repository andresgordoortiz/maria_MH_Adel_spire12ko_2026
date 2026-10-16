#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop('Usage: scripts/run_R.sh analysis/run.R legacy|RUN/csv/spire_all_events.csv OUTPUT_DIR')
.libPaths(c(normalizePath('analysis/renv/library', mustWork = FALSE), .libPaths()))
source('analysis/common.R')
source('config/analysis.R')
stopifnot(dataset == 'spire', as.character(getRversion()) == trimws(readLines('.R-version')))
out <- args[2]
if (file.exists(out)) stop('Output already exists; choose a new directory: ', out)
input <- readRDS(config$filtered)
samples <- names(config$sample_groups)
stopifnot(identical(names(input$PSI)[7:12], samples), !anyDuplicated(input$PSI$EVENT),
  nrow(input$PSI) > 0, all(is.finite(as.matrix(input$PSI[, samples]))),
  all(as.matrix(input$PSI[, samples]) >= 0), all(as.matrix(input$PSI[, samples]) <= 100))
source_files <- config$filtered
stats_cols <- c('FDR', 'invertedFDR', 'betasPsiA', 'betasPsiB', 'deltapsi')
if (args[1] == 'legacy') {
  events <- input$PSI
  ref_stats <- NULL
  alignment <- list()
  for (kind in c('exons', 'introns', 'alt')) {
    path <- paste0('input_tables/legacy/spire_pdiff_', kind, '.csv')
    old <- as.data.frame(fread(path, check.names = FALSE))
    identified <- !is.na(old$EVENT) & nzchar(old$EVENT)
    stopifnot(nrow(old) == nrow(events), all(stats_cols %in% names(old)), any(identified))
    for (col in names(events)) stopifnot(isTRUE(all.equal(old[identified, col],
      events[identified, col], check.attributes = FALSE)))
    if (is.null(ref_stats)) ref_stats <- old[, stats_cols]
    else stopifnot(identical(ref_stats, old[, stats_cols]))
    alignment[[kind]] <- data.frame(family = kind, rows = nrow(old), identified = sum(identified), mismatches = 0)
    source_files <- c(source_files, path)
  }
  events <- cbind(events, ref_stats)
  provenance_label <- 'PROVISIONAL: legacy betAS; Spire1/2 DKO vs Cont'
} else {
  path <- args[1]
  run_dir <- dirname(dirname(path))
  if (!file.exists(file.path(run_dir, 'COMPLETE'))) stop('New simulation must have its COMPLETE marker')
  run <- readRDS(file.path(run_dir, 'run.rds'))
  stopifnot(identical(run$mode, 'full-table'), run$rows == nrow(input$PSI),
    identical(run$filtered_md5, unname(tools::md5sum(config$filtered))),
    identical(run$config$sample_groups, config$sample_groups), run$config$min_reads == 10,
    identical(run$config$group_a, config$group_a), identical(run$config$group_b, config$group_b))
  events <- as.data.frame(fread(path, check.names = FALSE))
  stopifnot(all(c(names(input$PSI), stats_cols) %in% names(events)),
    nrow(events) == nrow(input$PSI), identical(events$EVENT, input$PSI$EVENT))
  for (col in names(input$PSI)) stopifnot(isTRUE(all.equal(events[[col]], input$PSI[[col]], check.attributes = FALSE)))
  source_files <- c(source_files, path, file.path(run_dir, c('COMPLETE', 'run.rds')))
  provenance_label <- 'Spire1/2 DKO vs Cont; new full-table betAS run'
}
events$direction <- factor(classify(events$FDR, events$deltapsi), levels = c('up', 'down', 'ndiff'))
events$event_type <- factor(ifelse(grepl('^MmuINT', events$EVENT), 'IR',
  ifelse(grepl('^MmuALTA', events$EVENT), 'Alt3', ifelse(grepl('^MmuALTD', events$EVENT), 'Alt5', 'EX'))),
  levels = c('EX', 'IR', 'Alt5', 'Alt3'))
stopifnot(!anyNA(events$event_type), all(abs(events$deltapsi - (events$betasPsiB - events$betasPsiA)) < 1e-8))
dir.create(file.path(out, 'figures'), recursive = TRUE)
dir.create(file.path(out, 'tables'))
dir.create(file.path(out, 'external_inputs'))
write_table(events, 'all_events')
write_table(events[events$direction != 'ndiff', ], 'differential_events')
if (exists('alignment')) write_table(bind_rows(alignment), 'legacy_alignment')
figure_names <- character()
main_plots <- list()
source('analysis/splicing.R')
status <- data.frame(step = 'core_splicing', status = 'complete')
annotation_path <- file.path('data/reference', 'spire_event_annotations.tsv')
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
    patchwork::plot_annotation(title = 'Spire1/2 DKO: splicing landscape', tag_levels = 'A')
  save_plot(combined, 'Figure2_Spire_splicing_panels', 14, if (length(names2) > 4) 15 else 10)
  names3 <- intersect(c('Fig3A_protein_impact', 'GO_paper_OrgDb_all_BP'), names(main_plots))
  if (length(names3) == 2) save_plot(patchwork::wrap_plots(main_plots[names3], ncol = 1) +
    patchwork::plot_annotation(title = 'Spire1/2 DKO: predicted functional impact', tag_levels = 'A'),
    'Figure3_Spire_functional_panels', 12, 14)
}
write_table(status, 'analysis_status')
source_files <- unique(c(source_files, list.files('analysis', pattern = '[.](R|py|sh|md|tsv)$', full.names = TRUE),
  'config/analysis.R', '.R-version', 'renv.lock', 'analysis/renv.lock',
  list.files('data/reference', pattern = 'sources.json$', full.names = TRUE)))
reference_sources <- list.files('data/reference', pattern = 'sources.json$', full.names = TRUE)
for (path in reference_sources) file.copy(path, file.path(out, 'tables', basename(path)))
file.copy('analysis/renv.lock', file.path(out, 'analysis_renv.lock'))
jsonlite::write_json(list(source = args[1], min_reads_all_samples = 10L,
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
html <- c('<!doctype html><html lang="en"><meta charset="utf-8"><title>Spire splicing analysis</title>',
  '<style>body{font:16px system-ui;max-width:1100px;margin:30px auto;padding:0 20px}img{width:100%}article{border-top:1px solid #ccc;margin:30px 0}a{color:#514a77}</style>',
  '<h1>Spire1/2 DKO splicing</h1>', paste0('<p><strong>', provenance_label, '</strong></p>'),
  sprintf('<p>N10: %s events; selected: %s events in %s genes. FPR ≤0.05, |ΔPSI| ≥0.10. Positive ΔPSI = greater inclusion/site use/retention in DKO.</p>',
    nrow(events), sum(events$direction != 'ndiff'), length(unique(events$GENE[events$direction != 'ndiff']))),
  '<p>betAS FDR is an empirical per-event FPR, not a BH-adjusted p-value. Historical simulation settings are unverified. See analysis/README.md and tables/analysis_status.csv for methods and pending panels.</p>',
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
report <- c('# Spire1/2 DKO splicing', '', provenance_label, '',
  sprintf('Six profiles (3 Cont, 3 DKO); independent animals, pooling and batches remain unverified. N10 retains %s events. The paper criteria select %s events in %s genes.',
    nrow(events), sum(events$direction != 'ndiff'), length(unique(events$GENE[events$direction != 'ndiff']))), '',
  'Open [the figure gallery](index.html). PDFs and PNGs are in figures/; plotted data and full statistical tables are in tables/.', '',
  'betAS FDR is an empirical per-event FPR, not a BH-adjusted FDR. Small fitted differences at PSI boundaries can receive zero FPR; the paper effect-size filter is applied, but does not calibrate that statistic.', '',
  'The source code, sample mapping, parameters, source hashes and R package versions are recorded. Legacy simulation parameters are unverified; replace the legacy statistics with a completed new full-table export before final biological interpretation.', '')
if (exists('go_results')) report <- c(report,
  sprintf('GO background sensitivity: %s BP terms with the paper OrgDb universe versus %s with the tested-N10 universe. The two analyses are exported separately.',
    bp$enriched_terms[bp$background == 'paper_OrgDb'], bp$enriched_terms[bp$background == 'tested_N10']), '')
if (exists('matched')) report <- c(report,
  sprintf('VASTDB EVENT_INFO lacks %s tested IDs (%s selected). An additional %s IDs have changed coordinates and are excluded from sequence/Matt/RBP inputs. They remain in splicing/GO results. Protein-impact annotations use available mm10 v3 rather than the unavailable v2.3.',
    sum(!matched), sum(!matched & events$direction != 'ndiff'), sum(changed_coords)), '')
report <- c(report, 'Panel status:', '', paste0('- ', status$step, ': ', status$status), '',
  'Drug intersections/convergence require multiple contrasts. Gene expression and experimental imaging/meiosis panels are outside the requested splicing-only analysis. RBP motifs require an external rMAPS2 run; inputs and the checked output parser are prepared. No upload was performed.')
writeLines(report, file.path(out, 'REPORT.md'))
stopifnot(all(file.exists(file.path(out, 'figures', paste0(figs, '.pdf')))),
  sum(table(events$direction)) == nrow(events))
writeLines(provenance_label, file.path(out, 'COMPLETE'))
cat(length(figs), 'figures generated in', out, '\n')
