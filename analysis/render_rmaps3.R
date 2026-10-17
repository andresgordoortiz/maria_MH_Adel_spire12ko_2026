.libPaths(c(normalizePath('analysis/renv/library'),.libPaths()))
options(warn=1)
source('analysis/common.R')
source('analysis/rbp_helpers.R')
rows <- list()
for (dataset in c('spire','fmn2')) {
  analysis_dir <- file.path('results',paste0(dataset,'_final'))
  raw <- file.path(analysis_dir,'rmaps3_outputs')
  stopifnot(file.exists(file.path(raw,'COMPLETE')))
  metadata <- jsonlite::read_json(file.path(raw,'provenance.json'))
  stopifnot(metadata$counter=='corrected_half_open_windows',metadata$stat_method=='mannwhitney_greater')
  df <- read_rmaps_motifs(raw,analysis_dir)
  df$dataset <- dataset
  stopifnot(nrow(df)==114*29*2)
  rows[[dataset]] <- df
}
motifs <- bind_rows(rows)
minima <- rbp_minima(motifs)
stopifnot(nrow(minima)==5278*2)
out <- 'results/rbp'
dir.create(file.path(out,'figures'),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(out,'tables'),recursive=TRUE,showWarnings=FALSE)
main_plots <- list();figure_names <- character()
provenance_label <- 'Corrected rMAPS3 counter; raw p-values; sparse-count approximation can be optimistic; original sample labels'
write_table(motifs,'RBP_all_motif_pvalues')
write_table(minima,'RBP_minimum_motif_pvalues')
source('analysis/rbp_paper_figures.R')
stopifnot(all(c('2H','S6A','S6B','S6C','S6D','S6E') %in% registry$paper_panel),nrow(registry)==20)
joined <- inner_join(filter(minima,dataset=='spire') %>% select(-dataset),
  filter(minima,dataset=='fmn2') %>% select(-dataset),
  by=c('RBP','region','type','direction'),suffix=c('_spire','_fmn2'))
stopifnot(nrow(joined)==5278,!anyDuplicated(joined[c('RBP','region','type','direction')]))
joined$raw_p_le_005_both <- joined$pvalue_spire<=.05 & joined$pvalue_fmn2<=.05
write_table(joined,'rMAPS3_shared_RBP_regions')
write_table(filter(joined,raw_p_le_005_both),'rMAPS3_shared_RBP_regions_raw_p005')
write_table(joined %>% group_by(type,direction) %>% summarise(RBPs_tested=n_distinct(RBP),
  RBPs_raw_p005_both=n_distinct(RBP[raw_p_le_005_both]),
  RBP_regions_raw_p005_both=sum(raw_p_le_005_both),.groups='drop'),'rMAPS3_shared_RBP_summary')
html <- c('<!doctype html><html lang="en"><meta charset="utf-8"><title>RBP paper figures</title>',
  '<style>body{font:17px/1.5 system-ui;max-width:1300px;margin:30px auto;padding:0 20px}img{max-width:100%}article{margin:50px 0}a{margin-right:10px}</style>',
  '<h1>Spire and FMN2: RBP figures</h1>',
  '<p>Figure 2H and all splicing-only Figure S6 panels, adapted side by side for Spire and FMN2. Alternative-site panels and count summaries follow below.</p>',
  '<p>Raw motif p-values and descriptive score distributions are exploratory. Sparse counts can give overly small asymptotic p-values; Spire sample identities remain unresolved.</p>',
  '<p><a href="REPORT.md">Panel coverage and methods</a><a href="tables/paper_panel_manifest.csv">Figure manifest</a><a href="../../analysis/RMAPS3.md">Reproduce</a><a href="../index.html">All splicing results</a></p>',
  '<p><a href="../spire_final/rmaps3_outputs/index.html">Spire raw maps</a><a href="../fmn2_final/rmaps3_outputs/index.html">FMN2 raw maps</a><a href="tables/rMAPS3_shared_RBP_summary.csv">Shared raw-p candidates</a></p>')
registry <- registry[order(match(registry$paper_panel,c('2H','S6A','S6B','S6C','S6D','S6E','Repository extension','Repository count summary')),registry$type),]
for (i in seq_len(nrow(registry))) {
  r <- registry[i,]
  html <- c(html,paste0('<article id="',r$figure,'"><h2>',r$paper_panel,' · ',r$type,'</h2><p>',r$figure,
    '</p><a href="figures/',r$figure,'.pdf">PDF</a><a href="figures/',r$figure,'.png">PNG</a><img loading="lazy" src="figures/',r$figure,'.png" alt="',r$figure,'"></article>'))
}
writeLines(c(html,'</html>'),file.path(out,'index.html'))
report <- c('# RBP paper figures: Spire and FMN2','',
  'The active RBP figures are in one [gallery](index.html). Six manuscript panels are present; alternative-site panels and the repository’s count summaries are explicit extensions. Full raw motif maps remain with each dataset.','',
  '| Paper panel | Adapted figure |','|---|---|',
  vapply(which(registry$paper_panel %in% c('2H','S6A','S6B','S6C','S6D','S6E')),function(i)
    paste0('| ',registry$paper_panel[i],' | [',registry$figure[i],'](figures/',registry$figure[i],'.pdf) |'),character(1)),'',
  'Figure 2H uses the reference code’s top 30 exon RBP labels across datasets; S6D uses its top 15 intron labels. Both use Euclidean/complete row clustering, Cividis colours saturating at −log10(p)=8, and columns ordered region → dataset → event direction. Annotation strips identify dataset and event outcome. Colours follow the shared project palette in analysis/palettes.R. The installed ggplot2/patchwork tools reproduce this layout without installing ComplexHeatmap.','',
  'Tracks follow the reference selections: top 15 labels per dataset/direction across regions (ties retained), top three distinct labels annotated per region, and jitter seed 42. The deposited LOESS span 0.3 gives singular local fits at these five/eight distinct region positions; this adaptation uses span 0.75. All 16 dataset/direction fits and prediction standard errors were finite without singularity warnings. The manuscript caption says top 15 per position; the deposited code selects across positions. Region positions and the descriptive curves are schematic distributions, not nucleotide binding coordinates.','',
  'Boxplots retain each motif with raw p≤0.05 in any region within its dataset/direction, as in the reference code, then show all its region scores. Wilcoxon stars compare these motif-score distributions between directions; they are raw descriptive comparisons, not tests of biological replicates or corrected RBP significance. Score tables and test tables accompany every panel. Zero p-values are retained with a plotting floor of 1e−300.','',
  'Shared candidates require raw p≤0.05 for the same RBP label, broad region and direction in both datasets. A label may have different winning motif variants and windows. This criterion is broad: 85 of 91 labels qualify for IR-up and 86 for IR-down, so it does not identify a specific common regulator. See [all paired results](tables/rMAPS3_shared_RBP_regions.csv) and the [summary](tables/rMAPS3_shared_RBP_summary.csv).','',
  'Direct sequence spot checks reproduce the CLI tests but expose sparse-count calibration limits. For FMN2 EX, 1 positive event among 61 versus 2 among 19,992 background events produces approximate p≈1.4×10⁻²⁵; an exact conditional allocation gives p≈0.0091 before accounting for all searched windows/motifs. Even shared SRSF1 intron scores have different winning motifs: the Spire window has 21/127 positive events, while the FMN2 window has 1/15. See [the recorded support counts](tables/rMAPS3_shared_SRSF1_IR_support.csv).','',
  '[Spire output audit](../spire_final/rmaps3_outputs/output_audit.json) · [FMN2 output audit](../fmn2_final/rmaps3_outputs/output_audit.json) · [Spire sparse-count audit](../spire_final/rmaps3_outputs/real_window_audit.json) · [FMN2 sparse-count audit](../fmn2_final/rmaps3_outputs/real_window_audit.json). All 741,000 window p-values per dataset, regional minima, FASTA identities/order and 456 PDF/PNG pairs were checked.','',
  'The counter correction and original failure reproducer are now maintained in the sibling rMAPS3 fork branch fix/motif-window-counting (commit 1076240). Completed scans used the same counting function. The fork distinguishes the demonstrated hit-loss errors from boundary conventions and the separate decision to suppress temporary matrices. See [methods](../../analysis/RMAPS3.md) and the [illustrated explanation](../../audit/rmaps3/explanation/README.md).','',
  'RBP expression/abundance plots remain outside the user’s splicing-only scope; matching expression data are unavailable. The paper’s manually added SF3B1 motif is absent from this catalogue. Spire sample identities remain unresolved; FMN2 controls are heterozygous, with two pools per group. Motif aliases, dependent events, sequence composition, dataset-specific backgrounds, sparse counts and multiple searching limit inference. No common causal RBP is established.','',
  'Previous figures, Fisher baselines and the exact source used for the scans are under legacy/rmaps3/. Preflights and the relocation checksum manifest are under audit/rmaps3/. Historical provenance paths are preserved as written; the relocation manifest maps them to their current files.','')
writeLines(report,file.path(out,'REPORT.md'))
replace_section <- function(page,id,block) {
  pattern <- paste0('<section id="',id,'">[\\s\\S]*?</section>')
  if (grepl(pattern,page,perl=TRUE)) sub(pattern,block,page,perl=TRUE)
  else {stopifnot(grepl('</html>',page,fixed=TRUE));sub('</html>',paste0(block,'\n</html>'),page,fixed=TRUE)}
}
for (dataset in c('spire','fmn2')) {
  path <- file.path('results',paste0(dataset,'_final'))
  raw <- file.path(path,'rmaps3_outputs')
  page <- c('<!doctype html><html lang="en"><meta charset="utf-8"><title>Raw rMAPS3 maps</title>',
    '<style>body{font:17px/1.5 system-ui;max-width:1100px;margin:30px auto;padding:0 20px}</style>',
    paste0('<h1>',dataset,' raw rMAPS3 motif maps</h1>'),
    '<p><a href="../../rbp/index.html">Paper figures and comparison</a> · <a href="../../rbp/REPORT.md">Methods and limitations</a> · <a href="provenance.json">As-run provenance</a> · <a href="event_counts.csv">Event counts</a> · <a href="real_window_audit.json">Sequence spot checks</a></p>')
  for (kind in c('EX','IR','Alt5','Alt3')) {
    pdfs <- sort(list.files(file.path(raw,kind,'maps'),pattern='[.]pdf$'))
    stopifnot(length(pdfs)==114)
    labels <- gsub('>','&gt;',gsub('<','&lt;',gsub('&','&amp;',pdfs,fixed=TRUE),fixed=TRUE),fixed=TRUE)
    links <- paste0(kind,'/maps/',vapply(pdfs,URLencode,character(1),reserved=TRUE))
    page <- c(page,paste0('<h2>',kind,'</h2><p><a href="',kind,'/pVal.up.vs.bg.RNAmap.txt">Up table</a> · <a href="',kind,'/pVal.dn.vs.bg.RNAmap.txt">Down table</a></p><details><summary>114 individual motif maps</summary><ul>'),
      paste0('<li><a href="',links,'">',labels,'</a></li>'),'</ul></details>')
  }
  writeLines(c(page,'</html>'),file.path(raw,'index.html'))
  gallery <- file.path(path,'index.html');page <- readChar(gallery,file.info(gallery)$size)
  page <- replace_section(page,'rmaps3-results','<section id="rmaps3-results"><h2>RBP motif analysis</h2><p><a href="../rbp/index.html">Figure 2H, all splicing Figure S6 panels and comparison</a> · <a href="rmaps3_outputs/index.html">Raw motif maps</a></p></section>')
  writeLines(page,gallery)
  panel_status <- as.data.frame(fread(file.path(path,'tables/analysis_status.csv')))
  panel_status$status[panel_status$step=='RBP_heatmaps_and_region_summaries'] <- 'complete; paper panels 2H and S6A-E in results/rbp'
  fwrite(panel_status,file.path(path,'tables/analysis_status.csv'))
  report_path <- file.path(path,'REPORT.md');text <- readLines(report_path)
  text[startsWith(text,'- RBP_heatmaps_and_region_summaries:')] <- '- RBP_heatmaps_and_region_summaries: complete; see [RBP paper panels](../rbp/index.html)'
  writeLines(text,report_path)
}
gallery <- 'results/comparison/index.html';page <- readChar(gallery,file.info(gallery)$size)
writeLines(replace_section(page,'rmaps3-comparison','<section id="rmaps3-comparison"><h2>RBP comparison</h2><p><a href="../rbp/index.html">Paper figures and side-by-side Spire–FMN2 comparison</a> · <a href="../rbp/REPORT.md">Panel coverage and limitations</a></p></section>'),gallery)
metadata <- list(R=R.version.string,source_md5=as.list(tools::md5sum(c('analysis/render_rmaps3.R',
  'analysis/rbp_helpers.R','analysis/rbp_paper_figures.R','analysis/common.R','analysis/palettes.R','analysis/renv.lock'))),
  counter='corrected_half_open_windows',stat_method='mannwhitney_greater',panels=registry)
jsonlite::write_json(metadata,file.path(out,'tables/render_provenance.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,'tables/sessionInfo.txt'))
cat('Rendered all six splicing-only manuscript RBP panels and 14 repository extensions.\n')
