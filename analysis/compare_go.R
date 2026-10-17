suppressPackageStartupMessages(library(org.Mm.eg.db))
mapping <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db, keys = gene_universe,
  keytype = 'SYMBOL', columns = 'ENTREZID'))
write_table(mapping, 'common_universe_Entrez_mapping')
entrez <- function(g) unique(na.omit(mapping$ENTREZID[mapping$SYMBOL %in% g]))
common_entrez <- entrez(gene_universe)
sets <- c(gene_sets, list(shared_genes = shared_genes, concordant_event_genes = valid_genes(common$GENE[same])))
go_summary <- go_all <- go_selected <- list()
for (ont in c('BP', 'CC', 'MF')) for (set in names(sets)) {
  message('Comparison GO ', ont, ' ', set)
  ids <- entrez(sets[[set]])
  ego <- if (length(ids) >= 5) suppressMessages(clusterProfiler::enrichGO(gene = ids,
    universe = common_entrez, OrgDb = org.Mm.eg.db, keyType = 'ENTREZID', ont = ont,
    pAdjustMethod = 'BH', pvalueCutoff = .05, qvalueCutoff = .2,
    minGSSize = 10, maxGSSize = 500, readable = TRUE)) else NULL
  empty <- data.frame(ID = character(), Description = character(), p.adjust = numeric(), geneID = character(), Count = integer())
  selected_go <- if (is.null(ego)) empty else as.data.frame(ego)
  all_go <- if (is.null(ego)) empty else ego@result
  key <- paste(set, ont, sep = '_')
  go_all[[key]] <- all_go; go_selected[[key]] <- selected_go
  write_table(selected_go, paste0('GO_common_universe_', key))
  write_table(all_go, paste0('GO_common_universe_', key, '_all_tested_terms'))
  saveRDS(ego, file.path(out, 'tables', paste0('GO_common_universe_', key, '.rds')))
  go_summary[[key]] <- data.frame(set = set, ontology = ont, genes = length(sets[[set]]), mapped = length(ids),
    universe_Entrez = length(common_entrez), tested = length(ids) >= 5, enriched_terms = nrow(selected_go))
  if (set %in% c('shared_genes', 'concordant_event_genes')) {
    if (nrow(selected_go)) {
      x <- head(selected_go[order(selected_go$p.adjust, -selected_go$Count), ], 15)
      x$Description <- factor(x$Description, levels = rev(x$Description))
      p <- ggplot(x, aes(-log10(p.adjust), Description)) + geom_point(aes(size = Count), color = '#222222') +
        labs(x = '-log10 BH-adjusted p', y = NULL, title = paste('Shared-target GO:', set, ont),
          subtitle = paste(length(common_entrez), 'Entrez genes in the common tested universe'))
    } else p <- empty_plot(paste(if (length(ids) < 5) 'Insufficient mapped genes:' else 'No enriched terms:', set, ont))
    save_plot(p, paste0('Comparison_GO_', key), 10, 6)
  }
}
write_table(bind_rows(go_summary), 'GO_common_universe_summary')
shared_terms <- list(); plot_terms <- list()
for (ont in c('BP', 'CC', 'MF')) {
  a <- go_selected[[paste0('spire_', ont)]]; b <- go_selected[[paste0('fmn2_', ont)]]
  for (id in intersect(a$ID, b$ID)) {
    x <- a[a$ID == id, ]; y <- b[b$ID == id, ]
    gx <- strsplit(x$geneID, '/', fixed = TRUE)[[1]]; gy <- strsplit(y$geneID, '/', fixed = TRUE)[[1]]
    shared_terms[[length(shared_terms) + 1L]] <- data.frame(ontology = ont, ID = id, Description = x$Description,
      spire_p_adjust = x$p.adjust, fmn2_p_adjust = y$p.adjust, spire_Count = x$Count, fmn2_Count = y$Count,
      shared_driver_genes = paste(intersect(gx, gy), collapse = ';'), n_shared_drivers = length(intersect(gx, gy)),
      spire_drivers = paste(gx, collapse = ';'), fmn2_drivers = paste(gy, collapse = ';'))
  }
  x <- go_all[[paste0('spire_', ont)]]; y <- go_all[[paste0('fmn2_', ont)]]
  ids <- unique(c(head(x$ID[order(x$p.adjust)], 8), head(y$ID[order(y$p.adjust)], 8)))
  for (ds in c('spire', 'fmn2')) {
    z <- go_all[[paste(ds, ont, sep = '_')]]
    z <- z[z$ID %in% ids, ]
    if (nrow(z)) plot_terms[[length(plot_terms) + 1L]] <- data.frame(dataset = ds, ontology = ont,
      ID = z$ID, Description = z$Description, p.adjust = z$p.adjust, Count = z$Count)
  }
}
shared_terms <- if (length(shared_terms)) bind_rows(shared_terms) else data.frame(ontology = character(), ID = character(), Description = character(),
  spire_p_adjust = numeric(), fmn2_p_adjust = numeric(), shared_driver_genes = character(), n_shared_drivers = integer())
write_table(shared_terms, 'shared_significant_GO_terms_and_drivers')
if (length(plot_terms)) {
  x <- bind_rows(plot_terms)
  x$score <- -log10(pmax(x$p.adjust, 1e-300))
  x$selected <- x$p.adjust <= .05
  write_table(x, 'GO_comparison_plotted_terms')
  save_plot(ggplot(x, aes(dataset, Description)) + geom_point(aes(size = Count, color = score, alpha = selected)) +
    facet_grid(ontology ~ ., scales = 'free_y', space = 'free_y') +
    scale_color_viridis_c(option = sequential_palette, direction = -1) +
    scale_alpha_manual(values = c('FALSE' = .25, 'TRUE' = 1)) +
    labs(x = NULL, y = NULL, color = '-log10 BH p', alpha = 'BH p ≤0.05',
      title = 'GO term profiles on the same tested-gene background', subtitle = 'Top 8 terms in either study per ontology; faded points are not significant'),
    'Comparison_GO_term_profiles', 12, max(8, length(unique(x$ID)) * .22))
}
