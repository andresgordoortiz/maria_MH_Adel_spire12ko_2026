# Fixed hypothesis families and background/selection diagnostics for the publication.
go_themes <- c('GO:0051301', 'GO:0051321', 'GO:0007059', 'GO:0007051', 'GO:0030036',
  'GO:0005635', 'GO:0006338', 'GO:0016607', 'GO:0008380')
engineered <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db, keys = c('Spire1', 'Spire2'),
  keytype = 'SYMBOL', columns = 'ENTREZID'))$ENTREZID
all_symbols <- valid_genes(c(unlist(lapply(events, function(x) x$GENE)), dge$gene_name[is.finite(dge$padj)]))
direct <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db,
  keys = intersect(all_symbols, AnnotationDbi::keys(org.Mm.eg.db, 'SYMBOL')), keytype = 'SYMBOL', columns = 'ENTREZID'))
direct <- unique(direct[!is.na(direct$ENTREZID), ])
missing <- setdiff(all_symbols, direct$SYMBOL)
alias <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db,
  keys = intersect(missing, AnnotationDbi::keys(org.Mm.eg.db, 'ALIAS')), keytype = 'ALIAS', columns = 'ENTREZID'))
alias <- unique(alias[!is.na(alias$ENTREZID), ]); names(alias)[1] <- 'SYMBOL'
alias$method <- 'unique alias'; direct$method <- 'direct symbol'
mapping_all <- bind_rows(direct, alias) %>% group_by(SYMBOL) %>%
  filter(n_distinct(ENTREZID) == 1) %>% ungroup() %>% distinct(SYMBOL, ENTREZID, .keep_all = TRUE)
export(mapping_all, 'S6_GO_resolved_symbol_mapping')
export(alias, 'S6_GO_alias_candidates')
export(left_join(data.frame(SYMBOL = all_symbols), mapping_all, by = 'SYMBOL'), 'S6_GO_all_symbol_mapping_audit')
membership <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db, keys = names(division_terms),
  keytype = 'GOALL', columns = c('ENTREZID', 'EVIDENCEALL')))
division <- merge(membership, mapping_all, by = 'ENTREZID')
division$term <- unname(division_terms[division$GOALL])
division_selected <- robust[robust$GENE %in% division$SYMBOL & !robust$GENE %in% c('Spire1', 'Spire2'), ]
export(division, 'S4_division_GO_membership')
export(division_selected, 'S4_division_selected_splicing')
go_revised <- go_matched <- go_mapping <- go_opportunities <- fpr_diagnostics <- list()
for (ds in names(studies)) {
  mapping <- mapping_all[mapping_all$SYMBOL %in% events[[ds]]$GENE, ]
  map <- function(g) setdiff(unique(mapping$ENTREZID[mapping$SYMBOL %in% g]), engineered)
  s <- sensitivity[[ds]]
  rules <- list(selected = s$GENE, consistent = s$GENE[s$loo_min_signed_pp >= 10],
    separated = s$GENE[s$complete_separation])
  for (rule in names(rules)) {
    symbols <- setdiff(valid_genes(rules[[rule]]), c('Spire1', 'Spire2'))
    go_mapping[[paste(ds, rule)]] <- data.frame(dataset = ds, selection = rule,
      selected_symbols = length(symbols), mapped_symbols = sum(symbols %in% mapping$SYMBOL),
      rescued_aliases = sum(symbols %in% mapping$SYMBOL[mapping$method == 'unique alias']),
      unique_Entrez = length(map(symbols)), unmapped = paste(setdiff(symbols, mapping$SYMBOL), collapse = ';'))
  }
  opportunity_genes <- events[[ds]] %>% filter(GENE %in% valid_genes(GENE), !GENE %in% c('Spire1', 'Spire2')) %>%
    group_by(GENE) %>% summarise(tested_events = n(), selected = any(direction != 'ndiff'), .groups = 'drop')
  opportunity_genes$stratum <- floor(log2(opportunity_genes$tested_events)); opportunity_genes$dataset <- ds
  go_opportunities[[ds]] <- opportunity_genes
  degrees <- merge(unique(events[[ds]][, c('GENE', 'EVENT')]), mapping, by.x = 'GENE', by.y = 'SYMBOL') %>%
    group_by(ENTREZID) %>% summarise(tested_events = n_distinct(EVENT), .groups = 'drop')
  for (ont in c('BP', 'CC', 'MF')) {
    objects <- lapply(c('tested_N10', 'paper_OrgDb'), function(bg) {
      path <- file.path('reports', studies[[ds]], 'tables', paste0('GO_', bg, '_all_', ont, '.rds'))
      input_files <<- union(input_files, path)
      object <- readRDS(path)
      old <- object@result
      check <- go_counts(object@universe, intersect(object@gene, object@universe), object@geneSets[old$ID])
      stopifnot(max(abs(check$pvalue[match(old$ID, check$ID)] - old$pvalue)) < 1e-12)
      object
    })
    names(objects) <- c('tested_N10', 'paper_OrgDb')
    primary <- objects$paper_OrgDb
    universe <- intersect(map(events[[ds]]$GENE), primary@universe)
    sizes <- lengths(lapply(primary@geneSets, intersect, universe))
    eligible <- names(sizes)[sizes >= 10 & sizes <= 500]
    terms <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = eligible, keytype = 'GOID', columns = 'TERM'))
    descriptions <- setNames(terms$TERM, terms$GOID)
    stopifnot(!anyNA(descriptions[eligible]))
    for (bg in names(objects)) for (rule in names(rules)) {
      object <- primary
      u <- if (bg == 'tested_N10') universe else setdiff(object@universe, engineered)
      query <- intersect(map(rules[[rule]]), u)
      stopifnot(all(eligible %in% names(object@geneSets)))
      x <- go_counts(u, query, object@geneSets[eligible])
      x$Description <- unname(descriptions[x$ID]); x$dataset <- ds; x$ontology <- ont
      x$background <- bg; x$selection <- rule
      x$geneID <- vapply(object@geneSets[eligible], function(g)
        paste(sort(intersect(g, query)), collapse = '/'), character(1))
      x$odds_ratio <- x$CI_low <- x$CI_high <- NA_real_
      for (i in which(x$ID %in% go_themes)) {
        k <- x$Count[i]; n <- x$background_term_genes[i]; q <- x$query_genes[i]; N <- x$universe_genes[i]
        test <- fisher.test(matrix(c(k, q-k, n-k, N-q-n+k), 2), alternative = 'two.sided')
        x$odds_ratio[i] <- unname(test$estimate); x$CI_low[i] <- test$conf.int[1]; x$CI_high[i] <- test$conf.int[2]
      }
      go_revised[[paste(ds, ont, bg, rule)]] <- x
    }
    ids <- intersect(go_themes, eligible)
    if (ont == 'BP') {
      division_union <- Reduce(base::union, primary@geneSets[names(division_terms)])
      selected_u <- universe %in% map(rules$selected)
      degree <- degrees$tested_events[match(universe, degrees$ENTREZID)]
      strata <- split(seq_along(universe), floor(log2(degree)))
      export(data.frame(dataset = ds, selected = sum(selected_u), universe = length(universe),
        observed = sum(selected_u & universe %in% division_union), background_annotated = sum(universe %in% division_union),
        expected_uniform = sum(selected_u)*mean(universe %in% division_union),
        expected_event_matched = sum(vapply(strata, function(ix)
          sum(selected_u[ix])*mean(universe[ix] %in% division_union), numeric(1)))),
        paste0('S6_', ds, '_division_union_counts'))
    }
    if (length(ids)) {
      membership <- vapply(primary@geneSets[ids], function(g) universe %in% g, logical(length(universe)))
      matched <- matched_category_counts(universe %in% map(rules$selected), membership,
        degrees$tested_events[match(universe, degrees$ENTREZID)])
      matched$summary$dataset <- ds; matched$summary$ontology <- ont
      matched$summary$universe_genes <- length(universe)
      go_matched[[paste(ds, ont)]] <- matched$summary
      export(as.data.frame(matched$null), paste0('S6_', ds, '_', ont, '_event_count_matched_null'))
    }
  }
  m <- as.matrix(events[[ds]][, names(groups[[ds]])])
  constant <- apply(m, 1, min) == apply(m, 1, max)
  fpr_diagnostics[[ds]] <- data.frame(dataset = ds, tested = nrow(m),
    constant_raw_PSI = sum(constant), constant_with_zero_FPR = sum(constant & events[[ds]]$FDR == 0),
    constant_selected = sum(constant & events[[ds]]$direction != 'ndiff'),
    zero_FPR = sum(events[[ds]]$FDR == 0), nsim = 1000,
    balanced_label_assignments = choose(ncol(m), ncol(m)/2),
    smallest_exact_two_sided_label_P = 2/choose(ncol(m), ncol(m)/2))
}
go_revised <- bind_rows(go_revised) %>% group_by(dataset, background, selection) %>%
  mutate(q_BH = p.adjust(pvalue, 'BH'), family_tests = n()) %>% ungroup()
go_matched <- bind_rows(go_matched)
go_matched$q_BH <- p.adjust(go_matched$pvalue, 'BH')
go_mapping <- bind_rows(go_mapping); go_opportunities <- bind_rows(go_opportunities)
export(go_revised, 'S6_GO_revised_all_tests')
export(go_matched, 'S6_GO_event_count_matched_themes')
export(go_mapping, 'S6_GO_mapping_audit')
export(go_opportunities, 'S6_GO_gene_event_opportunities')
export(bind_rows(fpr_diagnostics), 'S3_FPR_calibration_diagnostics')
expression_go <- bind_rows(lapply(c('BP', 'CC', 'MF'), function(ont) {
  object <- readRDS(file.path('reports/spire_final/tables', paste0('GO_paper_OrgDb_all_', ont, '.rds')))
  map_expression <- function(g) setdiff(unique(mapping_all$ENTREZID[mapping_all$SYMBOL %in% g]), engineered)
  u <- intersect(map_expression(dge$gene_name[is.finite(dge$padj)]), object@universe)
  q <- intersect(map_expression(dge$gene_name[dge$selected]), u)
  sizes <- lengths(lapply(object@geneSets, intersect, u))
  ids <- names(sizes)[sizes >= 10 & sizes <= 500]
  x <- go_counts(u, q, object@geneSets[ids])
  terms <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = ids, keytype = 'GOID', columns = 'TERM'))
  x$Description <- terms$TERM[match(x$ID, terms$GOID)]; x$ontology <- ont
  x
}))
expression_go$q_BH <- p.adjust(expression_go$pvalue, 'BH')
export(expression_go, 'S6_GO_expression_all_tests')

presence_root <- 'reports/regional_motif_presence'
stopifnot(file.exists(file.path(presence_root, 'COMPLETE')))
checks <- read_input(file.path(presence_root, 'input_checksums.csv'))
for (i in seq_len(nrow(checks))) stopifnot(identical(digest::digest(file = checks$file[i], algo = 'sha256'), checks$sha256[i]))
input_files <- union(input_files, checks$file)
presence <- read_input(file.path(presence_root, 'regional_presence.csv'))
stopifnot(all(presence$positive >= 0 & presence$positive <= presence$n), !anyDuplicated(presence[, c('dataset', 'event_type', 'group', 'region', 'RBP')]))
for (ds in names(studies)) for (kind in c('EX', 'IR', 'Alt5', 'Alt3')) for (group in c('up', 'dn', 'bg')) {
  prefix <- if (ds == 'spire') 'spire_fastq' else 'fmn2'
  path <- file.path('reports', studies[[ds]], 'rmaps3_outputs/inputs', paste0(prefix, '_', kind, '_',
    c(up = 'up', dn = 'down', bg = 'ndiff')[[group]], '.txt'))
  input <- read_input(path)
  n <- unique(presence$n[presence$dataset == ds & presence$event_type == kind & presence$group == group])
  stopifnot(identical(n, nrow(input)))
}
background <- presence[presence$group == 'bg', c('dataset', 'event_type', 'region', 'RBP', 'n', 'positive')]
names(background)[5:6] <- c('background_n', 'background_positive')
motif_tests <- merge(presence[presence$group != 'bg', ], background, by = c('dataset', 'event_type', 'region', 'RBP'))
motif_tests$direction <- ifelse(motif_tests$group == 'up', 'up', 'down')
motif_tests$pvalue <- with(motif_tests, phyper(positive-1, positive+background_positive,
  n+background_n-positive-background_positive, n, lower.tail = FALSE))
motif_tests <- motif_tests %>% group_by(dataset) %>% mutate(q_BH = p.adjust(pvalue, 'BH'), family_tests = n()) %>% ungroup()
motif_tests$log2_OR <- with(motif_tests, log2((positive+.5)*(background_n-background_positive+.5)/
  ((n-positive+.5)*(background_positive+.5))))
stopifnot(nrow(motif_tests) == 2*91*(8+5+8+8)*2)
export(motif_tests, 'S8_RBP_regional_presence_exact_BH')
export(presence, 'S8_RBP_regional_presence_counts')

overlap_inference <- tab('comparison', 'overlap_enrichment')[, c('unit', 'p.value')]
overlap_inference <- rbind(overlap_inference,
  data.frame(unit = 'shared_event_direction', p.value = binom.test(sum(shared$spire_deltapsi*shared$fmn2_deltapsi > 0), nrow(shared))$p.value),
  data.frame(unit = 'event_count_matched_gene_overlap', p.value = degree_overlap$pvalue))
overlap_inference$q_BH <- p.adjust(overlap_inference$p.value, 'BH')
export(overlap_inference, 'S5_overlap_inference_BH')
ledger <- data.frame(family = c('SPIRE expression', 'SPIRE expression GO', 'SPIRE GO primary', 'FMN2 GO primary',
  'GO opportunity sensitivity', 'SPIRE Matt features', 'FMN2 Matt features',
  'SPIRE regional RBP presence', 'FMN2 regional RBP presence', 'Cross-study overlap', 'Splicing selection'),
  tests = c(sum(is.finite(dge$padj)), nrow(expression_go), sum(go_revised$dataset == 'spire' & go_revised$background == 'tested_N10' & go_revised$selection == 'selected'),
    sum(go_revised$dataset == 'fmn2' & go_revised$background == 'tested_N10' & go_revised$selection == 'selected'), nrow(go_matched),
    sum(matt$dataset == 'spire' & is.finite(matt$p.value)), sum(matt$dataset == 'fmn2' & is.finite(matt$p.value)),
    sum(motif_tests$dataset == 'spire'), sum(motif_tests$dataset == 'fmn2'), nrow(overlap_inference), NA),
  correction = c(rep('BH within the complete stated family', 10), 'Effect/FPR selection only; no calibrated FDR claim'),
  notes = c('Saved DESeq2 independent filtering retained', 'BP+CC+MF pooled; finite-padj expression universe', 'BP+CC+MF pooled; engineered loci excluded',
    'BP+CC+MF pooled; engineered loci excluded', 'All 18 question-driven tests; event-count matched',
    'All exported features and three contrasts', 'All exported features and three contrasts',
    'All 91 RBPs, four event types, 29 regions and both directions',
    'All 91 RBPs, four event types, 29 regions and both directions',
    'Seven original overlaps plus direction and matched-gene tests', 'Simulation-tail calibration and biological replication remain limitations'))
export(ledger, 'S1_inference_families')
