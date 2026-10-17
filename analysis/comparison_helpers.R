overlap_test <- function(a, b, universe) {
  a <- unique(a); b <- unique(b); universe <- unique(universe)
  stopifnot(length(universe) > 0, all(a %in% universe), all(b %in% universe))
  shared <- length(intersect(a, b))
  tab <- matrix(c(shared, length(a) - shared, length(b) - shared,
    length(universe) - length(union(a, b))), nrow = 2)
  test <- fisher.test(tab, alternative = 'greater')
  data.frame(universe = length(universe), spire = length(a), fmn2 = length(b), shared = shared,
    expected = length(a) * length(b) / length(universe), odds_ratio = unname(test$estimate), p.value = test$p.value)
}
join_common_events <- function(a, b) {
  stopifnot(!anyDuplicated(a$EVENT), !anyDuplicated(b$EVENT))
  ids <- intersect(a$EVENT, b$EVENT)
  stopifnot(length(ids) > 0)
  a <- a[match(ids, a$EVENT), ]; b <- b[match(ids, b$EVENT), ]
  for (field in c('GENE', 'EVENT', 'COORD', 'COMPLEX')) stopifnot(isTRUE(all.equal(a[[field]], b[[field]], check.attributes = FALSE)))
  ans <- a[, c('GENE', 'EVENT', 'COORD', 'COMPLEX', 'event_type')]
  for (name in c('FDR', 'betasPsiA', 'betasPsiB', 'deltapsi')) {
    ans[[paste0('spire_', name)]] <- a[[name]]
    ans[[paste0('fmn2_', name)]] <- b[[name]]
  }
  ans$spire_selected <- classify(ans$spire_FDR, ans$spire_deltapsi) != 'ndiff'
  ans$fmn2_selected <- classify(ans$fmn2_FDR, ans$fmn2_deltapsi) != 'ndiff'
  both <- ans$spire_selected & ans$fmn2_selected
  ans$comparison <- ifelse(both, ifelse(ans$spire_deltapsi * ans$fmn2_deltapsi > 0,
    'Both: same direction', 'Both: opposite direction'), ifelse(ans$spire_selected,
      'Spire only', ifelse(ans$fmn2_selected, 'FMN2 only', 'Neither')))
  ans$locus <- sub('-[0-9]+/[0-9]+$', '', ans$EVENT)
  ans
}
