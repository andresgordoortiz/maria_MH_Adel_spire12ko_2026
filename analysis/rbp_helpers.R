read_rmaps_motifs <- function(rmaps_dir, analysis_dir) {
  stopifnot(identical(read.csv(file.path(rmaps_dir,'input_checksums.csv')),
    read.csv(file.path(analysis_dir,'external_inputs','input_checksums.csv'))))
  rows <- list()
  for (kind in c('EX','IR','Alt5','Alt3')) for (direction in c('up','down')) {
    name <- if (direction=='up') 'pVal.up.vs.bg.RNAmap.txt' else 'pVal.dn.vs.bg.RNAmap.txt'
    path <- file.path(rmaps_dir,kind,name)
    if (!file.exists(path)) next
    df <- as.data.frame(data.table::fread(path))
    regions <- if (kind=='EX') 8L else if (kind=='IR') c(4L,5L) else c(6L,8L)
    stopifnot((ncol(df)-1L) %in% regions)
    names(df) <- c('motif',paste0('R',seq_len(ncol(df)-1L)))
    for (column in names(df)[-1])
      if (is.logical(df[[column]]) && all(is.na(df[[column]]))) df[[column]] <- as.numeric(df[[column]])
    stopifnot(!anyNA(df$motif),all(nzchar(df$motif)),!anyDuplicated(df$motif),
      all(vapply(df[-1],is.numeric,logical(1))))
    df <- tidyr::pivot_longer(df,-motif,names_to='region',values_to='pvalue')
    if (any(!is.na(df$pvalue) & (!is.finite(df$pvalue) | df$pvalue<0 | df$pvalue>1)))
      stop('Invalid rMAPS p-value')
    df$RBP <- sub('[.].*','',df$motif)
    df$type <- kind; df$direction <- direction
    rows[[length(rows)+1L]] <- dplyr::filter(df,!is.na(pvalue))
  }
  if (!length(rows)) stop('No supported rMAPS output tables')
  dplyr::bind_rows(rows)
}

rbp_minima <- function(motifs) {
  motifs %>% group_by(RBP,region,type,direction,dataset) %>%
    summarise(pvalue=min(pvalue),.groups='drop')
}

rbp_regions <- function(kind) {
  labels <- switch(kind,
    EX=c('Upstream exon\n3′ end','Upstream intron\n5′ end','Upstream intron\n3′ end',
      'Target exon\n5′ end','Target exon\n3′ end','Downstream intron\n5′ end',
      'Downstream intron\n3′ end','Downstream exon\n5′ end'),
    IR=c('Upstream exon\n5′ end','Upstream exon\n3′ end','Intron body\n(scaled)',
      'Downstream exon\n5′ end','Downstream exon\n3′ end'),
    Alt5=c('Common exon\n5′ end','Common exon\n3′ end','Alternative segment\n5′ end',
      'Alternative segment\n3′ end','Intron after\nlong exon','Intron before\nflank',
      'Flanking exon\n5′ end','Flanking exon\n3′ end'),
    Alt3=c('Flanking exon\n5′ end','Flanking exon\n3′ end','Following intron\n5′ end',
      'Following intron\n3′ end','Alternative segment\n5′ end','Alternative segment\n3′ end',
      'Common exon\n5′ end','Common exon\n3′ end'))
  widths <- switch(kind,EX=c(50,250,250,50,50,250,250,50),IR=c(50,50,250,50,50),
    Alt5=c(50,50,50,50,250,250,50,50),Alt3=c(50,50,250,250,50,50,50,50))
  data.frame(region=paste0('R',seq_along(labels)),label=labels,width=widths,
    end=cumsum(widths),mid=cumsum(widths)-widths/2)
}

rbp_panel_data <- function(motifs,kind) {
  regions <- rbp_regions(kind)
  df <- filter(motifs,type==kind) %>% group_by(dataset,direction,motif) %>%
    filter(any(pvalue<=.05)) %>% ungroup() %>% left_join(regions,by='region')
  stopifnot(!anyNA(df$mid))
  df$score <- -log10(pmax(df$pvalue,1e-300))
  df$region <- factor(df$region,levels=regions$region)
  df$dataset <- factor(df$dataset,levels=c('spire','fmn2'))
  df$direction <- factor(df$direction,levels=c('up','down'))
  df$outcome <- factor(as.character(df$direction),levels=c('up','down'),
    labels=switch(kind,EX=c('Included','Skipped'),IR=c('Retained','Excised'),
      c('Higher site usage','Lower site usage')))
  df
}

rbp_box_tests <- function(df) {
  df %>% group_by(dataset,region) %>% summarise(
    n_up=sum(direction=='up'),n_down=sum(direction=='down'),
    pvalue=if (all(c('up','down') %in% direction))
      suppressWarnings(wilcox.test(score[direction=='up'],score[direction=='down'],exact=FALSE)$p.value)
      else NA_real_,
    y=max(score)*1.1+.5,.groups='drop') %>% mutate(
      label=case_when(is.na(pvalue)~'NA',pvalue<.0001~'****',pvalue<.001~'***',
        pvalue<.01~'**',pvalue<.05~'*',TRUE~'ns'))
}
