.libPaths(c(normalizePath('analysis/renv/library'),.libPaths()))
source('analysis/common.R')
source('analysis/rbp_helpers.R')
for (kind in c('EX','IR','Alt5','Alt3')) {
  r <- rbp_regions(kind)
  stopifnot(nrow(r)==if (kind=='IR') 5 else 8,all(diff(r$mid)>0),
    sum(r$width)==switch(kind,EX=1200,IR=450,800))
}
fixture <- expand_grid(dataset=c('spire','fmn2'),direction=c('up','down'),
  motif=c('A.one','A.two','B.one'),region=paste0('R',1:8))
fixture$type <- 'EX';fixture$RBP <- sub('[.].*','',fixture$motif);fixture$pvalue <- .8
fixture$pvalue[fixture$motif=='A.one' & fixture$region=='R1'] <- .01
fixture$pvalue[fixture$motif=='B.one' & fixture$region=='R2'] <- 0
selected <- rbp_panel_data(fixture,'EX')
stopifnot(nrow(selected)==2*2*2*8,!any(selected$motif=='A.two'),
  all(selected$score[selected$pvalue==0]==300),all(is.finite(selected$score)))
collapsed <- rbp_minima(fixture)
stopifnot(nrow(collapsed)==2*2*2*8,
  all(collapsed$pvalue[collapsed$RBP=='A' & collapsed$region=='R1']==.01))
tests <- rbp_box_tests(selected)
stopifnot(nrow(tests)==16,all(tests$n_up==2),all(tests$n_down==2),
  all(is.na(tests$pvalue) | tests$pvalue>=0 & tests$pvalue<=1))
cat('PASS: region geometry, motif eligibility, retained zero p-values, alias minima and boxplot denominators.\n')
