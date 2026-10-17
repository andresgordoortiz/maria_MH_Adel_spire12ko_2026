.libPaths(c(normalizePath('analysis/renv/library'),.libPaths()))
source('analysis/common.R')
out <- 'results/biological_review'
main_plots <- list();figure_names <- character()
provenance_label <- 'Original sample labels retained; target-locus patterns require genotype/sample-sheet verification'
x <- fread(file.path(out,'tables/target_locus_read_support.csv'))
ex <- x[GENE %in% c('Spire1','Spire2','Fmn2') & COMPLEX %in% c('S','C1','C2','C3','MIC')]
sample_key <- rbind(
 cbind(dataset='spire',fread('metadata/samples.csv')),
 cbind(dataset='fmn2',fread('metadata/fmn2/samples.csv')))
stopifnot(nrow(sample_key)==10L,!anyDuplicated(sample_key$sample),
 setequal(ex$sample,sample_key$sample),!anyDuplicated(ex[,.(dataset,EVENT,sample)]),
 all(is.finite(ex$corrected_total)),all(ex$corrected_total>=0),
 all(abs(ex$corrected_total-ex$corrected_inclusion-ex$corrected_exclusion)<1e-8))
ex <- merge(ex,sample_key[,.(dataset,sample,group,replicate_label)],by=c('dataset','sample'),sort=FALSE)
expected_events <- c(Spire1=13L,Spire2=6L,Fmn2=15L)
event_counts <- ex[,.(events=uniqueN(EVENT)),by=.(GENE,sample)]
stopifnot(nrow(event_counts)==30L,
 all(event_counts$events==expected_events[event_counts$GENE]),
 all(ex[,.(profiles=uniqueN(sample)),by=EVENT]$profiles==10L))
ex$start <- as.integer(sub('.*:([0-9]+)-.*','\\1',ex$COORD))
ex$ordered_start <- ifelse(ex$GENE=='Spire1',-ex$start,ex$start)
setorder(ex,GENE,ordered_start,dataset,sample)
ex[,support_z := {
 spread <- sd(corrected_total)
 if (spread==0) rep(0,.N) else (corrected_total-mean(corrected_total))/spread
},by=.(dataset,EVENT)]
stopifnot(all(is.finite(ex$support_z)),
 all(abs(ex[,.(mean_z=mean(support_z)),by=.(dataset,EVENT)]$mean_z)<1e-10))
summary <- ex[,.(nonANN_exons=.N,exons_N10=sum(corrected_total>=10),
 median_corrected_event_support=median(corrected_total)),by=.(dataset,GENE,sample,group,replicate_label)]
for (ds in c('spire','fmn2')) {
 genes <- if (ds=='spire') c('Spire1','Spire2') else 'Fmn2'
 write_table(ex[GENE %in% genes],paste0(ds,'_exon_support_by_coordinate'))
 write_table(summary[GENE %in% genes],paste0(ds,'_exon_support_summary'))
}
ex$sample_label <- factor(paste(ex$group,ex$replicate_label,sep='_'),
 levels=paste(sample_key$group,sample_key$replicate_label,sep='_'))
ex$dataset <- factor(ex$dataset,levels=c('spire','fmn2'),labels=c('Spire experiment','FMN2 experiment'))
ex$GENE <- factor(ex$GENE,levels=c('Spire1','Spire2','Fmn2'))
ord <- unique(ex[order(GENE,ordered_start),.(GENE,EVENT,COORD)])
ex$event_label <- factor(paste(ex$GENE,ex$EVENT,ex$COORD),levels=rev(paste(ord$GENE,ord$EVENT,ord$COORD)))
for (name in c('Spire','Fmn2')) {
 plotted <- ex[GENE %in% if (name=='Spire') c('Spire1','Spire2') else 'Fmn2']
 save_plot(ggplot(plotted,aes(sample_label,event_label,fill=support_z)) +
  geom_tile(color='white') +
  geom_text(aes(label=sprintf('%g',round(corrected_total,1))),size=3.2) +
  facet_grid(GENE~dataset,scales='free',space='free') +
  scale_fill_gradient2(name='Within-dataset\nexon z-score',low=signed_colors[1],mid=signed_colors[2],high=signed_colors[3],
   midpoint=0,limits=c(-2.5,2.5),breaks=-2:2) +
  labs(x=NULL,y=NULL,title=paste(if (name=='Spire') 'Spire1 and Spire2' else 'Formin2 (Fmn2)',
   'exon support across both datasets'),
   subtitle='Numbers = absolute corrected read support; colour = exon z-score within each dataset; rows 5′ to 3′'),
  paste0(name,'_exon_support_by_coordinate'),16,if (name=='Spire') 10 else 8)
}
summary$sample_label <- factor(paste(summary$group,summary$replicate_label,sep='_'),
 levels=paste(sample_key$group,sample_key$replicate_label,sep='_'))
summary$dataset <- factor(summary$dataset,levels=c('spire','fmn2'),labels=c('Spire experiment','FMN2 experiment'))
p1 <- ggplot(summary,aes(sample_label,GENE,fill=median_corrected_event_support)) + geom_tile(color='white') +
 geom_text(aes(label=sprintf('%g',round(median_corrected_event_support,1)),
  color=median_corrected_event_support>max(median_corrected_event_support)*.55),size=4) +
 scale_color_manual(values=c('FALSE'='black','TRUE'='white'),guide='none') +
 facet_grid(~dataset,scales='free_x',space='free_x') +
 scale_fill_viridis_c(option=sequential_palette,direction=-1,
  name='Median corrected\nread support',limits=c(0,max(summary$median_corrected_event_support))) +
 labs(x=NULL,y=NULL,title='Target-locus exon support across both datasets',
  subtitle='Median per measured non-ANN exon event; raw support counts, not normalized gene expression')
save_plot(p1,'Target_locus_read_support',15,4)
x <- x[dataset=='spire']
samples <- unique(x$sample)
labels <- sub('Oocytes_FG_Spire12_','',samples)
x$sample <- factor(x$sample,levels=samples,labels=labels)
ir <- x[EVENT=='MmuINT0151614']
stopifnot(nrow(ir)==6,identical(as.numeric(ir$PSI),c(4.49,6.93,100,76,100,5.66)))
ir$declared_group <- ifelse(grepl('^Cont',ir$sample),'Cont','DKO')
p2 <- ggplot(ir,aes(sample,PSI,fill=declared_group)) + geom_col(width=.65) +
 geom_text(aes(label=sprintf('%.2f',PSI)),vjust=-.4) + scale_fill_manual(values=condition_colors) +
 coord_cartesian(ylim=c(0,110)) + labs(x=NULL,y='Raw intron PSI (%)',fill='Declared group',
 title='Spire1 MmuINT0151614: the replicates are heterogeneous',
 subtitle='betAS fitted delta: +70.1 pp; raw mean delta: +23.4 pp; Pdiff: 0.618')
save_plot(p2,'Spire1_retained_intron_sample_PSI',10,5)
# Label-blind variable-event PCA is descriptive and does not assign genotypes.
all <- fread('results/spire_final/tables/all_events.csv')
mat <- as.matrix(all[,..samples])
means <- rowMeans(mat); variance <- rowMeans((mat-means)^2)
ix <- head(order(variance,decreasing=TRUE),5000)
z <- sweep(mat[ix,,drop=FALSE],1,means[ix],'-') / sqrt(variance[ix])
stopifnot(all(is.finite(z)),ncol(z)==6)
pc <- prcomp(t(z),center=TRUE,scale.=FALSE)
pct <- 100*pc$sdev^2/sum(pc$sdev^2)
points <- data.frame(sample=labels,group=sub('_.*','',labels),PC1=pc$x[,1],PC2=pc$x[,2])
write_table(points,'spire_label_blind_PCA_coordinates')
write_table(data.frame(EVENT=all$EVENT[ix],PSI_variance=variance[ix]),'spire_label_blind_PCA_events')
save_plot(ggplot(points,aes(PC1,PC2,color=group,label=sample)) + geom_point(size=3) +
 ggrepel::geom_text_repel(seed=42,max.overlaps=Inf) + scale_color_manual(values=condition_colors) +
 labs(title='Spire sample structure without using group labels to select events',
 subtitle='Top 5,000 PSI-variance events, centered and scaled per event; PCA cannot establish sample genotype',
 x=sprintf('PC1 (%.1f%%)',pct[1]),y=sprintf('PC2 (%.1f%%)',pct[2])), 'Spire_label_blind_PCA',10,6)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
cat('Review figures generated; original analysis files and sample labels untouched.\n')
