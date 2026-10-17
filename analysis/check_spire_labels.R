.libPaths(c(normalizePath('analysis/renv/library'),.libPaths()))
source('analysis/common.R')
out <- 'results/spire_identity_check'
main_plots <- list(); figure_names <- character()
provenance_label <- 'Diagnostic tests only; original labels retained; RNA phenotypes do not establish DNA genotype'
select_markers <- function(m, gene, id, A, B, cutoff) {
 a <- rowMeans(m[,A,drop=FALSE]); b <- rowMeans(m[,B,drop=FALSE])
 gap <- pmax(apply(m[,B,drop=FALSE],1,min)-apply(m[,A,drop=FALSE],1,max),
   apply(m[,A,drop=FALSE],1,min)-apply(m[,B,drop=FALSE],1,max))
 keep <- which(!is.na(gene) & !gene %in% c('','NA','NaN','Spire1','Spire2') & gap>0 & abs(b-a)>=cutoff)
 z <- data.frame(row=keep,GENE=gene[keep],EVENT=id[keep],control_reference=a[keep],KO_reference=b[keep],delta=b[keep]-a[keep],anchor_gap=gap[keep])
 z <- z[order(-z$anchor_gap,z$EVENT),]
 z[!duplicated(z$GENE),]
}
project <- function(m,markers) {
 stopifnot(nrow(markers)>0,all(abs(markers$delta)>0))
 sweep(sweep(m[markers$row,,drop=FALSE],1,markers$control_reference,'-'),1,markers$delta,'/')
}
# Synthetic swap, duplicate-gene filtering and independence from held-out PSI.
fixture <- rbind(c(20,22,82,80,82,21),c(80,78,18,20,18,79),c(50,50,50,50,50,50),c(0,0,100,100,100,0),c(20,21,65,60,61,20))
gene <- c('G1','G2','G3','Spire1','G1');id <- paste0('E',1:5)
f <- select_markers(fixture,gene,id,1:2,4:5,10)
stopifnot(identical(f$EVENT,c('E1','E2')),all(project(fixture,f)[,3]>.9),all(project(fixture,f)[,6]<.1))
changed <- fixture;changed[,c(3,6)] <- -999
stopifnot(identical(f,select_markers(changed,gene,id,1:2,4:5,10)))

x <- as.data.frame(fread('results/spire_final/tables/all_events.csv'))
samples <- names(x)[grepl('^Oocytes_FG_Spire12_',names(x))]
labels <- sub('Oocytes_FG_Spire12_','',samples)
stopifnot(identical(labels,c('Cont_a','Cont_b','Cont_c','DKO_a','DKO_b','DKO_c')))
m <- as.matrix(x[,samples]);stopifnot(all(is.finite(m)),!anyDuplicated(x$EVENT))
primary <- select_markers(m,x$GENE,x$EVENT,1:2,4:5,10)
write_table(primary,'anchor_training_markers')
score <- project(m,primary)
primary_summary <- data.frame(sample=labels,marker_genes=nrow(primary),median_anchor_score=apply(score,2,median),
 fraction_closer_to_KO=colMeans(score>.5),used_in_training=seq_along(samples) %in% c(1,2,4,5))
write_table(primary_summary,'anchor_prediction_summary')
write_table(data.frame(GENE=primary$GENE,EVENT=primary$EVENT,Cont_c=score[,3],DKO_c=score[,6]),'held_out_c_gene_scores')
sens <- list();plot_scores <- list()
for (cutoff in c(10,20)) for (kind in c('all','EX','IR','Alt5','Alt3','exclude_chr8_chr18')) {
 rows <- if (kind=='all') seq_len(nrow(x)) else if (kind=='exclude_chr8_chr18') which(!grepl('^chr(8|18):',x$COORD)) else which(x$event_type==kind)
 z <- select_markers(m[rows,,drop=FALSE],x$GENE[rows],x$EVENT[rows],1:2,4:5,cutoff)
 if (!nrow(z)) next
 p <- project(m[rows,,drop=FALSE],z)
 sens[[length(sens)+1L]] <- data.frame(min_anchor_delta_pp=cutoff,event_subset=kind,marker_genes=nrow(z),
 Cont_c_score=median(p[,3]),DKO_c_score=median(p[,6]),Cont_c_KO_vote=mean(p[,3]>.5),DKO_c_KO_vote=mean(p[,6]>.5))
 if(kind=='all')for (j in c(3,6))plot_scores[[length(plot_scores)+1L]]<-data.frame(sample=labels[j],score=p[,j],cutoff=paste0('Anchor difference ≥',cutoff,' pp'))
}
write_table(bind_rows(sens),'anchor_marker_sensitivity')
ps <- bind_rows(plot_scores)
save_plot(ggplot(ps,aes(sample,score,fill=sample))+geom_violin(scale='width',trim=TRUE)+
 stat_summary(fun=median,geom='point',size=2)+geom_hline(yintercept=c(0,1),color='grey50')+geom_hline(yintercept=.5,linetype='dashed',color='grey50')+
 facet_wrap(~cutoff)+coord_cartesian(ylim=c(-1,2))+
 scale_fill_manual(values=c(Cont_c=condition_colors[['Cont']],DKO_c=condition_colors[['DKO']]))+
 labs(x=NULL,y='Projection: control reference=0; KO reference=1',title='Other genes do not confidently identify both c samples',
 subtitle='Markers selected using a/b samples only; one event per gene; Spire1/2 excluded. This is not a genotype probability.'),
 'Held_out_c_marker_scores',12,6)

cv <- list()
for (A in 1:2) {
 B <- A+3L;test <- setdiff(c(1,2,4,5),c(A,B))
 z <- select_markers(m,x$GENE,x$EVENT,A,B,10);p <- project(m,z)
 for (j in test)cv[[length(cv)+1L]]<-data.frame(training_pair=letters[A],held_out_sample=labels[j],marker_genes=nrow(z),
 median_score=median(p[,j]),declared_group=if(j<4)'Cont' else 'DKO',predicted_group=if(median(p[,j])>.5)'DKO' else 'Cont')
}
cv <- bind_rows(cv);cv$agrees_with_declared_group<-cv$predicted_group==cv$declared_group
write_table(cv,'anchor_cross_validation')
save_plot(ggplot(cv,aes(held_out_sample,median_score,color=declared_group))+geom_point(size=3)+
 geom_hline(yintercept=.5,linetype='dashed')+facet_wrap(~training_pair,scales='free_x')+
 scale_color_manual(values=condition_colors)+
 labs(x=NULL,y='Held-out median score',color='Declared group',title='The anchor classifier fails its held-out check',
 subtitle='Train on a control/KO pair and test the other pair; threshold 0.5. Training labels are themselves assumed, not verified.'),
 'Anchor_cross_validation',12,5)

# Label-blind, gene-balanced grouping: enumerate all ten unique 3-versus-3 partitions.
means <- rowMeans(m);variance <- rowMeans((m-means)^2)
good <- which(variance>=25 & !is.na(x$GENE) & !x$GENE %in% c('','NA','NaN','Spire1','Spire2'))
good <- good[order(-variance[good],x$EVENT[good])];good <- good[!duplicated(x$GENE[good])]
base <- sweep(m[good,,drop=FALSE],1,means[good],'-')
write_table(data.frame(GENE=x$GENE[good],EVENT=x$EVENT[good],PSI_variance=variance[good]),'label_blind_global_features')
partitions <- combn(2:6,2,simplify=FALSE)
ranks <- list();distances <- list()
for (representation in c('raw_PSI','unit_variance_per_gene')) {
 z <- if(representation=='raw_PSI')base else base/sqrt(variance[good])
 pc <- prcomp(t(z),center=FALSE,scale.=FALSE)
 fraction <- pc$sdev^2/sum(pc$sdev^2)
 stopifnot(max(abs(pc$x %*% t(pc$rotation)-t(z)))<1e-7,
   abs(sum(fraction)-1)<1e-10,max(abs(as.matrix(dist(pc$x))-as.matrix(dist(t(z)))))<1e-7)
 coordinates <- data.frame(sample=labels,declared_group=ifelse(seq_along(labels)<4,'Cont','DKO'),pc$x)
 write_table(coordinates,paste0('PCA_coordinates_',representation))
 write_table(data.frame(component=paste0('PC',seq_along(fraction)),variance_fraction=fraction),paste0('PCA_variance_',representation))
 write_table(data.frame(GENE=x$GENE[good],EVENT=x$EVENT[good],pc$rotation),paste0('PCA_loadings_',representation))
 save_plot(ggplot(coordinates,aes(PC1,PC2,color=declared_group,label=sample))+
  geom_hline(yintercept=0,color='grey90')+geom_vline(xintercept=0,color='grey90')+
  geom_point(size=3)+ggrepel::geom_text_repel(seed=42,max.overlaps=Inf,min.segment.length=0)+
  scale_color_manual(values=condition_colors)+coord_equal()+
  scale_x_continuous(expand=expansion(mult=.15))+scale_y_continuous(expand=expansion(mult=.15))+
  labs(title='Spire samples: label-blind splicing PCA',
   subtitle=paste0(format(nrow(z),big.mark=','),' genes; one variable N10 event per gene; Spire1/2 excluded; ',
    if(representation=='raw_PSI')'centered PSI, unscaled' else 'centered PSI, scaled per gene'),
   x=sprintf('PC1 (%.1f%%)',100*fraction[1]),y=sprintf('PC2 (%.1f%%)',100*fraction[2]),color='Original label'),
  paste0('Spire_PCA_',representation),10,7)
 tss <- sum(z^2)
 for (p in partitions) {
  A <- c(1L,p);B <- setdiff(1:6,A)
  wss <- sum(sweep(z[,A,drop=FALSE],1,rowMeans(z[,A,drop=FALSE]),'-')^2)+sum(sweep(z[,B,drop=FALSE],1,rowMeans(z[,B,drop=FALSE]),'-')^2)
  flag <- if(identical(A,1:3))'Original labels' else if(identical(A,c(1L,2L,6L)))'Swap c labels' else 'Other partition'
  ranks[[length(ranks)+1L]]<-data.frame(representation=representation,group_A=paste(labels[A],collapse=' + '),group_B=paste(labels[B],collapse=' + '),assignment=flag,
   genes=nrow(z),variance_explained=1-wss/tss)
 }
 d<-as.matrix(dist(t(z)));dimnames(d)<-list(labels,labels)
 write_table(data.frame(sample=labels,d,check.names=FALSE),paste0('sample_distances_',representation))
 distances[[representation]]<-data.frame(representation=representation,sample=labels,mean_distance=rowSums(d)/5)
}
ranks <- bind_rows(ranks)%>%group_by(representation)%>%mutate(rank=rank(-variance_explained,ties.method='min'))%>%ungroup()
stopifnot(nrow(ranks)==20,sum(ranks$assignment=='Original labels')==2,sum(ranks$assignment=='Swap c labels')==2,all(ranks$variance_explained>=0 & ranks$variance_explained<=1))
write_table(ranks,'all_balanced_partition_scores');write_table(bind_rows(distances),'label_blind_outlier_distances')
save_plot(ggplot(ranks,aes(reorder(group_A,variance_explained),100*variance_explained,fill=assignment))+geom_col()+
 coord_flip()+facet_wrap(~representation,scales='free_y')+
 scale_fill_manual(values=c('Original labels'='#222222','Swap c labels'='#E69F00','Other partition'='grey75'))+
 labs(x='One group; the remaining three samples form the other',y='Variance explained by grouping (%)',
 title='How the original and c-swapped groupings compare',subtitle='All ten balanced partitions; features selected without group labels; ranks are descriptive, not permutation p-values.'),
 'All_balanced_partitions',14,7)

position <- fread(file.path(out,'tables/exon_position_support_by_gene.csv'))
position$sample<-factor(position$sample,levels=samples,labels=labels)
stopifnot(length(unique(position$GENE))==3738,all(table(position$sample)==3738))
save_plot(ggplot(position,aes(sample,log2_five_to_three_support))+geom_boxplot(outlier.shape=NA)+
 geom_hline(yintercept=0,linetype='dashed',color='grey50')+coord_cartesian(ylim=c(-3,3))+
 labs(x=NULL,y='log2(5′ exon support / 3′ exon support)',title='No extreme global 3′ exon-support bias in the c samples',
 subtitle='3,738 genes expressed in all four anchor samples; first/last three observed S exon events; this is not BAM gene-body coverage QC.'),
 'Exon_position_support_QC',12,5)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
cat('Primary markers:',nrow(primary),'\n');print(primary_summary)
cat('Cross-validation agrees:',sum(cv$agrees_with_declared_group),'of',nrow(cv),'held-out samples\n')
print(ranks%>%filter(assignment!='Other partition'))
cat('Diagnostic checks complete; original labels, filtering and simulations untouched.\n')
