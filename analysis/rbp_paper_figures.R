# Adapted from Figure2H_S4A_B.R and FigureS4C_E.R; see SOURCES.md.
registry <- data.frame(figure=character(),paper_panel=character(),type=character())
record_panel <- function(plot,name,panel,kind,width,height) {
  save_plot(plot,name,width,height)
  registry <<- rbind(registry,data.frame(figure=name,paper_panel=panel,type=kind))
}
for (kind in c('EX','IR','Alt5','Alt3')) {
  regions <- rbp_regions(kind)
  df <- rbp_panel_data(motifs,kind)
  outcome_colors <- setNames(unname(direction_colors[c('up','down')]),levels(df$outcome))
  write_table(df,paste0('RBP_',kind,'_paper_motif_scores'))
  top <- df %>% group_by(dataset,direction,RBP) %>% summarise(score=max(score),.groups='drop') %>%
    group_by(dataset,direction) %>% slice_max(score,n=15) %>% ungroup()
  tracks <- semi_join(df,top,by=c('dataset','direction','RBP'))
  labels <- tracks %>% group_by(dataset,direction,outcome,region,mid,RBP) %>%
    summarise(score=max(score),.groups='drop') %>% group_by(dataset,direction,region) %>%
    slice_max(score,n=3,with_ties=FALSE) %>% ungroup()
  write_table(top,paste0('RBP_',kind,'_track_selection'))
  p <- ggplot(tracks,aes(mid,score,color=outcome,fill=outcome))+
    geom_point(size=1.5,alpha=.45,position=position_jitter(width=8,height=0,seed=42))+
    # The original .3 span is singular with only five/eight distinct region positions.
    geom_smooth(method='loess',formula=y~x,span=.75,alpha=.15,linewidth=1,se=TRUE)+
    geom_vline(xintercept=head(regions$end,-1),linetype='dotted',color='grey65',linewidth=.25)+
    ggrepel::geom_text_repel(data=labels,aes(label=RBP),size=2.4,fontface='bold',
      max.overlaps=Inf,seed=42,min.segment.length=.1,box.padding=.25,point.padding=.2)+
    geom_hline(yintercept=-log10(.05),linetype='dashed',color='grey45',linewidth=.35)+
    facet_wrap(~dataset,ncol=2,scales='free_y')+
    scale_x_continuous(breaks=regions$mid,labels=regions$label,expand=expansion(mult=.03),guide=guide_axis(n.dodge=2))+
    scale_color_manual(values=outcome_colors)+scale_fill_manual(values=outcome_colors)+
    labs(x=NULL,y='−log10 raw motif p',color='Event usage',fill='Event usage',
      title=paste(kind,'RBP motif scores across RNA regions'),
      subtitle='Top 15 labels per dataset/direction (ties retained); top 3 labelled per region; schematic regions')+
    theme(axis.text.x=element_text(size=7),legend.position='top',panel.spacing=grid::unit(1.2,'lines'))
  name <- switch(kind,EX='FigS6A_EX_RBP_tracks',IR='FigS6C_IR_RBP_tracks',paste0('Supplement_',kind,'_RBP_tracks'))
  record_panel(p,name,switch(kind,EX='S6A',IR='S6C','Repository extension'),kind,18,7)

  tests <- rbp_box_tests(df)
  write_table(tests,paste0('RBP_',kind,'_direction_wilcoxon'))
  p <- ggplot(df,aes(dataset,score,fill=outcome))+
    geom_boxplot(alpha=.7,outlier.shape=NA,linewidth=.4,position=position_dodge(.8))+
    geom_point(aes(color=outcome),alpha=.25,size=.6,
      position=position_jitterdodge(jitter.width=.16,dodge.width=.8,seed=42))+
    geom_text(data=tests,aes(dataset,y,label=label),inherit.aes=FALSE,size=3)+
    geom_hline(yintercept=-log10(.05),linetype='dashed',color='grey55',linewidth=.3)+
    facet_wrap(~region,nrow=1,scales='free_y',labeller=as_labeller(setNames(regions$label,regions$region)))+
    scale_fill_manual(values=outcome_colors)+scale_color_manual(values=outcome_colors)+
    labs(x=NULL,y='−log10 raw motif p',fill='Event usage',color='Event usage',
      title=paste(kind,'motif-score distributions by direction'),
      subtitle='Paper selection: motifs with raw p ≤0.05 anywhere; raw Wilcoxon stars compare score distributions')+
    theme(legend.position='top',axis.text.x=element_text(angle=45,hjust=1),strip.text=element_text(size=8))
  name <- switch(kind,EX='FigS6B_EX_RBP_boxplots',IR='FigS6E_IR_RBP_boxplots',paste0('Supplement_',kind,'_RBP_boxplots'))
  record_panel(p,name,switch(kind,EX='S6B',IR='S6E','Repository extension'),kind,18,6)

  collapsed <- df %>% group_by(RBP,dataset,direction,region) %>%
    summarise(pvalue=min(pvalue),score=max(score),.groups='drop')
  top <- collapsed %>% group_by(RBP) %>% summarise(score=max(score),.groups='drop') %>%
    arrange(desc(score),RBP) %>% slice_head(n=if (kind=='IR') 15 else 30) %>% pull(RBP)
  h <- expand_grid(RBP=top,dataset=c('spire','fmn2'),direction=c('up','down'),region=regions$region) %>%
    left_join(collapsed,by=c('RBP','dataset','direction','region')) %>%
    mutate(score=replace_na(score,0),pvalue=replace_na(pvalue,1),
      column=paste(dataset,direction,region,sep=':'))
  columns <- unlist(lapply(regions$region,function(r)
    paste(rep(c('spire','fmn2'),each=2),rep(c('up','down'),2),r,sep=':')))
  wide <- h %>% select(RBP,column,score) %>% pivot_wider(names_from=column,values_from=score)
  matrix <- as.matrix(as.data.frame(wide)[,columns,drop=FALSE]);rownames(matrix) <- wide$RBP
  order <- if (nrow(matrix)>1) rownames(matrix)[hclust(dist(matrix),method='complete')$order] else rownames(matrix)
  h$RBP <- factor(h$RBP,levels=order);h$column <- factor(h$column,levels=columns)
  h$region <- factor(h$region,levels=regions$region)
  write_table(h,paste0('RBP_',kind,'_paper_heatmap_scores'))
  ann <- h %>% distinct(dataset,direction,region,column)
  ann$outcome <- levels(df$outcome)[match(ann$direction,c('up','down'))]
  annotations <- bind_rows(transmute(ann,column,region,track='Dataset',value=dataset),
    transmute(ann,column,region,track='Event usage',value=outcome))
  annotations$track <- factor(annotations$track,levels=c('Event usage','Dataset'))
  annotation_plot <- ggplot(annotations,aes(column,track,fill=value))+
    geom_tile(color='white',linewidth=.1)+facet_grid(~region,scales='free_x',space='free_x',
      labeller=as_labeller(setNames(regions$label,regions$region)))+
    scale_fill_manual(values=c(dataset_colors,outcome_colors),name=NULL)+
    theme_void()+theme(strip.text=element_text(size=8),axis.text.y=element_text(size=8),legend.position='top')
  heat <- ggplot(h,aes(column,RBP,fill=score))+geom_tile(color='white',linewidth=.15)+
    facet_grid(~region,scales='free_x',space='free_x')+
    scale_fill_viridis_c(option=sequential_palette,direction=-1,
      limits=c(0,8),oob=scales::squish,name='−log10 raw p')+
    labs(x='Columns within each region: Spire up/down, FMN2 up/down',y=NULL)+
    theme(axis.text.x=element_blank(),axis.ticks=element_blank(),axis.line=element_blank(),
      strip.text=element_blank(),panel.spacing.x=grid::unit(.15,'lines'),legend.position='right')
  p <- patchwork::wrap_plots(annotation_plot,heat,ncol=1,heights=c(1,6))+
    patchwork::plot_annotation(title=paste(kind,'RBP motif scores across both datasets'),
      subtitle='Paper layout: region → dataset → direction; Euclidean/complete row clustering; colour saturates at 8')
  name <- switch(kind,EX='Fig2H_EX_RBP_heatmap',IR='FigS6D_IR_RBP_heatmap',paste0('Supplement_',kind,'_RBP_heatmap'))
  record_panel(p,name,switch(kind,EX='2H',IR='S6D','Repository extension'),kind,16,if (kind=='IR') 7 else 10)

  counts <- minima %>% filter(type==kind,pvalue<=.05) %>% count(dataset,direction,region,name='n_RBP_labels') %>%
    complete(dataset=c('spire','fmn2'),direction=c('up','down'),region=regions$region,fill=list(n_RBP_labels=0))
  counts$region <- factor(counts$region,levels=regions$region)
  write_table(counts,paste0('RBP_',kind,'_region_counts'))
  p <- ggplot(counts,aes(region,n_RBP_labels,fill=dataset,alpha=direction))+
    geom_col(position=position_dodge(.8))+
    scale_fill_manual(values=dataset_colors)+scale_alpha_manual(values=c(up=1,down=.6))+
    scale_x_discrete(labels=setNames(regions$label,regions$region))+
    labs(x=NULL,y='RBP labels with raw p ≤0.05',title=paste(kind,'RBP label counts by region'))+
    theme(legend.position='top',axis.text.x=element_text(size=8))
  record_panel(p,paste0('Supplement_',kind,'_RBP_region_counts'),'Repository count summary',kind,14,5)
  totals <- minima %>% filter(type==kind,pvalue<=.05) %>% group_by(dataset,direction) %>%
    summarise(n_RBP_labels=n_distinct(RBP),.groups='drop')
  write_table(totals,paste0('RBP_',kind,'_total_counts'))
  p <- ggplot(totals,aes(dataset,n_RBP_labels,fill=direction))+
    geom_col(position='dodge')+scale_fill_manual(values=direction_colors)+
    labs(x=NULL,y='RBP labels with raw p ≤0.05 anywhere',title=paste(kind,'RBP label counts across regions'))
  record_panel(p,paste0('Supplement_',kind,'_RBP_total_counts'),'Repository count summary',kind,7,5)
}
write_table(registry,'paper_panel_manifest')
