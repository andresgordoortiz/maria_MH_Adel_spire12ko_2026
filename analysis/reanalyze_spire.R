# Compare the validated FASTQ reanalysis with the preserved historical analysis.
.libPaths(c(normalizePath('analysis/renv/library'), .libPaths()))
source('analysis/common.R')
source('config/analysis.R')
stopifnot(dataset == 'spire_fastq', file.exists('reports/spire_final/COMPLETE'))
out <- 'reports/spire_review'
if (file.exists(out)) stop('Output already exists: ', out)
dir.create(file.path(out, 'tables'), recursive = TRUE)
dir.create(file.path(out, 'figures'))
figure_names <- character(); main_plots <- list()
provenance_label <- 'SPIRE FASTQ reanalysis; original FASTQ labels; 2 November 2026'
old_path <- 'legacy/pre_fastq_2026-11-02/results/spire_final/tables/all_events.csv'
new_path <- 'reports/spire_final/tables/all_events.csv'
old <- as.data.frame(fread(old_path)); new <- as.data.frame(fread(new_path))
samples <- names(config$sample_groups)
old_samples <- grep('^Oocytes_FG_Spire12_', names(old), value = TRUE)
stopifnot(length(old_samples) == 6, !anyDuplicated(old$EVENT), !anyDuplicated(new$EVENT))
common <- merge(old, new, by = 'EVENT', suffixes = c('_old', '_new'), sort = FALSE)
same_annotation <- Reduce(`&`, lapply(c('GENE','COORD','LENGTH','FullCO','COMPLEX'), function(k) {
  x <- common[[paste0(k,'_old')]]; y <- common[[paste0(k,'_new')]]
  (is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & x == y)
}))
write_table(common[!same_annotation, ], 'excluded_annotation_changes')
common <- common[same_annotation, ]
stopifnot(nrow(common) > 0)
write_table(data.frame(EVENT = setdiff(old$EVENT, new$EVENT)), 'old_events_not_tested_now')
write_table(data.frame(EVENT = setdiff(new$EVENT, old$EVENT)), 'new_events_not_tested_before')

# A planted permutation checks the distance orientation and sample matching.
sample_rmse <- function(a, b) t(vapply(seq_len(ncol(a)), function(i)
  sqrt(colMeans((b - a[,i])^2)), numeric(ncol(b))))
fixture <- cbind(c(0,20,60),c(80,10,40),c(30,90,0))
stopifnot(identical(max.col(-sample_rmse(fixture, fixture[,c(3,1,2)])), c(2L,3L,1L)))
a <- as.matrix(common[,old_samples]); b <- as.matrix(common[,samples])
variance <- apply(cbind(a,b), 1, var)
ix <- which(variance >= 25 & common$GENE_new %in% valid_genes(common$GENE_new) &
  !common$GENE_new %in% c('Spire1','Spire2'))
ix <- ix[order(-variance[ix],common$EVENT[ix])]
ix <- ix[!duplicated(common$GENE_new[ix])]
stopifnot(length(ix) > 100)
write_table(common[ix,c('EVENT','GENE_new')], 'sample_matching_features')
distances <- sample_rmse(a[ix,], b[ix,]); correlations <- cor(a[ix,], b[ix,], method = 'spearman')
best <- max.col(-distances)
stopifnot(!anyDuplicated(best), identical(best,max.col(correlations)))
matching <- data.frame(old_sample = old_samples, new_sample = samples[best],
  old_group = rep(c('Cont','DKO'), each = 3), new_group = unname(config$sample_groups[best]),
  genes = length(ix), RMSE_pp = distances[cbind(1:6,best)],
  Spearman = correlations[cbind(1:6,best)],
  next_best_RMSE_pp = apply(distances,1,function(x) sort(x)[2]))
matching$condition_changed <- matching$old_group != matching$new_group
write_table(matching, 'sample_correspondence')
pairs <- expand.grid(old_sample = old_samples, new_sample = samples)
pairs$RMSE_pp <- as.vector(distances); pairs$Spearman <- as.vector(correlations)
write_table(pairs, 'all_sample_pairs')
save_plot(ggplot(pairs,aes(sub('^SPIRE_','',new_sample),sub('Oocytes_FG_Spire12_','',old_sample),fill=Spearman))+
  geom_tile()+geom_text(aes(label=sprintf('%.4f',Spearman)),size=3)+
  scale_fill_viridis_c()+labs(x='New FASTQ sample',y='Historical label',
    title='Historical and FASTQ sample correspondence',subtitle=paste(length(ix),'variable genes; one event per gene; Spire1/2 excluded')),
  'Sample_correspondence',9,5)

common$change <- ifelse(common$direction_old == 'ndiff',
  ifelse(common$direction_new == 'ndiff','Neither','New only'),
  ifelse(common$direction_new == 'ndiff','Old only',
    ifelse(common$direction_old == common$direction_new,'Retained same direction','Reversed direction')))
common$raw_delta_old <- (rowMeans(a[,4:6])-rowMeans(a[,1:3]))/100
common$raw_delta_new <- (rowMeans(b[,4:6])-rowMeans(b[,1:3]))/100
write_table(common, 'old_new_common_events')
write_table(common[common$change != 'Neither', ], 'changed_or_retained_selected_events')
changes <- as.data.frame(table(common$change)); names(changes) <- c('change','events')
write_table(changes, 'selection_changes')
counts <- bind_rows(lapply(list(historical=old,FASTQ=new),function(x)
  data.frame(tested=nrow(x),selected=sum(x$direction!='ndiff'),
    genes=length(valid_genes(x$GENE[x$direction!='ndiff'])),up=sum(x$direction=='up'),down=sum(x$direction=='down'))),.id='dataset')
write_table(counts, 'splicing_summary')
write_table(bind_rows(lapply(list(historical=old,FASTQ=new),function(x)
  x %>% count(event_type,direction)),.id='dataset'), 'splicing_counts_by_type')
save_plot(ggplot(common,aes(100*deltapsi_old,100*deltapsi_new,color=change))+
  geom_point(size=.5,alpha=.5)+geom_hline(yintercept=0,color='grey60')+
  geom_vline(xintercept=0,color='grey60')+labs(x='Historical fitted delta PSI (pp)',y='New fitted delta PSI (pp)',
    title='Differential splicing before and after reprocessing',color=NULL), 'Old_new_splicing',9,6)

plot_pca <- function(m,name,title) {
  pc <- prcomp(t(m),center=TRUE,scale.=FALSE)
  fraction <- pc$sdev^2/sum(pc$sdev^2)
  coords <- data.frame(sample=colnames(m),group=unname(config$sample_groups[colnames(m)]),pc$x)
  write_table(coords,paste0(name,'_coordinates'))
  write_table(data.frame(component=seq_along(fraction),fraction=fraction),paste0(name,'_variance'))
  save_plot(ggplot(coords,aes(PC1,PC2,color=group,label=sub('^SPIRE_','',sample)))+
    geom_point(size=3)+ggrepel::geom_text_repel(seed=42,max.overlaps=Inf)+
    scale_color_manual(values=condition_colors)+labs(title=title,color='FASTQ group',
      x=sprintf('PC1 (%.1f%%)',100*fraction[1]),y=sprintf('PC2 (%.1f%%)',100*fraction[2])),name,9,6)
}
plot_pca(b[ix,],'Splicing_PCA','SPIRE splicing PCA: one variable event per gene')

dge_dir <- 'data/inputs/spire_fastq/dge/tables'
dge_path <- file.path(dge_dir,'differential/rnaseq/DKO_vs_Cont_SPIRE_DKO_deseq2.annotated.tsv')
dge <- as.data.frame(fread(dge_path))
stopifnot(!anyDuplicated(dge$gene_id),all(is.na(dge$padj) | (dge$padj>=0 & dge$padj<=1)))
dge$selected <- !is.na(dge$padj) & dge$padj<=.05
write_table(dge, 'DGE_all_genes')
write_table(dge[dge$selected, ][order(dge$padj[dge$selected]), ], 'DGE_padj005')
dge_summary <- data.frame(input_genes=nrow(dge),finite_padj=sum(is.finite(dge$padj)),
  padj005=sum(dge$selected),up=sum(dge$selected & dge$log2FoldChange>0),
  down=sum(dge$selected & dge$log2FoldChange<0),
  padj005_abs_log2FC_ge1=sum(dge$selected & abs(dge$log2FoldChange)>=1))
write_table(dge_summary, 'DGE_summary')
vst_path <- file.path(dge_dir,'processed_abundance/rnaseq/all.vst.tsv')
vst <- as.data.frame(fread(vst_path)); v <- as.matrix(vst[,samples])
stopifnot(!anyDuplicated(vst$gene_id),all(is.finite(v)))
eligible <- which(!vst$gene_id %in% dge$gene_id[dge$gene_name %in% c('Spire1','Spire2')])
top <- head(eligible[order(apply(v[eligible,],1,var),decreasing=TRUE)],500)
write_table(data.frame(gene_id=vst$gene_id[top]),'Expression_PCA_features')
plot_pca(v[top,],'Expression_PCA','SPIRE expression PCA: 500 variable VST genes, Spire1/2 excluded')
save_plot(ggplot(dge,aes(log2FoldChange,-log10(pmax(padj,1e-300)),color=selected))+
  geom_point(size=.7,alpha=.5,na.rm=TRUE)+geom_hline(yintercept=-log10(.05),linetype='dashed')+
  scale_color_manual(values=c('FALSE'='grey65','TRUE'=condition_colors[['DKO']]))+
  labs(x='Exported log2 fold change (DKO / Cont)',y='−log10 adjusted p',color='padj ≤0.05',title='New SPIRE differential expression'),
  'DGE_volcano',9,6)
norm_path <- file.path(dge_dir,'processed_abundance/rnaseq/all.normalised_counts.tsv')
norm <- as.data.frame(fread(norm_path)); stopifnot(!anyDuplicated(norm$gene_id))
target <- dge[dge$gene_name %in% c('Spire1','Spire2','Fmn2'), ]
write_table(target,'Target_gene_DGE')
expression <- merge(target[,c('gene_id','gene_name')],norm,by='gene_id') %>%
  pivot_longer(all_of(samples),names_to='sample',values_to='normalized_count') %>%
  mutate(group=unname(config$sample_groups[sample]))
write_table(expression,'Target_gene_expression')
save_plot(ggplot(expression,aes(sub('^SPIRE_','',sample),normalized_count+1,color=group))+
  geom_point(size=3)+facet_wrap(~gene_name,scales='free_y')+scale_y_log10()+
  scale_color_manual(values=condition_colors)+theme(axis.text.x=element_text(angle=45,hjust=1))+
  labs(x=NULL,y='Normalized counts + 1 (log scale)',title='Target loci in the original FASTQ groups',color=NULL),
  'Target_gene_expression',12,5)
write_table(merge(new[new$GENE %in% c('Spire1','Spire2','Fmn2'), ],
  old[,c('EVENT','deltapsi','FDR','direction')],by='EVENT',all.x=TRUE,suffixes=c('_new','_old')), 'Target_locus_splicing')
raw <- fread(config$input,select=c('GENE','EVENT','COORD','COMPLEX',samples,paste0(samples,'-Q')))
raw <- as.data.frame(raw[GENE %in% c('Spire1','Spire2','Fmn2') & COMPLEX %in% c('S','C1','C2','C3','MIC')])
corrected_total <- function(q) vapply(strsplit(sub('.*@','',q),',',fixed=TRUE),function(x) {
  y <- as.numeric(x); stopifnot(length(y)==2)
  if (any(!is.finite(y) | y<0)) NA_real_ else sum(y)
},numeric(1))
stopifnot(identical(corrected_total(c('N,N@3.25,6.75','N,N@-1,-1')),c(10,NA_real_)))
exons <- bind_rows(lapply(samples,function(s) data.frame(raw[,1:4],sample=s,
  group=config$sample_groups[[s]],PSI=raw[[s]],corrected_total=corrected_total(raw[[paste0(s,'-Q')]]))))
write_table(exons,'Target_exon_support_all_samples')
exon_summary <- exons %>% group_by(GENE,EVENT,COORD) %>% summarise(
  min_control_support=min(corrected_total[group=='Cont']),max_DKO_support=max(corrected_total[group=='DKO']),.groups='drop') %>%
  mutate(supported_all_controls_lost_all_DKO=min_control_support>=10 & max_DKO_support<10)
write_table(exon_summary,'Target_exon_support_summary')
save_plot(ggplot(exons,aes(factor(sub('^SPIRE_','',sample),levels=sub('^SPIRE_','',samples)),EVENT,fill=log10(corrected_total+1)))+
  geom_tile()+facet_wrap(~GENE,scales='free_y')+scale_fill_viridis_c(na.value='grey80')+
  theme(axis.text.x=element_text(angle=45,hjust=1),axis.text.y=element_text(size=7))+
  labs(x=NULL,y=NULL,fill='log10 support + 1',title='Target exon support, including events excluded by all-sample N10',
    subtitle='VAST corrected junction support; overlapping events are not gene expression counts'),
  'Target_exon_support',14,7)
qc_path <- 'data/inputs/spire_fastq/rnaseq/multiqc/star_salmon/multiqc_report_data/multiqc_star.txt'
qc <- fread(qc_path); stopifnot(setequal(qc$Sample,samples))
write_table(qc,'RNAseq_alignment_QC')

source_files <- c(old_path,new_path,dge_path,vst_path,norm_path,qc_path,config$input,'config/spire_fastq.tsv',
  'config/spire_contrasts.csv','analysis/reanalyze_spire.R')
write_table(data.frame(file=source_files,sha256=vapply(source_files,function(p)
  digest::digest(file=p,algo='sha256'),character(1))), 'source_checksums')
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
report <- c('# Reanalisis de SPIRE desde los FASTQ originales','',
  '2 de noviembre de 2026. Se conservan las etiquetas de los FASTQ; no se modifica ninguna tabla original.','',
  '## Correspondencia entre muestras','',
  sprintf('Emparejamiento independiente de la condicion con %s genes variables (un evento por gen, sin Spire1/2). Los seis mejores pares son univocos y coinciden usando RMSE y correlacion de Spearman.',length(ix)), '',
  '| Etiqueta antigua | Muestra FASTQ mas parecida | Spearman | RMSE (pp) | Cambio de condicion |',
  '|---|---|---:|---:|---|',
  vapply(seq_len(6),function(i) sprintf('| %s | %s | %.5f | %.3f | %s |',
    sub('Oocytes_FG_Spire12_','',matching$old_sample[i]),matching$new_sample[i],matching$Spearman[i],matching$RMSE_pp[i],matching$condition_changed[i]),character(1)), '',
  'La correspondencia identifica una inconsistencia de etiquetas entre tablas. No constituye un genotipado de ADN ni demuestra en que paso se produjo el cambio.','',
  '## Splicing','',
  sprintf('Antes: %s eventos evaluados, %s seleccionados en %s genes. Ahora: %s evaluados, %s seleccionados en %s genes.',
    counts$tested[1],counts$selected[1],counts$genes[1],counts$tested[2],counts$selected[2],counts$genes[2]),'',
  sprintf('La comparacion usa %s eventos compartidos con anotacion identica; se excluyen %s discrepancias de anotacion.',nrow(common),sum(!same_annotation)), '',
  paste0('- ',changes$change,': ',changes$events), '',
  'Mismos criterios: N>=10 y PSI completa en las seis muestras; FPR empirico betAS<=0.05 y |delta PSI ajustado|>=0.10. El campo FDR de betAS es FPR por evento, no FDR corregida por pruebas multiples. Los cambios mezclan el nuevo procesamiento con la asignacion correcta de condiciones; no se atribuyen exclusivamente al etiquetado.','',
  '## Expresion diferencial y QC','',
  sprintf('DGE suministrado, contraste DKO/Cont: %s genes con padj<=0.05 (%s aumentan, %s disminuyen); %s cumplen ademas |log2FC|>=1. No se ha vuelto a ejecutar DESeq2.',
    dge_summary$padj005,dge_summary$up,dge_summary$down,dge_summary$padj005_abs_log2FC_ge1),'',
  vapply(seq_len(nrow(target)),function(i) sprintf('- %s: log2FC=%+.3f; padj=%.4g.',target$gene_name[i],target$log2FoldChange[i],target$padj[i]),character(1)), '',
  vapply(c('Spire1','Spire2'),function(g) sprintf('- %s: %s de %s eventos exonicos observados tienen soporte>=10 en todos los controles y <10 en todos los DKO. Se revisan antes del filtro N10, que los excluye del contraste de splicing.',g,
    sum(exon_summary$supported_all_controls_lost_all_DKO[exon_summary$GENE==g],na.rm=TRUE),sum(exon_summary$GENE==g)),character(1)), '',
  sprintf('STAR: %.2f-%.2f%% de lecturas alineadas de forma unica. PCA de expresion sobre los 500 genes VST mas variables, sin Spire1/2; PCA de splicing sobre los genes usados para emparejar, sin usar seleccion diferencial.',min(qc$uniquely_mapped_percent),max(qc$uniquely_mapped_percent)), '',
  'El archivo normalizado permite revisar cada muestra. La expresion total de Spire1 no valida ni invalida por si sola una delecion de exones: el transcrito residual puede seguir siendo cuantificable. El investigador confirmó el 2 de noviembre de 2026 que las seis muestras proceden de animales independientes. Los nombres de FASTQ no establecen efectos de lote. Cont_7 lleva fecha 14mars; las otras cinco muestras, 12mars.','',
  '[Figuras](index.html) · [Splicing completo](../spire_final/index.html) · [Comparacion actual con FMN2](../comparison/index.html) · [Historico](../../legacy/pre_fastq_2026-11-02/results/index.html)','')
writeLines(report,file.path(out,'REPORT.md'))
html <- c('<!doctype html><html lang="es"><meta charset="utf-8"><title>Revision SPIRE FASTQ</title>',
  '<style>body{font:17px system-ui;max-width:1100px;margin:30px auto;padding:0 20px}img{max-width:100%}article{margin-top:40px}</style>',
  '<h1>SPIRE: reanalisis desde FASTQ</h1><p><a href="REPORT.md">Resultados y limites</a> · <a href="tables/sample_correspondence.csv">Correspondencia de muestras</a> · <a href="tables/selection_changes.csv">Cambios de seleccion</a> · <a href="../index.html">Inicio</a></p>',
  sprintf('<p>%s eventos seleccionados antes; %s ahora. DGE: %s genes con padj ≤0.05.</p>',counts$selected[1],counts$selected[2],dge_summary$padj005))
for (name in figure_names) html <- c(html,sprintf('<article><h2>%s</h2><a href="figures/%s.pdf">PDF</a><img src="figures/%s.png" alt="%s" loading="lazy"></article>',name,name,name,name))
writeLines(c(html,'</html>'),file.path(out,'index.html'))
stopifnot(all(file.exists(file.path(out,'figures',paste0(figure_names,'.pdf')))))
writeLines(provenance_label,file.path(out,'COMPLETE'))
cat(paste(report,collapse='\n'),'\n')
