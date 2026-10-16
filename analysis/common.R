# Adapted from the oocyte paper scripts; see SOURCES.md and LICENSE.reference.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})
theme_set(theme_classic(base_size = 11) + theme(plot.title = element_text(face = 'bold'),
  strip.background = element_blank(), legend.position = 'bottom'))
direction_colors <- c(up = '#D4A574', down = '#B85450', ndiff = '#808080')
classify <- function(fdr, delta) {
  if (any(!is.finite(fdr) | fdr < 0 | fdr > 1 | !is.finite(delta) | abs(delta) > 1))
    stop('Invalid betAS statistics')
  ifelse(fdr <= .05 & delta >= .1, 'up', ifelse(fdr <= .05 & delta <= -.1, 'down', 'ndiff'))
}
robust_scale <- function(mat, limit) {
  med <- apply(mat, 1, median)
  spread <- apply(mat, 1, mad)
  spread[spread == 0] <- 1
  pmax(pmin((mat - med) / spread, limit), -limit)
}
cluster_rows <- function(mat) {
  if (nrow(mat) == 1) return(list(order = 1L, cluster = 1L))
  stopifnot(nrow(mat) > 1, all(apply(mat, 1, sd) > 0))
  hc <- hclust(as.dist(1 - cor(t(mat), method = 'spearman')), method = 'ward.D2')
  list(order = hc$order, cluster = cutree(hc, k = 2))
}
save_plot <- function(p, name, width = 7, height = 5) {
  if (grepl('^Fig[23]|^GO_paper_OrgDb_all_BP$', name)) main_plots[[name]] <<- p
  p <- p + labs(caption = provenance_label) + theme(plot.caption = element_text(size = 8, color = 'grey35'))
  for (ext in c('pdf', 'png')) ggsave(file.path(out, 'figures', paste0(name, '.', ext)),
    p, device = if (ext == 'pdf') grDevices::cairo_pdf else 'png',
    width = width, height = height, dpi = 300, bg = 'white', limitsize = FALSE)
  figure_names <<- c(figure_names, name)
}
write_table <- function(x, name) fwrite(x, file.path(out, 'tables', paste0(name, '.csv')), na = 'NA')
empty_plot <- function(text) ggplot() + annotate('text', x = 0, y = 0, label = text) + theme_void()
