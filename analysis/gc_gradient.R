# GC binning functions from the pinned paper repository (MIT; LICENSE.reference).
# ---- bin layouts ----
# Exon layout (5 regions, 68 bins total)
EXON_LAYOUT <- list(
  upExon    = list(end = "right", n_bins =  8, pad_nt = 200),
  upIntron  = list(mode = "full",  n_bins = 20),
  target    = list(mode = "full",  n_bins = 12),
  doIntron  = list(mode = "full",  n_bins = 20),
  downExon  = list(end = "left",  n_bins =  8, pad_nt = 200)
)

# Two intron layouts — picked by sample size (see file header).
make_intron_layout <- function(level = c("fine", "coarse")) {
  if (level == "fine") {
    list(
      upExon        = list(end = "right", n_bins = 16, pad_nt = 400),
      upIntron5SS   = list(end = "left",  n_bins =  8, pad_nt =  80),
      intronBody    = list(mode = "full",  n_bins = 40),
      doIntron3SS   = list(end = "right", n_bins =  8, pad_nt =  80),
      downExon      = list(end = "left",  n_bins = 16, pad_nt = 400)
    )
  } else {
    list(
      upExon        = list(end = "right", n_bins =  6, pad_nt = 150),
      upIntron5SS   = list(end = "left",  n_bins =  3, pad_nt =  80),
      intronBody    = list(mode = "full",  n_bins = 12),
      doIntron3SS   = list(end = "right", n_bins =  3, pad_nt =  80),
      downExon      = list(end = "left",  n_bins =  6, pad_nt = 150)
    )
  }
}

# Region order used by the plot facets
region_order_exon  <- c("upExon", "upIntron", "Exon", "downIntron", "downExon")
region_order_intron <- c("upExon", "Intron 5' flank", "Intron body",
                         "Intron 3' flank", "downExon")

is_empty_str <- function(s) {
  if (length(s) == 0) return(TRUE)
  if (is.na(s[1]))   return(TRUE)
  nchar(s[1]) == 0
}

# ---- per-bin GC ----
# Average GC of a character vector of nucleotides (skipping N).
gc_of_vec <- function(s) {
  if (length(s) == 0) return(NA_real_)
  s <- s[s %in% c("A","T","G","C")]
  if (length(s) == 0) return(NA_real_)
  100 * sum(s %in% c("G","C")) / length(s)
}

# Break a DNA string into exactly n_bins contiguous chunks (as equal in
# size as possible) and return the GC of each chunk. Returns NA for
# chunks with no sequence.
bin_to_gc <- function(seq, n_bins) {
  if (is_empty_str(seq)) return(rep(NA_real_, n_bins))
  chars <- strsplit(toupper(seq[1]), "")[[1]]
  chars <- chars[chars %in% c("A","T","G","C")]
  n <- length(chars)
  if (n == 0) return(rep(NA_real_, n_bins))
  if (n < n_bins) {
    out <- rep(NA_real_, n_bins)
    out[seq_len(n)] <- vapply(seq_len(n), function(i) gc_of_vec(chars[i]), numeric(1))
    return(out)
  }
  sizes <- rep(floor(n / n_bins), n_bins)
  sizes[seq_len(n %% n_bins)] <- sizes[seq_len(n %% n_bins)] + 1
  ends   <- cumsum(sizes)
  starts <- c(1, head(ends, -1) + 1)
  unname(mapply(function(s, e) gc_of_vec(chars[s:e]), starts, ends))
}

# Take an end-slice of a DNA string (the LAST n nt, or the FIRST n nt).
left_end  <- function(seq, n) {
  if (is_empty_str(seq)) return(NA_character_)
  L <- nchar(seq)
  if (L <= n) return(seq)
  substr(seq, L - n + 1, L)
}
right_end <- function(seq, n) {
  if (is_empty_str(seq)) return(NA_character_)
  L <- nchar(seq)
  if (L <= n) return(seq)
  substr(seq, 1, n)
}

# Peal off n_each_end nucleotides from each end of an intron (used to
# recover the body after the 5'/3' flanks have been binned separately).
strip_intron_ends <- function(seq, n_each_end) {
  if (is_empty_str(seq)) return(NA_character_)
  L <- nchar(seq)
  if (L <= 2 * n_each_end) return(NA_character_)
  substr(seq, n_each_end + 1, L - n_each_end)
}

# Per-event gradient builders: each returns the per-bin GC vector for
# that event (concatenated across the regions).
build_event_exon_gradient <- function(seq_upExon, seq_upIntron, seq_exon,
                                      seq_doIntron, seq_doExon) {
  c(bin_to_gc(left_end(seq_upExon,   EXON_LAYOUT$upExon$pad_nt),
              EXON_LAYOUT$upExon$n_bins),
    bin_to_gc(seq_upIntron,  EXON_LAYOUT$upIntron$n_bins),
    bin_to_gc(seq_exon,      EXON_LAYOUT$target$n_bins),
    bin_to_gc(seq_doIntron,  EXON_LAYOUT$doIntron$n_bins),
    bin_to_gc(right_end(seq_doExon, EXON_LAYOUT$downExon$pad_nt),
              EXON_LAYOUT$downExon$n_bins))
}

build_event_intron_gradient <- function(seq_upExon, seq_intron, seq_doExon,
                                        layout) {
  peel <- layout$upIntron5SS$pad_nt
  c(bin_to_gc(left_end(seq_upExon, layout$upExon$pad_nt),
              layout$upExon$n_bins),
    # Correct the reference script's reversed intron 5'/3' end slices.
    bin_to_gc(right_end(seq_intron, peel),
              layout$upIntron5SS$n_bins),
    bin_to_gc(strip_intron_ends(seq_intron, peel),
              layout$intronBody$n_bins),
    bin_to_gc(left_end(seq_intron, peel),
              layout$doIntron3SS$n_bins),
    bin_to_gc(right_end(seq_doExon, layout$downExon$pad_nt),
              layout$downExon$n_bins))
}

# Region labels for the per-bin index produced by the gradient builders.
x_anchors_exon <- c(rep("upExon",   EXON_LAYOUT$upExon$n_bins),
                    rep("upIntron", EXON_LAYOUT$upIntron$n_bins),
                    rep("Exon",     EXON_LAYOUT$target$n_bins),
                    rep("downIntron",EXON_LAYOUT$doIntron$n_bins),
                    rep("downExon", EXON_LAYOUT$downExon$n_bins))

x_anchors_intron_from_layout <- function(layout) {
  c(rep("upExon",          layout$upExon$n_bins),
    rep("Intron 5' flank", layout$upIntron5SS$n_bins),
    rep("Intron body",     layout$intronBody$n_bins),
    rep("Intron 3' flank", layout$doIntron3SS$n_bins),
    rep("downExon",        layout$downExon$n_bins))
}

for (kind in c('EX', 'IR')) {
  df <- matt_tables[[kind]]
  bg <- which(df$DATASET == 'ndiff')
  set.seed(42)
  if (length(bg) > 1500) bg <- sample(bg, 1500)
  df <- df[c(which(df$DATASET != 'ndiff'), bg), ]
  write_table(df[, c('EVENT', 'DATASET')], paste0('GC_gradient_', kind, '_event_manifest'))
  for (level in if (kind == 'IR') c('fine', 'coarse') else 'fine') {
    if (kind == 'EX') {
      required <- c('SEQ_UPEXON', 'SEQ_UPINTRON', 'SEQ_EXON', 'SEQ_DOINTRON', 'SEQ_DOEXON')
      anchors <- x_anchors_exon
      build <- function(i) do.call(build_event_exon_gradient, unname(as.list(df[i, required])))
    } else {
      required <- c('SEQ_LONGESTUPEXON', 'SEQ_INTRON', 'SEQ_LONGESTDOEXON')
      layout <- make_intron_layout(level)
      anchors <- x_anchors_intron_from_layout(layout)
      build <- function(i) do.call(build_event_intron_gradient, c(unname(as.list(df[i, required])), list(layout = layout)))
    }
    stopifnot(all(required %in% names(df)))
    values <- t(vapply(seq_len(nrow(df)), build, numeric(length(anchors))))
    long <- data.frame(EVENT = rep(df$EVENT, each = length(anchors)),
      direction = rep(df$DATASET, each = length(anchors)),
      bin = rep(seq_along(anchors), times = nrow(df)),
      region = factor(rep(anchors, times = nrow(df)), levels = unique(anchors)), GC = as.vector(t(values)))
    long <- long[is.finite(long$GC), ]
    for (view in c('directions', 'differential_vs_ndiff')) {
      long$category <- if (view == 'directions') long$direction else ifelse(long$direction == 'ndiff', 'ndiff', 'differential')
      agg <- long %>% group_by(category, bin, region) %>% summarise(mean_GC = mean(GC),
        se = sd(GC) / sqrt(n()), n = n(), .groups = 'drop')
      write_table(agg, paste('GC_gradient', kind, level, view, sep = '_'))
      colors <- c(direction_colors, differential = '#222222')
      save_plot(ggplot(agg, aes(bin, mean_GC, color = category, fill = category)) +
        geom_ribbon(aes(ymin = mean_GC - 1.96 * se, ymax = mean_GC + 1.96 * se), alpha = .15, color = NA) +
        geom_line(linewidth = .6) + facet_grid(~region, scales = 'free_x', space = 'free_x') +
        scale_color_manual(values = colors) + scale_fill_manual(values = colors) +
        labs(x = "Proportional bins along the 5′ → 3′ sequence", y = 'GC (%)',
          title = paste(kind, 'GC sequence gradient'), subtitle = paste(level, 'layout; mean ±1.96 SE; ndiff capped at 1,500')) +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()),
        paste('Supplement_GC_gradient', kind, level, view, sep = '_'), 11, 4)
    }
  }
}
