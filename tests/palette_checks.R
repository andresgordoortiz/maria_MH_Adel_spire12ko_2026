.libPaths(c(normalizePath('analysis/renv/library'), .libPaths()))
source('analysis/common.R')
stopifnot(identical(names(dataset_colors), c('spire', 'fmn2')),
  identical(names(direction_colors), c('up', 'down', 'ndiff')),
  signed_colors[1] == direction_colors['down'], signed_colors[3] == direction_colors['up'],
  condition_colors['Cont'] == condition_colors['Het'], condition_colors['DKO'] == condition_colors['KO'],
  overlap_colors['Spire only'] == dataset_colors['spire'], overlap_colors['FMN2 only'] == dataset_colors['fmn2'],
  impact_colors['ORF disrupted upon inclusion'] == direction_colors['up'],
  impact_colors['ORF disrupted upon exclusion'] == direction_colors['down'])
# High magnitudes have lower luminance, including when printed in greyscale.
rgb <- grDevices::col2rgb(viridisLite::cividis(32, direction = -1)) / 255
linear <- ifelse(rgb <= .04045, rgb / 12.92, ((rgb + .055) / 1.055)^2.4)
luminance <- drop(c(.2126, .7152, .0722) %*% linear)
stopifnot(all(diff(luminance) < 0))
cat('PASS: shared colour meanings and monotonic sequential luminance.\n')
