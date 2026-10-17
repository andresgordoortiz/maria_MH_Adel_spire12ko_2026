# Okabe–Ito categories; Cividis magnitudes; blue–neutral–orange signed effects.
dataset_colors <- c(spire = '#CC79A7', fmn2 = '#009E73')
direction_colors <- c(up = '#D55E00', down = '#0072B2', ndiff = '#B3B3B3')
condition_colors <- c(Cont = '#777777', Het = '#777777', DKO = '#222222', KO = '#222222')
signed_colors <- c(direction_colors['down'], neutral = '#F7F7F7', direction_colors['up'])
sequential_palette <- 'cividis'
impact_colors <- c('ORF disrupted upon inclusion' = '#D55E00',
  'ORF disrupted upon exclusion' = '#0072B2', 'Alternative isoform' = '#009E73',
  Regulatory = '#CC79A7', 'Unknown/noncoding' = '#B3B3B3')
overlap_colors <- c('Both: same direction' = '#222222', 'Both: opposite direction' = '#E69F00',
  'Spire only' = dataset_colors[['spire']], 'FMN2 only' = dataset_colors[['fmn2']], Neither = '#D9D9D9')
