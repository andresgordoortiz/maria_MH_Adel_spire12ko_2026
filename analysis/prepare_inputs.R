# Validate a completed simulation and prepare the sequence-tool inputs.
source('analysis/load.R')
annotation_path <- 'data/reference/event_annotations.tsv'
source('analysis/annotation_input.R')
status <- data.frame(step = character(), status = character())
prepare_only <- TRUE
source('analysis/external.R')
write_table(status, 'analysis_status')
writeLines('Validated inputs for Matt and rMAPS2; figure analysis is a separate step.', file.path(out, 'COMPLETE'))
