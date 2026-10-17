message('Sequence features, protein impact and GO')
annotation <- as.data.frame(fread(annotation_path, na.strings = c('', 'NA')))
stopifnot(all(c('EVENT', 'GENE', 'COORD_o') %in% names(annotation)), !anyDuplicated(annotation$EVENT))
matched <- events$EVENT %in% annotation$EVENT
annotation <- annotation[match(events$EVENT, annotation$EVENT), ]
stopifnot(identical(annotation$EVENT[matched], events$EVENT[matched]))
changed_coords <- matched & !is.na(annotation$COORD_o) & annotation$COORD_o != events$COORD
write_table(data.frame(EVENT = events$EVENT[changed_coords], input_coord = events$COORD[changed_coords],
  reference_coord = annotation$COORD_o[changed_coords]), 'VASTDB_changed_event_coordinates')
annotation[changed_coords, setdiff(names(annotation), c('EVENT', 'GENE'))] <- NA
renamed <- matched & !is.na(annotation$GENE) & annotation$GENE != events$GENE
write_table(data.frame(EVENT = events$EVENT[renamed], input_gene = events$GENE[renamed],
  reference_gene = annotation$GENE[renamed], direction = events$direction[renamed]), 'VASTDB_gene_name_changes')
write_table(events[!matched, c('EVENT', 'GENE', 'COMPLEX', 'direction')], 'missing_VASTDB_event_annotations')
annotation$EVENT <- events$EVENT
annotation$reference_GENE <- annotation$GENE
annotation$GENE <- events$GENE
source_files <- c(source_files, annotation_path, 'data/reference/PROT_IMPACT-mm10-v3.tab.gz')
