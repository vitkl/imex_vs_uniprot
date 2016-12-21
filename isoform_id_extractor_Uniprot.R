## the function to extract isoforms given
## reference_proteome or whole_proteome

## Regular expression of uniprot isoforms ID is ([[:alnum:]]+-[[:digit:]]+)

isoform_id_extractor_Uniprot <- function(proteome) {

isoforms <- as.character(proteome$Alternative.products..isoforms.)
isoform_ids <- unlist(gsubfn::strapplyc(isoforms,"IsoId=([[:alnum:]]+-[[:digit:]]+);",simplify = T))

return(isoform_ids)

}