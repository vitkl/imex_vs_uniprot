## the function to extract isoforms given
## reference_proteome or whole_proteome

## Regular expression of uniprot isoforms ID is ([[:alnum:]]+-[[:digit:]]+)

variant_id_extractor_Uniprot <- function(proteome) {
  proteome$Cross.reference..dbSNP.
variants <- as.character(proteome$Cross.reference..dbSNP.)
variant_ids <- unlist(gsubfn::strapply(X = variants,pattern = "rs[[:digit:]]+ ",FUN = length, simplify = T))
variant_number = sum(variant_ids)

return(variant_number)

}