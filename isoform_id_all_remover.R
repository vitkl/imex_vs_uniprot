## function to remove isoform one "-1" XXXXXX-1 from IDs
isoform_id_all_remover <- function(isoform_ids) {
  isoform_id_all_removed <- gsub("-[[:digit:]]+$","",isoform_ids)
}