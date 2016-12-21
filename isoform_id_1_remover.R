## function to remove isoform one "-1" XXXXXX-1 from IDs
isoform_id_1_remover <- function(isoform_ids_1) {
isoform_ids_1_removed <- gsub("-1$","",isoform_ids_1)
}