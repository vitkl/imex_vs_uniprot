## intact_vs_swissprot_vs_missing_swissprot analysis

## human only
SPECIES_NAME = c("Homo sapiens")

## Swiss-Prot only
reviewed = c(2)

## ## Distinguish between isoforms or use only generic Uniprot IDs: TRUE / FALSE?
isoforms = c(FALSE)

source("SPECIES_NAME_TO_ID.R")
SPECIES_IDs = SPECIES_NAME_TO_ID(SPECIES_NAME)

SPECIES_ID = SPECIES_IDs$SPECIES_ID;
Proteome_ID = SPECIES_IDs$Proteome_ID;


## Query Uniprot for the reference proteome for SPECIES_ID
source("download_reference_proteome_q.R")
reference_proteome_query = download_reference_proteome_q(SPECIES_ID, Proteome_ID)

## Use only Swissprot data if reviewed == 2 
if(reviewed == 2) {
  reference_proteome_query = dplyr::filter(reference_proteome_query, Status == "reviewed")
}
## Use only TrEMBL data if reviewed == 3 
if(reviewed == 3) {
  reference_proteome_query = dplyr::filter(reference_proteome_query, Status == "unreviewed")
}

table(reference_proteome_query$Protein.existence)
length(reference_proteome_query$Protein.existence)

## Query Uniprot for the whole proteome for SPECIES_ID
source("download_whole_proteome.R")
whole_proteome_query = download_whole_proteome(SPECIES_ID)

## Use only Swissprot data if reviewed == 2
if(reviewed == 2) {
  whole_proteome_query = dplyr::filter(whole_proteome_query, Status == "reviewed")
}
## Use only TrEMBL data if reviewed == 3
if(reviewed == 3) {
  whole_proteome_query = dplyr::filter(whole_proteome_query, Status == "unreviewed")
}

table(whole_proteome_query$Protein.existence)
length(whole_proteome_query$Protein.existence)

## In case isoform argument to the function is TRUE  - to extract isoform IDs 
## from queried proteomes and merge with generic IDs
if(isoforms == TRUE){
  source("isoform_id_extractor_Uniprot.R")
  reference_proteome_isoforms = isoform_id_extractor_Uniprot(reference_proteome_query)
  reference_proteome = as.character(c(reference_proteome_isoforms, reference_proteome_query$Entry))
  whole_proteome_isoforms = isoform_id_extractor_Uniprot(whole_proteome_query)
  whole_proteome = as.character(c(whole_proteome_isoforms, whole_proteome_query$Entry))
}

## In case isoform argument to the function is FALSE  - to extract generic IDs 
## from queried proteomes
if(isoforms == FALSE){
  reference_proteome = reference_proteome_query$Entry
  whole_proteome = whole_proteome_query$Entry
}

## In case isoform argument to the function is TRUE - REMOVE "-1" XXXXXX-1 from IDs
if(isoforms == TRUE){
  source("isoform_id_1_remover.R")
  all_interactors_SPECIES_ID_only$interactor_IDs = isoform_id_1_remover(all_interactors_SPECIES_ID_only$interactor_IDs)
  reference_proteome = isoform_id_1_remover(reference_proteome)
  whole_proteome = isoform_id_1_remover(whole_proteome)
}

type = character(5)
type[1] = "Evidence_at_transcript_level"
type[2] = "Inferred_from_homology"
type[3] = "Predicted"
type[4] = "Uncertain"
type[5] = "Evidence_at_protein_level"

url = character(5)
url[1] = "https://api.nextprot.org/export/entries.txt?sparql=select%20distinct%20%3Fentry%20where%20%7B%0A%3Fentry%20%3Aexistence%20%3AEvidence_at_transcript_level.%0A%7D"
url[2] = "https://api.nextprot.org/export/entries.txt?sparql=select%20distinct%20%3Fentry%20where%20%7B%0A%3Fentry%20%3Aexistence%20%3AInferred_from_homology.%0A%7D%0A"
url[3] = "https://api.nextprot.org/export/entries.txt?sparql=select%20distinct%20%3Fentry%20where%20%7B%0A%3Fentry%20%3Aexistence%20%3APredicted.%0A%7D%0A"
url[4] = "https://api.nextprot.org/export/entries.txt?sparql=select%20distinct%20%3Fentry%20where%20%7B%0A%3Fentry%20%3Aexistence%20%3AUncertain.%0A%7D%0A"
url[5] = "https://api.nextprot.org/export/entries.txt?sparql=select%20distinct%20%3Fentry%20where%20%7B%0A%3Fentry%20%3Aexistence%20%3AEvidence_at_protein_level.%0A%7D%0A"

missing_swissprot_filename = character(5)
for(i in 1:4){
  missing_swissprot_filename[i] <- paste("./Data/","missing_swissprot_speciesID","_",SPECIES_ID,"_",type[i], "2016-11-28.txt", sep = "")
}

Evidence_at_transcript_level = data.frame(read.delim(missing_swissprot_filename[1],stringsAsFactors = F, quote = ""), evidence = type[1])
colnames(Evidence_at_transcript_level)[1] = "entry"
Inferred_from_homology = data.frame(read.delim(missing_swissprot_filename[2],stringsAsFactors = F, quote = ""), evidence = type[2])
colnames(Inferred_from_homology)[1] = "entry"
Predicted = data.frame(read.delim(missing_swissprot_filename[3],stringsAsFactors = F, quote = ""), evidence = type[3])
colnames(Predicted)[1] = "entry"
Uncertain = data.frame(read.delim(missing_swissprot_filename[4],stringsAsFactors = F, quote = ""), evidence = type[4])
colnames(Uncertain)[1] = "entry"
Evidence_at_protein_level = data.frame(read.delim(missing_swissprot_filename[5],stringsAsFactors = F, quote = ""), evidence = type[5])
colnames(Evidence_at_protein_level)[1] = "entry"

missing_swissprot = rbind(Evidence_at_transcript_level, Inferred_from_homology, Predicted, Uncertain)
all_swissprot = rbind(missing_swissprot, Evidence_at_protein_level)
table(all_swissprot$evidence)
length(all_swissprot$evidence)
missing_swissprot$entry = gsub("NX_", "", missing_swissprot$entry)
all_swissprot$entry = gsub("NX_", "", all_swissprot$entry)

missing_swissprot2 = data.frame(missing_swissprot, missing_swissprot = 1, stringsAsFactors = F)
all_swissprot2 = data.frame(all_swissprot, all_swissprot_nextprot = 1, stringsAsFactors = F)
swissprot_from_uniprot = data.frame(whole_proteome_query$Entry, whole_proteome_query$Protein.existence, swissprot_from_uniprot = 1, stringsAsFactors = F)

x = merge(all_swissprot2, swissprot_from_uniprot, all = TRUE, by.x = "entry", by.y = "whole_proteome_query.Entry" )
x2=x
x2$all_swissprot_nextprot[is.na(x2$all_swissprot_nextprot)] <- 0
x2$swissprot_from_uniprot[is.na(x2$swissprot_from_uniprot)] <- 0

sum((x2$all_swissprot_nextprot == 1) & (x2$swissprot_from_uniprot == 1))
dim(x2)

table(x2$evidence)
table(x2$whole_proteome_query.Protein.existence)
x2$whole_proteome_query.Protein.existence = gsub(" ", "_", x2$whole_proteome_query.Protein.existence)

sum((x2$evidence != x2$whole_proteome_query.Protein.existence) & (x2$evidence == "Evidence_at_protein_level"), na.rm = T)
