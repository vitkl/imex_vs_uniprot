## Function downloads reference proteome list from Uniprot, saves it,
## finds and returns Proteome_ID and SPECIES_ID for given SPECIES_NAME
## default SPECIES_NAME is "Homo sapiens"

SPECIES_NAME_TO_ID <- function(SPECIES_NAME = "Homo sapiens") {
url <-  "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/README"

UniProtSpecies_filename <- paste("./Data/","UniProtSpecies_", Sys.Date(), ".txt", sep = "")
if(!file.exists(UniProtSpecies_filename)) {downloader::download(url, destfile = UniProtSpecies_filename)}

UniProtSpecies_temp = as.character(readLines(UniProtSpecies_filename))
start = grep("Proteome_ID", UniProtSpecies_temp, ignore.case = TRUE)
end = grep("Gene mapping files", UniProtSpecies_temp, ignore.case = TRUE)
UniProtSpecies = UniProtSpecies_temp[start:end]

## type in species name 

SPECIES_NAME = SPECIES_NAME

SPECIES_extracted = grep(SPECIES_NAME, UniProtSpecies, ignore.case = TRUE)
SPECIES_UNIPROT_details = UniProtSpecies[SPECIES_extracted]
print(SPECIES_UNIPROT_details)

# if multiple strains of the same species
if(length(SPECIES_UNIPROT_details)>1){
  SPECIES_UNIPROT_details = SPECIES_UNIPROT_details[1]
}

## code to extract Proteome_ID and SPECIES_ID
Proteome_ID = gsub(" .*$","",SPECIES_UNIPROT_details)
x = gsub(paste(Proteome_ID,""),"",SPECIES_UNIPROT_details)
SPECIES_ID = gsub(" .*$","",x)


ID = data.frame(Proteome_ID = Proteome_ID, SPECIES_ID = SPECIES_ID, stringsAsFactors = F)

print(ID)
return(ID)
}
