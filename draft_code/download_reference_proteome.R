## the function to download the reference proteome (isoforms included) from Uniprot given:
### Proteome_ID (default is ID for human)
### SPECIES_ID (default is ID for human)
### domain (default is "Eukaryota")
## the function saves, reads and returns the gene2acc data.frame 
## in which 2nd column contains Uniprot identifiers

## in contrast to download_reference_proteome_q downloads file from ftp instead 
## of querying the database and includes isoforms

download_reference_proteome <- function(SPECIES_ID = "9606", Proteome_ID = "UP000005640", domain = "Eukaryota") {

url <-  "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/"
url2 <- paste(url,domain,"/", Proteome_ID,"_",SPECIES_ID,".gene2acc.gz", sep = "")
UniProtLIST_filename <- paste("./Data/","Reference_proteome_speciesID","_",SPECIES_ID,"_", Sys.Date(),".gene2acc.gz", sep = "")
UniProtLIST_read <- paste("./Data/","Reference_proteome_speciesID","_",SPECIES_ID,"_", Sys.Date(),".gene2acc", sep = "")

if(!file.exists(UniProtLIST_read)) {
downloader::download(url2, destfile = UniProtLIST_filename)
R.utils::gunzip(UniProtLIST_filename)
}

## read the reference proteome
UniProtLIST = as.data.frame(read.table(UniProtLIST_read))
print(paste("Number of proteins (isoforms included, from ftp) in the reference proteome for SPECIES_ID", SPECIES_ID))
print(dim(UniProtLIST)[1])
print(dim(UniProtLIST)[1])

## return the reference proteome
return(UniProtLIST)
}


