## the function to download the reference proteome 
## (isoforms in Alternative.products..isoforms.) from Uniprot given:
### SPECIES_ID (default is ID for human)
### Proteome_ID (default is ID for human)
## the function saves the query result, reads it and returns a data.frame 
## in which the 1st column contains Uniprot identifiers

## in contrast to download_reference_proteome does a query instead 
## of downloading file from ftp

download_reference_proteome_q <- function(SPECIES_ID = "9606", Proteome_ID = "UP000005640", date = Sys.Date()) {
  
  print("Querying Uniprot for the reference proteome")
  
  url = paste("http://www.uniprot.org/uniprot/?sort=score&desc=&compress=no&query=proteome:",Proteome_ID, "&fil=&force=no&preview=true&format=tab&columns=id,reviewed,length,organism,organism-id,mass,database(dbSNP),comment(ALTERNATIVE%20PRODUCTS),annotation%20score,existence,protein%20names",sep = "")
  
  UniProtLIST_filename <- paste("./Data/","Reference_proteome_speciesID","_",SPECIES_ID,"_",Proteome_ID,"_", date, sep = "")
  
  ## If the proteome has not been downloaded -> 
  ## download it from url and save into UniProtLIST_filename
  if(!file.exists(UniProtLIST_filename)) { 
    downloader::download(url, destfile = UniProtLIST_filename) 
    print("Querying Uniprot for the reference proteome...loaded from uniprot")}
  if(file.size(UniProtLIST_filename) == 0) {}
    
  ## read the reference proteome, some columns may include quotes
  UniProtLIST = as.data.frame(read.delim(UniProtLIST_filename,stringsAsFactors = F, quote = ""))
  
  print(paste("Number of proteins (not including isoforms, but isoforms dowloaded, from query) in the reference proteome for SPECIES_ID", SPECIES_ID,"total/unique"))
  print(table(UniProtLIST$Organism.ID))
  print(length(unique(UniProtLIST$Entry)))
  
  print("Querying Uniprot for the reference proteome...done")
  
  ## return the reference proteome
  return(UniProtLIST)
}