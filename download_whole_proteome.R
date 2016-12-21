## the function to download the whole proteome 
## (isoforms in Alternative.products..isoforms.) from Uniprot given:
### SPECIES_ID (default is ID for human)
## the function saves the query result, reads it and returns a data.frame 
## in which the 1st column contains Uniprot identifiers

download_whole_proteome <- function(SPECIES_ID = "9606", date = Sys.Date()) {
  
  print("Querying Uniprot for the whole proteome")
  
  url = paste("http://www.uniprot.org/uniprot/?sort=score&desc=&compress=no&query=organism:",SPECIES_ID, "&fil=&force=no&preview=true&format=tab&columns=id,reviewed,length,organism,organism-id,mass,database(dbSNP),comment(ALTERNATIVE%20PRODUCTS),annotation%20score,existence,protein%20names",sep = "")
  # organism:ID - not taxonomy:ID - because the last one returns denisovan 
  # and neanderthal proteins in addition to modern human
  
  UniProtLIST_filename <- paste("./Data/","Whole_proteome_speciesID","_",SPECIES_ID,"_", date, sep = "")
  
  ## If the proteome has not been downloaded -> 
  ## download it from url and save into UniProtLIST_filename
  if(!file.exists(UniProtLIST_filename)) { 
    downloader::download(url, destfile = UniProtLIST_filename) 
    print("Querying Uniprot for the whole proteome...loaded from uniprot")
          }
  
  ## read the whole proteome, escape unexpected quotes
  UniProtLIST = as.data.frame(read.delim(UniProtLIST_filename,stringsAsFactors = F, quote = ""))
  
  print(paste("Number of proteins (not including isoforms, but isoforms dowloaded, from query) in the whole proteome for SPECIES_ID", SPECIES_ID,"total/unique"))
  print(table(UniProtLIST$Organism.ID))
  print(length(unique(UniProtLIST$Entry)))
  
  print("Querying Uniprot for the whole proteome...done")
  
  ## return the whole proteome
  return(UniProtLIST)
}