## the function to query PSICQUIC for interactions given:
### SPECIES_NAME
### SPECIES_ID (default NA)
### databases
### date - option allows to read already saved files
### detmethod - optional
### pmethod - optional

query_PSICQUIC_for_interactions = function(SPECIES_ID = NA, SPECIES_NAME = NA, databases = NA, date = Sys.Date(), detmethod = NA, pmethod = NA) {

## checks if databases argument was provided, if not - sets default - all IMEX databases
if(is.na(databases)[1]){ 
  databases = c("IntAct", "MINT", "bhf-ucl", "MPIDB", "MatrixDB", "HPIDb","I2D-IMEx","InnateDB-IMEx", "MolCon", "UniProt", "MBInfo")
}

## Converts SPECIES_NAME to SPECIES_ID if SPECIES_ID is not stated
if(is.na(SPECIES_ID) & !is.na(SPECIES_NAME)) { 
  source("SPECIES_NAME_TO_ID.R")
  SPECIES_ID = SPECIES_NAME_TO_ID(SPECIES_NAME)$SPECIES_ID }
else {SPECIES_ID = SPECIES_ID}
##============================================================================##
  ## Constructs database filename and PSICQUIC_query depending on 
  ## whether detmethod and pmethod were provided as arguments
if(!is.na(detmethod)){
  if(!is.na(pmethod)){
    database.name <- paste("./Data/","databaseName_", databases[1], "...", databases[length(databases)],"_","speciesID_",SPECIES_ID,"_",SPECIES_NAME,"_detmethod_",detmethod, "_pmethod_", pmethod, "_", date, sep = "")
    PSICQUIC_query = paste("species:",SPECIES_ID," AND ","detmethod:",detmethod," AND ","pmethod:",pmethod, sep = "")
  }
  if(is.na(pmethod)){
    database.name <- paste("./Data/","databaseName_", databases[1], "...", databases[length(databases)],"_","speciesID_",SPECIES_ID,"_",SPECIES_NAME,"_detmethod_",detmethod, "_", date, sep = "")
    PSICQUIC_query = paste("species:",SPECIES_ID," AND ","detmethod:",detmethod, sep = "")
  }
}
if(is.na(detmethod)){
    database.name <- paste("./Data/","databaseName_", databases[1], "...", databases[length(databases)],"_","speciesID_",SPECIES_ID,"_",SPECIES_NAME,"_",date, sep = "")
    PSICQUIC_query = paste("species:",SPECIES_ID, sep = "")
}
## Checks if databases have been queried today, if not - sends query to the database
if(!file.exists(database.name)) {
  print("dowloaded using PSICQUIC")
## Load PSICQUIC functionality
library(PSICQUIC)
psicquic <- PSICQUIC()
providers <- providers(psicquic)

# query databases for all known SPECIES_ID protein interactions
SPECIES_ID_interactome = data.frame()
NO_SPECIES_ID_interactome = character(length = length(databases))
for(indices in 1:length(databases)) {
  if(databases[indices] %in% providers) {
    ## Query
    SPECIES_ID_interactome_d <- rawQuery(psicquic, databases[indices], PSICQUIC_query)
    ## add query results to the single table if the table format is compatible
    if(ncol(SPECIES_ID_interactome_d)==ncol(SPECIES_ID_interactome) | (indices==1)) {
      SPECIES_ID_interactome <- rbind(SPECIES_ID_interactome, SPECIES_ID_interactome_d)
    }
    ## if the query finds no proteins or the table format is not compatible
    ## save that information
    else {
      NO_SPECIES_ID_interactome[indices] = databases[indices]
    }
  }
}
##============================================================================##
## Save dowloaded query result into file 
save(SPECIES_ID_interactome, file = database.name)
##============================================================================##
## Show what's found
print(paste("interactions for ",SPECIES_NAME,", detmethod(",detmethod,"), pmethod(",pmethod,"): ", sep = ""), quote = F)
print(paste0("total number: ", nrow(SPECIES_ID_interactome)), quote = F)
print(paste("there is no interactions in the databases: ", sep = ""), quote = F)
print(NO_SPECIES_ID_interactome[,1], quote = F)
print(paste("the number of interactions per database ", sep = ""), quote = F)
dbs = as.data.frame(table(SPECIES_ID_interactome$V13, useNA = "ifany"))
colnames(dbs) = c("database", "N of interactions")
print(dbs, quote = F)
##============================================================================##
query_log_filename = paste("./Data/logs/","there is no interactions for ",SPECIES_NAME," in the databases ",date, sep = "", ".txt")
write.table(NO_SPECIES_ID_interactome, query_log_filename, col.names=T,row.names=F,sep="\t",quote=F)
##============================================================================##
return(SPECIES_ID_interactome)
}
##============================================================================##
## If file exists  - load, show what's found and return
if(file.exists(database.name)) {
  print("loaded from file")
load(database.name)
##============================================================================##  
  query_log_filename = paste("./Data/logs/","there is no interactions for ",SPECIES_NAME," in the databases ",date, sep = "", ".txt")
  NO_SPECIES_ID_interactome = read.table(file = query_log_filename, header =T,sep="\t", stringsAsFactors = F)

  ## Show what's found
  print(paste("interactions for ",SPECIES_NAME,", detmethod(",detmethod,"), pmethod(",pmethod,"): ", sep = ""), quote = F)
  print(paste0("total number: ", nrow(SPECIES_ID_interactome)), quote = F)
  print(paste("there is no interactions in the databases: ", sep = ""), quote = F)
  print(NO_SPECIES_ID_interactome[,1], quote = F)
  print(paste("the number of interactions per database ", sep = ""), quote = F)
  dbs = as.data.frame(table(SPECIES_ID_interactome$V13, useNA = "ifany"))
  colnames(dbs) = c("database", "N of interactions")
  print(dbs, quote = F)
##============================================================================##
return(SPECIES_ID_interactome)
}
}