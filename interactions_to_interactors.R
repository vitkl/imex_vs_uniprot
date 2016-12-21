## function clears interactors IDs from interaction databases
## returns a data frame which contains:
## interactor_IDs
## interactor_IDs_databases - the database where ID come from (like, uniprotkb)
## interactor_SPECIES_ID - the SPECIES_ID of the interacting protein
## interactor_databases - the IMEX database where interaction information is stored

interactions_to_interactors = function(all_interactions){
temp1 = c(all_interactions$V1,all_interactions$V2)
temp1_1 = c(all_interactions$V10, all_interactions$V11)
temp2 = c(all_interactions$V13,all_interactions$V13)
all_interactors = data.frame(interactor_IDs = temp1, interactor_IDs_databases = temp1, interactor_databases = temp2, interactor_SPECIES_ID = temp1_1, stringsAsFactors = F)
all_interactors$interactor_IDs = gsub("^.*:","",all_interactors$interactor_IDs)
all_interactors$interactor_IDs_databases = gsub(":.*$","",all_interactors$interactor_IDs_databases)
all_interactors$interactor_SPECIES_ID = gsub("^.*:","",all_interactors$interactor_SPECIES_ID)
all_interactors$interactor_SPECIES_ID = gsub("[[:punct:]]{1}([[:print:]]+)[[:punct:]]{1}$","",all_interactors$interactor_SPECIES_ID)
all_interactors$interactor_databases = gsub("^.*MI:","",all_interactors$interactor_databases)
all_interactors = unique(all_interactors)
return(all_interactors)
}