### The function generates a logic table given:
## reference and whole proteomes and interactome
### logic table has a 1 when protein(by rows) is present in a database (by columns)
### logic table has a 0 when protein(by rows) is absent in a database (by columns)
## logic tables for reference proteome and whole proteome are saved
## logic tables for whole proteome (contains values of reference proteome) is returned

proteome_vs_interactome_logic_table_generator = function(reference_proteome_unique, whole_proteome_unique, all_interactors_SPECIES_ID_only, unique_interactors_SPECIES_ID_only) {

if(nrow(reference_proteome_unique)==0 | nrow(whole_proteome_unique)==0) {
  message(paste("no reference or whole proteome"))
  return (0)
  }
  
  # add columns of ones 
reference_proteome_unique = as.data.frame(cbind(reference_proteome_unique, rep(1, length(reference_proteome_unique))))
colnames(reference_proteome_unique)[2] = c("reference_proteome_Uniprot")
whole_proteome_unique = as.data.frame(cbind(whole_proteome_unique, rep(1, length(whole_proteome_unique))))
colnames(whole_proteome_unique)[2] = c("whole_proteome_Uniprot")
all_interactors_SPECIES_ID_only = as.data.frame(cbind(all_interactors_SPECIES_ID_only, rep(1, dim(all_interactors_SPECIES_ID_only)[1])))
IMEx = data.frame(interactor_IDs = unique_interactors_SPECIES_ID_only$interactor_IDs, rep(1, length(unique_interactors_SPECIES_ID_only$interactor_IDs)))
colnames(IMEx)[length(IMEx)] = "IMEx"

all_interactors_SPECIES_ID_only.trans = tidyr::spread(data = all_interactors_SPECIES_ID_only, key = interactor_databases, value = `rep(1, dim(all_interactors_SPECIES_ID_only)[1])`, fill = 0) 
all_interactors_SPECIES_ID_only.trans.imex = merge(x = all_interactors_SPECIES_ID_only.trans, 
                                          y = IMEx, by.x="interactor_IDs", 
                                          by.y="interactor_IDs", 
                                          all.x = T, all.y = T)
unique_interactors_SPECIES_ID_only.trans.imex = unique(all_interactors_SPECIES_ID_only.trans.imex)


## Add whole reference_proteome column of ones
reference_proteome_vs_interactome = merge(x = reference_proteome_unique, 
                                          y = unique_interactors_SPECIES_ID_only.trans.imex, by.x="Reference_proteome_IDs", 
                                          by.y="interactor_IDs", 
                                          all.x = T, all.y = T)
## Add whole whole_proteome column of ones
proteome_vs_interactome = merge(x = whole_proteome_unique, 
                                y = reference_proteome_vs_interactome, by.x="whole_proteome_IDs", 
                                by.y="Reference_proteome_IDs", 
                                all.x = T, all.y = T)

## substitute NAs for zeros in whole proteome-containing data.frame
proteome_vs_interactome_f =proteome_vs_interactome
proteome_vs_interactome_f[is.na(proteome_vs_interactome_f)] <- 0

return(proteome_vs_interactome_f)
}
