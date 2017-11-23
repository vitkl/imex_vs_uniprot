biogrid_from_mentha_vs_proteome_vs_imex = function(SPECIES_NAME, reviewed, isoforms, date = Sys.Date()){

source("SPECIES_NAME_TO_ID.R")
SPECIES_IDs = SPECIES_NAME_TO_ID(SPECIES_NAME)

SPECIES_ID = SPECIES_IDs$SPECIES_ID;
Proteome_ID = SPECIES_IDs$Proteome_ID;

databases <- c("mentha")
## Query PSICQUIC for interactions, get MI-TAB-2.5, save, return
source("query_PSICQUIC_for_interactions.R")
mentha = query_PSICQUIC_for_interactions(SPECIES_ID = SPECIES_ID, SPECIES_NAME = SPECIES_NAME, databases=databases, date)

library(dplyr)
biogrid_from_mentha = filter(mentha, V13 == "psi-mi:MI:0463(biogrid)")

## the function extracts interactor IDs from interaction databases
source("interactions_to_interactors.R")
biogrid_from_mentha_all = interactions_to_interactors(biogrid_from_mentha)

## filter interactor for uniprotkb only indentifiers 
## filter for SPECIES_ID only proteins
source("uniprotkb_and_SPECIES_ID_interactor_selector.R")
biogrid_from_mentha_SPECIES_ID_only = uniprotkb_and_SPECIES_ID_interactor_selector(all_interactors=biogrid_from_mentha_all, SPECIES_ID = SPECIES_ID)

#============================================================================#
## In case isoform argument to the function is TRUE - REMOVE "-1" XXXXXX-1 from IDs
if(isoforms == TRUE){
  source("isoform_id_1_remover.R")
  biogrid_from_mentha_SPECIES_ID_only$interactor_IDs = isoform_id_1_remover(biogrid_from_mentha_SPECIES_ID_only$interactor_IDs)
}

## In case isoform argument to the function is FALSE - REMOVE All isoform IDs (XXXXXX-X+) from IDs
if(isoforms == FALSE){
  source("isoform_id_all_remover.R")
  biogrid_from_mentha_SPECIES_ID_only$interactor_IDs = isoform_id_all_remover(biogrid_from_mentha_SPECIES_ID_only$interactor_IDs)
}

#============================================================================#
biogrid_from_mentha_SPECIES_ID_only_unique = unique(biogrid_from_mentha_SPECIES_ID_only[c("interactor_IDs","interactor_IDs_databases", "interactor_SPECIES_ID")])
biogrid_from_mentha_temp = data.frame(interactor_IDs = biogrid_from_mentha_SPECIES_ID_only_unique$interactor_IDs, rep(1, length(biogrid_from_mentha_SPECIES_ID_only_unique$interactor_IDs)),stringsAsFactors = FALSE)
colnames(biogrid_from_mentha_temp)[length(biogrid_from_mentha_temp)] = "BioGRID_from_Mentha"

## Summary of how many interactors have non-Uniprot identifiers 
## and how many interactors come from the species other than queried
source("uniprotkb_and_SPECIES_ID_interactor_summary.R")
uniprotkb_and_SPECIES_ID_interactor_summary = uniprotkb_and_SPECIES_ID_interactor_summary(biogrid_from_mentha_all, SPECIES_ID, SPECIES_NAME)
filename = paste("./summaries/","uniprotKB_IDs_and_",SPECIES_NAME,"_biogrid_from_mentha_interactors_summary", "_isoforms_",isoforms,"_", date,".txt", sep = "")
write.table(uniprotkb_and_SPECIES_ID_interactor_summary,filename,col.names=T,row.names=F,sep="\t",quote=F)
print(uniprotkb_and_SPECIES_ID_interactor_summary)

#============================================================================#
## Combine the "logic table" proteome_vs_interactome_f with BioGRIG
filename_vs_2 = paste("./analysis/","proteome_vs_interactome_f_", SPECIES_ID,"_reviewed_",reviewed,"_isoforms_",isoforms,"_", date,".txt", sep = "")
proteome_vs_interactome_o = as.data.frame(read.delim(filename_vs_2, stringsAsFactors = F, quote = ""))

biogrid_from_mentha_vs_proteome_vs_imex = merge(x = proteome_vs_interactome_o, 
                                                y = biogrid_from_mentha_temp, by.x="whole_proteome_IDs", 
                                                by.y="interactor_IDs", 
                                                all.x = T, all.y = T)
biogrid_from_mentha_vs_proteome_vs_imex = unique(biogrid_from_mentha_vs_proteome_vs_imex)
## substitute NAs for zeros in whole proteome-containing data.frame
biogrid_from_mentha_vs_proteome_vs_imex_f = biogrid_from_mentha_vs_proteome_vs_imex
biogrid_from_mentha_vs_proteome_vs_imex_f[is.na(biogrid_from_mentha_vs_proteome_vs_imex_f)] <- 0

filename_vs_3 = paste("./analysis/","proteome_vs_interactome_vs_BioGRID_f_", SPECIES_ID,"_reviewed_",reviewed,"_isoforms_",isoforms,"_", date,".txt", sep = "")
write.table(biogrid_from_mentha_vs_proteome_vs_imex_f,filename_vs_3,col.names=T,row.names=F,sep="\t",quote=F)

}
