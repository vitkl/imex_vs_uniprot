## Function to filter out NOT uniprotkb indentifiers and NOT SPECIES_ID proteins
## Function prints a summary of how much proteins were kept/discarded

uniprotkb_and_SPECIES_ID_interactor_summary = function(all_interactors, SPECIES_ID, SPECIES_NAME) {

library(dplyr)
## filter for only uniprotkb indentifiers 
all_interactors_uniprotkb_only = filter(all_interactors, interactor_IDs_databases == "uniprotkb")
all_interactors_other_only = filter(all_interactors, interactor_IDs_databases != "uniprotkb")
## filter for only SPECIES_ID proteins
all_interactors_SPECIES_ID_only = filter(all_interactors_uniprotkb_only, interactor_SPECIES_ID == SPECIES_ID)
all_interactors_not_SPECIES_ID_only = filter(all_interactors_uniprotkb_only, interactor_SPECIES_ID != SPECIES_ID)

summary_intact_vs_uniprot = data.frame(
  "all interactors" = length(unique(all_interactors$interactor_IDs)),
  "interactors with the UniprotKB identifier" = length(unique(all_interactors_uniprotkb_only$interactor_IDs)),
  "interactors with the other identifier" = length(unique(all_interactors_other_only$interactor_IDs)),
  "SPECIES_ID" = SPECIES_ID,
  "SPECIES_NAME" = SPECIES_NAME,
  "SPECIES_ID interactors" = length(unique(all_interactors_SPECIES_ID_only$interactor_IDs)),
  "interactors from other species" = length(unique(all_interactors_not_SPECIES_ID_only$interactor_IDs))
)
summary_intact_vs_uniprot2 = (summary_intact_vs_uniprot)
return(summary_intact_vs_uniprot2)
}