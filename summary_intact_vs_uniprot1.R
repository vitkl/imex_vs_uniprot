## Summary function - provides summary for:
# whole_proteome_unique
# reference_proteome_unique
# reference_proteome_ftp
# all_interactors_SPECIES_ID_only

summary_intact_vs_uniprot1 = function(whole_proteome_unique, reference_proteome_unique, reference_proteome_ftp, all_interactors_SPECIES_ID_only) {

  summary_intact_vs_uniprot = list(
      reference_proteome = length(reference_proteome_unique$Reference_proteome_IDs),
      #reference_proteome_from_ftp = length(unique(reference_proteome_ftp$V2)),
      whole_proteome = length(whole_proteome_unique$whole_proteome_IDs),
      interactome = dim(unique(all_interactors_SPECIES_ID_only))[1],
      interactome_unique = length(unique(all_interactors_SPECIES_ID_only$interactor_IDs))
      # interactome_db_overlap = length(all_interactors_SPECIES_ID_only$interactor_IDs) - length(unique(all_interactors_SPECIES_ID_only$interactor_IDs)) 
      )

print(paste("number of proteins (unique, including isoforms) in the reference proteome is: ", summary_intact_vs_uniprot$reference_proteome))
#print(paste("number of proteins (unique, including isoforms, from ftp) in the reference proteome is: ", summary_intact_vs_uniprot$reference_proteome_from_ftp))
print(paste("number of proteins (unique, including isoforms) in the total (whole) proteome is: ", summary_intact_vs_uniprot$whole_proteome))
print(paste("number of proteins (including isoforms, incl. overlap between databases) included in the interactome by IMEX: ", summary_intact_vs_uniprot$interactome))
print(paste("number of proteins (unique, including isoforms, isoform -1 not counted separately) included in the interactome by IMEX: ", summary_intact_vs_uniprot$interactome_unique))
# print(paste("overlap between interaction databases(IMEX), not necessarily unique proteins: ", summary_intact_vs_uniprot$interactome_db_overlap))

summary_intact_vs_uniprot2 = unlist(summary_intact_vs_uniprot)

return(summary_intact_vs_uniprot2)

}
