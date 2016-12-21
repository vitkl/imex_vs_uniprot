## Function which compares proteome with interactome from database
## Generates 

## Default database is all of IMEx

proteome_vs_interactome_summary =  function(proteome_vs_interactome_o, database = "IMEx", SPECIES_NAME, reviewed, isoforms) {

proteome_vs_interactome_summary = data.frame(
  "species name" = SPECIES_NAME,
  "database" = database,
  "reviewed" = reviewed,
  "isoforms" = isoforms,
  "whole proteome (Uniprot)" = sum((proteome_vs_interactome_o$whole_proteome_Uniprot == 1)),
  "uncovered whole proteome" = sum((proteome_vs_interactome_o$whole_proteome_Uniprot == 1)) - sum((proteome_vs_interactome_o$whole_proteome_Uniprot == 1) & (1 == proteome_vs_interactome_o[database])),
  "whole proteome, interactome available" = sum((proteome_vs_interactome_o$whole_proteome_Uniprot == 1) & (1 == proteome_vs_interactome_o[database])),
  "interactome, but not in Uniprot" = sum((proteome_vs_interactome_o$whole_proteome_Uniprot == 0) & (1 == proteome_vs_interactome_o[database])),
  "reference proteome (Uniprot)" = sum((proteome_vs_interactome_o$reference_proteome_Uniprot == 1)),
  "uncovered reference proteome" = sum((proteome_vs_interactome_o$reference_proteome_Uniprot == 1)) - sum((proteome_vs_interactome_o$reference_proteome_Uniprot == 1) & (1 == proteome_vs_interactome_o[database])),
  "reference proteome, interactome available" = sum((proteome_vs_interactome_o$reference_proteome_Uniprot == 1) & (1 == proteome_vs_interactome_o[database])),
  "interactome, but not in Uniprot" = sum((proteome_vs_interactome_o$reference_proteome_Uniprot == 0) & (1 == proteome_vs_interactome_o[database])),
  
  "Overlap between interactome and the reference proteome, fraction" = mean((proteome_vs_interactome_o$reference_proteome_Uniprot == 1) & (1 == proteome_vs_interactome_o[database]))/mean(proteome_vs_interactome_o$reference_proteome_Uniprot == 1),
  "Overlap between interactome and the whole proteome, fraction" = mean((proteome_vs_interactome_o$whole_proteome_Uniprot == 1) & (1 == proteome_vs_interactome_o[database]))
)
proteome_vs_interactome_summary2 = (proteome_vs_interactome_summary)
return(proteome_vs_interactome_summary2)

}