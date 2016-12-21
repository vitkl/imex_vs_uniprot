## ## unique interactors vs proteome function

unique_vs_proteome = function(all_interactors_unique, proteome_unique) {

y1 = unique(all_interactors_unique$interactor_IDs)
y2 = rep("IMEX", length(unique(all_interactors_unique$interactor_IDs)))
interactome = data.frame(interactor_IDs = y1, interactor_databases = y2, stringsAsFactors = F)
y1 = unique(proteome_unique[1])
y2 = rep("Uniprot", length(unique(proteome_unique[1])))
proteome = data.frame(y1, y2, stringsAsFactors = F)
colnames(proteome) = c("protein_IDs","Uniprot")
proteome_vs_interactome = merge(x = proteome, 
                                y = interactome, by.x="protein_IDs", 
                                by.y="interactor_IDs", 
                                all.x = T, all.y = F)
proteome_vs_interactome_summary = as.data.frame(table(proteome_vs_interactome$interactor_databases,useNA="ifany"))
# proteome_vs_interactome_summary2 = as.data.frame(table(proteome_vs_interactome$Uniprot,useNA="ifany"))
colnames(proteome_vs_interactome_summary) = c("database", "number of proteins")
# colnames(proteome_vs_interactome_summary2) = c("database", "number of proteins")
print(proteome_vs_interactome_summary)
# print(proteome_vs_interactome_summary2)
return(proteome_vs_interactome_summary)
}