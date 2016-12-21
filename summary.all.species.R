summary.all.species = function(SPECIES_NAME) {
  filename.summary = paste("proteome_vs_interactome_summary_",SPECIES_NAME,"_",Sys.Date(),".txt", sep = "")
  proteome_vs_interactome_summary_o = as.data.frame(read.delim(filename_vs_2, stringsAsFactors = F))
  
  proteome_vs_interactome_summary.all = rbind(proteome_vs_interactome_summary, proteome_vs_interactome_summary.all)
  
  return(proteome_vs_interactome_summary.all)
}