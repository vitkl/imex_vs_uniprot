## The function which does overlap calculations

intact_vs_uniprot_overlap = function(SPECIES_NAME, SPECIES_ID, reviewed, isoforms, date = Sys.Date()){
  
## read saved logic table
filename_vs_2 = paste("./analysis/","proteome_vs_interactome_f_", SPECIES_ID,"_reviewed_",reviewed,"_isoforms_",isoforms,"_", date,".txt", sep = "")

if(!file.exists(filename_vs_2)){
  message(paste(filename_vs_2, "doesn't exist"))
  return(0)
}

  proteome_vs_interactome_o = as.data.frame(read.delim(filename_vs_2, stringsAsFactors = F, quote = ""))
  proteome_vs_interactome_f = proteome_vs_interactome_o

databases_MI = colnames(proteome_vs_interactome_f)
if(databases_MI=="x") {
  message(paste(filename_vs_2, "is empty"))
  return(0)
}

print(str(proteome_vs_interactome_f))

## Summary of the overlap between the interactome (from all IMEx or IntAct only)
## and reference and whole proteomes
source("proteome_vs_interactome_summary.R")
IMEx = proteome_vs_interactome_summary(proteome_vs_interactome_f, database = "IMEx", SPECIES_NAME, reviewed, isoforms) 
#IntAct = proteome_vs_interactome_summary(proteome_vs_interactome_f, database = "X0469.IntAct.", SPECIES_NAME, reviewed, isoforms) 
proteome_vs_interactome_summary = IMEx#rbind(IMEx, IntAct)


source("A_vs_B_vs_C_overlap.R")
whole_proteome_Uniprot_vs_IMEx_vs_protein.exist = A_vs_B_vs_C_overlap(proteome_vs_interactome_f, 
                                                                                 A = "whole_proteome_Uniprot", 
                                                                                 B = "IMEx", 
                                                                                 C = "missing_protein_evidence", 
                                                                                 SPECIES_NAME = n, reviewed = r, isoforms = i)
#whole_proteome_Uniprot_vs_IntAct_vs_protein.exist = A_vs_B_vs_C_overlap(proteome_vs_interactome_f, 
#                                                                      A = "whole_proteome_Uniprot", 
#                                                                      B = "X0469.IntAct.", 
#                                                                      C = "missing_protein_evidence", 
#                                                                      SPECIES_NAME = n, reviewed = r, isoforms = i)


filename.summary = paste("./summaries/","proteome_vs_interactome_summary_",SPECIES_NAME,"_reviewed_",reviewed,"_isoforms_",isoforms,"_",date,".txt", sep = "")
write.table(proteome_vs_interactome_summary, filename.summary, col.names=T,row.names=F,sep="\t",quote=F)

filename.summary2 = paste("./summaries/","whole_proteome_Uniprot_vs_IMEx_vs_protein.exist_summary_",SPECIES_NAME,"_reviewed_",reviewed,"_isoforms_",isoforms,"_",date,".txt", sep = "")
write.table(whole_proteome_Uniprot_vs_IMEx_vs_protein.exist, filename.summary2, col.names=T,row.names=F,sep="\t",quote=F)

#filename.summary3 = paste("./summaries/","whole_proteome_Uniprot_vs_IntAct_vs_protein.exist_summary_",SPECIES_NAME,"_reviewed_",reviewed,"_isoforms_",isoforms,"_",date,".txt", sep = "")
#write.table(whole_proteome_Uniprot_vs_IntAct_vs_protein.exist, filename.summary3, col.names=T,row.names=F,sep="\t",quote=F)

print(proteome_vs_interactome_summary)

return(SPECIES_NAME)
}